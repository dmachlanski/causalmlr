# Evaluation metrics for causal estimators.

#' Absolute error of an ATE estimate
#'
#' Measures the absolute difference between an estimated and the true average
#' treatment effect:
#' \deqn{\epsilon_{ATE} = | \widehat{ATE} - ATE |}
#' Requires the true ATE, so it is only usable with (semi-)synthetic
#' benchmark data such as [sodium] or [ihdp_train].
#'
#' @param ate_true True ATE (a single number).
#' @param ate_pred Estimated ATE. Either a number or a `"causalmlr_ate"`
#'   object as returned by [ate_naive()] and friends.
#'
#' @return A single non-negative number.
#'
#' @examples
#' data(sodium)
#' est <- ate_naive(sodium, outcome = "bp", treatment = "sodium")
#' eps_ate(1.05, est)
#'
#' @seealso [pehe()], [r_loss()]
#' @export
eps_ate <- function(ate_true, ate_pred) {
  if (inherits(ate_true, "causalmlr_ate")) ate_true <- ate_true$estimate
  if (inherits(ate_pred, "causalmlr_ate")) ate_pred <- ate_pred$estimate
  abs(ate_true - ate_pred)
}

#' Precision in estimation of heterogeneous effects (PEHE)
#'
#' The root mean squared error between true and predicted conditional
#' average treatment effects:
#' \deqn{PEHE = \sqrt{\frac{1}{n} \sum_i (\tau(x^{(i)}) -
#'   \hat\tau(x^{(i)}))^2}}
#' Requires the true CATEs, so it is only usable with (semi-)synthetic
#' benchmark data such as [synth_test] or [ihdp_test].
#'
#' @param tau_true Numeric vector of true CATEs.
#' @param tau_pred Numeric vector of predicted CATEs (same length).
#'
#' @return A single non-negative number.
#'
#' @examples
#' pehe(c(1, 2, 3), c(1.1, 1.8, 3.2))
#'
#' @references
#' Hill, J. L. (2011). Bayesian nonparametric modeling for causal inference.
#' *Journal of Computational and Graphical Statistics*, 20(1), 217-240.
#'
#' @seealso [eps_ate()], [r_loss()]
#' @export
pehe <- function(tau_true, tau_pred) {
  if (length(tau_true) != length(tau_pred)) {
    stop("`tau_true` and `tau_pred` must have the same length.", call. = FALSE)
  }
  sqrt(mean((tau_true - tau_pred)^2))
}

#' Nuisance models for the R-Loss
#'
#' Estimates the two nuisance functions required by the R-Loss ([r_loss()]):
#' the conditional outcome mean \eqn{m(x) = E[Y | X = x]} and the propensity
#' score \eqn{e(x) = P(T = 1 | X = x)}, and returns the corresponding
#' residuals \eqn{Y - \hat m(X)} and \eqn{T - \hat e(X)}. By default the
#' nuisance predictions are cross-fitted (out-of-fold), which avoids
#' overfitting bias; set `folds = 1` for simple in-sample estimates.
#'
#' Fitting the nuisance models once and reusing the resulting object to
#' score many candidate CATE models (e.g. across a hyperparameter grid) is
#' both faster and methodologically cleaner, since all candidates are
#' compared against the same residuals.
#'
#' @param data A `data.frame` with the outcome, treatment and covariates.
#'   Typically a held-out validation set that was not used to train the
#'   CATE models being evaluated.
#' @param outcome Name of the outcome column. Default `"y"`.
#' @param treatment Name of the binary (0/1) treatment column. Default `"t"`.
#' @param outcome_learner An mlr3 regression learner for `m(x)` (the
#'   treatment is excluded from the features).
#' @param ps_learner An mlr3 classification learner for `e(x)`.
#' @param covariates Optional character vector of covariate columns. Defaults
#'   to all columns except the outcome and the treatment.
#' @param folds Number of cross-fitting folds. Default `5`; `1` disables
#'   cross-fitting.
#' @param ps_trim Optional trimming for the estimated propensity scores;
#'   see [ate_ipw()].
#'
#' @return An object of class `"rloss_nuisance"` with elements `y_res`,
#'   `t_res`, `m_hat` and `e_hat`, to be passed to [r_loss()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' set.seed(1)
#' nuis <- rloss_nuisance(synth_train, outcome = "y", treatment = "t",
#'                        outcome_learner = lrn("regr.rpart"),
#'                        ps_learner = lrn("classif.rpart"))
#' m <- s_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"))
#' r_loss(nuis, predict(m, synth_train))
#' }
#'
#' @seealso [r_loss()]
#' @export
rloss_nuisance <- function(data, outcome = "y", treatment = "t",
                           outcome_learner, ps_learner, covariates = NULL,
                           folds = 5, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  fold_id <- make_folds(nrow(data), folds)
  e_hat <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id,
                       ps_trim)
  m_hat <- crossfit_outcome_mean(data, outcome, covariates, outcome_learner,
                                 fold_id)
  structure(
    list(y_res = y - m_hat, t_res = t - e_hat, m_hat = m_hat, e_hat = e_hat,
         n = nrow(data), folds = folds),
    class = "rloss_nuisance"
  )
}

#' R-Loss: an observable score for CATE models
#'
#' Computes the R-Loss (Nie & Wager, 2021), also known as
#' \eqn{\tau\text{-risk}_R}, of a vector of CATE predictions:
#' \deqn{R\text{-}Loss = \frac{1}{n} \sum_i \left[ (Y^{(i)} -
#'   \hat m(X^{(i)})) - (T^{(i)} - \hat e(X^{(i)})) \,
#'   \hat\tau(X^{(i)}) \right]^2}
#' Unlike [pehe()], the R-Loss does not require the true treatment effects,
#' so it can be computed on real observational data. This makes it suitable
#' for hyperparameter tuning and model selection of causal estimators:
#' lower values indicate better CATE models.
#'
#' @param nuisance An `"rloss_nuisance"` object created with
#'   [rloss_nuisance()], holding the outcome and treatment residuals.
#' @param tau_pred Numeric vector of CATE predictions for the same rows that
#'   `nuisance` was computed on.
#'
#' @return A single non-negative number (lower is better).
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' set.seed(1)
#' nuis <- rloss_nuisance(synth_train, outcome = "y", treatment = "t",
#'                        outcome_learner = lrn("regr.rpart"),
#'                        ps_learner = lrn("classif.rpart"))
#'
#' # compare two candidate CATE models by their R-Loss
#' m_s <- s_learner(synth_train, outcome = "y", treatment = "t",
#'                  learner = lrn("regr.rpart"))
#' m_t <- t_learner(synth_train, outcome = "y", treatment = "t",
#'                  learner = lrn("regr.rpart"))
#' r_loss(nuis, predict(m_s, synth_train))
#' r_loss(nuis, predict(m_t, synth_train))
#' }
#'
#' @references
#' Nie, X., & Wager, S. (2021). Quasi-oracle estimation of heterogeneous
#' treatment effects. *Biometrika*, 108(2), 299-319.
#'
#' Machlanski, D., Samothrakis, S., & Clarke, P. (2023). Hyperparameter
#' tuning and model evaluation in causal effect estimation. *arXiv preprint
#' arXiv:2303.01412*.
#'
#' @seealso [rloss_nuisance()], [pehe()], [eps_ate()]
#' @export
r_loss <- function(nuisance, tau_pred) {
  if (!inherits(nuisance, "rloss_nuisance")) {
    stop("`nuisance` must be created with rloss_nuisance().", call. = FALSE)
  }
  if (length(tau_pred) != nuisance$n) {
    stop(sprintf("`tau_pred` has length %d but `nuisance` was computed on %d rows.",
                 length(tau_pred), nuisance$n), call. = FALSE)
  }
  mean((nuisance$y_res - nuisance$t_res * tau_pred)^2)
}

#' @export
print.rloss_nuisance <- function(x, ...) {
  cat("R-Loss nuisance residuals\n")
  cat("  N:     ", x$n, "\n", sep = "")
  cat("  Folds: ", x$folds, "\n", sep = "")
  cat("Pass this object to r_loss() together with CATE predictions.\n")
  invisible(x)
}
