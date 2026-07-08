# Average treatment effect (ATE) estimators.

new_ate <- function(method, estimate, se, n, details = list()) {
  structure(
    list(method = method, estimate = estimate, se = se, n = n,
         details = details),
    class = "causalmlr_ate"
  )
}

#' Naive ATE estimator (difference in means)
#'
#' Estimates the average treatment effect as the difference between the mean
#' outcome of the treated and the mean outcome of the controls. This estimator
#' is unbiased only when the treatment is randomised; under confounding it is
#' typically biased and serves as a baseline for the adjusted estimators
#' ([ate_ipw()], [ate_dr()], [ate_dml()]).
#'
#' @param data A `data.frame` with the outcome, treatment and covariates.
#' @param outcome Name of the outcome column. Default `"y"`.
#' @param treatment Name of the binary (0/1) treatment column. Default `"t"`.
#'
#' @return An object of class `"causalmlr_ate"` with elements `estimate`,
#'   `se` (standard error), `n` and `method`. Use [coef()] to extract the
#'   point estimate and [confint()] for a confidence interval.
#'
#' @examples
#' data(sodium)
#' ate_naive(sodium, outcome = "bp", treatment = "sodium")
#'
#' @seealso [ate_ipw()], [ate_dr()], [ate_dml()]
#' @export
ate_naive <- function(data, outcome = "y", treatment = "t") {
  assert_data(data, outcome, treatment)
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  y1 <- y[t == 1L]
  y0 <- y[t == 0L]
  estimate <- mean(y1) - mean(y0)
  se <- sqrt(stats::var(y1) / length(y1) + stats::var(y0) / length(y0))
  new_ate("naive", estimate, se, nrow(data))
}

#' Inverse propensity weighting (IPW) ATE estimator
#'
#' Estimates the average treatment effect by weighting outcomes with the
#' inverse of the estimated propensity score
#' \eqn{e(x) = P(T = 1 | X = x)}:
#' \deqn{\widehat{ATE} = \frac{1}{n} \sum_i \frac{T_i Y_i}{\hat e(X_i)} -
#'   \frac{1}{n} \sum_i \frac{(1 - T_i) Y_i}{1 - \hat e(X_i)}}
#' The propensity score is estimated with any mlr3 classification learner.
#'
#' @inheritParams ate_naive
#' @param ps_learner An mlr3 classification learner used as the propensity
#'   score model, e.g. `mlr3::lrn("classif.log_reg")`. Its predict type is
#'   set to `"prob"` automatically.
#' @param covariates Optional character vector of covariate columns. Defaults
#'   to all columns except the outcome and the treatment.
#' @param folds Number of cross-fitting folds for the propensity model.
#'   `1` (the default) fits and predicts in-sample; values greater than 1
#'   use out-of-fold predictions, which reduces overfitting bias when using
#'   flexible learners.
#' @param ps_trim Optional trimming for the estimated propensity scores.
#'   Either a single value `a` (clips to `[a, 1 - a]`) or a length-2 range.
#'   `NULL` (default) applies no trimming.
#'
#' @return An object of class `"causalmlr_ate"`; see [ate_naive()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(sodium)
#' ate_ipw(sodium, outcome = "bp", treatment = "sodium",
#'         ps_learner = lrn("classif.rpart"))
#' }
#'
#' @seealso [ate_dr()], [ate_dml()]
#' @export
ate_ipw <- function(data, outcome = "y", treatment = "t", ps_learner,
                    covariates = NULL, folds = 1, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  fold_id <- make_folds(nrow(data), folds)
  e <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id, ps_trim)
  psi <- (y * t) / e - (y * (1L - t)) / (1 - e)
  new_ate("ipw", mean(psi), stats::sd(psi) / sqrt(length(psi)), nrow(data),
          details = list(ps = e, folds = folds))
}

#' Doubly robust (AIPW) ATE estimator
#'
#' Combines an outcome model \eqn{\mu(x, t) = E[Y | X = x, T = t]} with a
#' propensity score model \eqn{e(x) = P(T = 1 | X = x)} in the augmented
#' inverse propensity weighting (AIPW) estimator:
#' \deqn{\widehat{ATE} = \frac{1}{n} \sum_i \left[ \hat\mu_1(X_i) -
#'   \hat\mu_0(X_i) + \frac{T_i (Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} -
#'   \frac{(1 - T_i)(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)} \right]}
#' The estimator is consistent if either the outcome model or the propensity
#' model is correctly specified (hence "doubly robust"). With `folds > 1`
#' both nuisance models are cross-fitted, which yields the cross-fitted AIPW
#' estimator (also known as DML for the interactive/fully heterogeneous
#' model).
#'
#' @inheritParams ate_ipw
#' @param outcome_learner An mlr3 regression learner for the outcome model,
#'   e.g. `mlr3::lrn("regr.ranger")`. The treatment indicator is included as
#'   a feature and potential outcomes are predicted by setting it to 0 and 1.
#'
#' @return An object of class `"causalmlr_ate"`; see [ate_naive()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(sodium)
#' ate_dr(sodium, outcome = "bp", treatment = "sodium",
#'        outcome_learner = lrn("regr.rpart"),
#'        ps_learner = lrn("classif.rpart"),
#'        folds = 5, ps_trim = 0.01)
#' }
#'
#' @references
#' Robins, J. M., Rotnitzky, A., & Zhao, L. P. (1994). Estimation of
#' regression coefficients when some regressors are not always observed.
#' *Journal of the American Statistical Association*, 89(427), 846-866.
#'
#' @seealso [ate_ipw()], [ate_dml()]
#' @export
ate_dr <- function(data, outcome = "y", treatment = "t", outcome_learner,
                   ps_learner, covariates = NULL, folds = 1, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  fold_id <- make_folds(nrow(data), folds)
  e <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id, ps_trim)
  po <- crossfit_outcome_po(data, outcome, treatment, covariates,
                            outcome_learner, fold_id)
  psi <- po$mu1 - po$mu0 +
    (t * (y - po$mu1)) / e -
    ((1L - t) * (y - po$mu0)) / (1 - e)
  new_ate("dr", mean(psi), stats::sd(psi) / sqrt(length(psi)), nrow(data),
          details = list(ps = e, mu0 = po$mu0, mu1 = po$mu1, folds = folds))
}

#' Double machine learning (DML) ATE estimator
#'
#' Estimates the average treatment effect in the partially linear model
#' \eqn{Y = \theta T + g(X) + \varepsilon} using the residual-on-residual
#' (orthogonalised) estimator of Chernozhukov et al. (2018). Two nuisance
#' functions are estimated with machine learning and cross-fitting:
#' the conditional outcome mean \eqn{m(x) = E[Y | X = x]} and the propensity
#' score \eqn{e(x) = P(T = 1 | X = x)}. The ATE is then
#' \deqn{\hat\theta = \frac{\sum_i (T_i - \hat e(X_i))(Y_i - \hat m(X_i))}
#'   {\sum_i (T_i - \hat e(X_i))^2}}
#'
#' Note that unlike [ate_dr()], the outcome model here regresses `Y` on the
#' covariates only (the treatment is excluded), and the model assumes a
#' constant (homogeneous) treatment effect.
#'
#' @inheritParams ate_ipw
#' @param outcome_learner An mlr3 regression learner for the conditional
#'   outcome mean `m(x) = E[Y | X]` (treatment excluded from the features).
#' @param folds Number of cross-fitting folds. Default `5`. Cross-fitting is
#'   a core ingredient of DML; `folds = 1` (in-sample nuisance estimates) is
#'   allowed but not recommended with flexible learners.
#'
#' @return An object of class `"causalmlr_ate"`; see [ate_naive()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(sodium)
#' set.seed(1)
#' ate_dml(sodium, outcome = "bp", treatment = "sodium",
#'         outcome_learner = lrn("regr.rpart"),
#'         ps_learner = lrn("classif.rpart"),
#'         ps_trim = 0.01)
#' }
#'
#' @references
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
#' Newey, W., & Robins, J. (2018). Double/debiased machine learning for
#' treatment and structural parameters. *The Econometrics Journal*, 21(1),
#' C1-C68.
#'
#' @seealso [ate_dr()], [ate_ipw()]
#' @export
ate_dml <- function(data, outcome = "y", treatment = "t", outcome_learner,
                    ps_learner, covariates = NULL, folds = 5, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  n <- nrow(data)
  fold_id <- make_folds(n, folds)
  e <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id, ps_trim)
  m <- crossfit_outcome_mean(data, outcome, covariates, outcome_learner,
                             fold_id)
  y_res <- y - m
  t_res <- t - e
  estimate <- sum(t_res * y_res) / sum(t_res^2)
  psi <- (y_res - estimate * t_res) * t_res
  se <- sqrt(mean(psi^2) / n) / mean(t_res^2)
  new_ate("dml", estimate, se, n,
          details = list(ps = e, m = m, y_res = y_res, t_res = t_res,
                         folds = folds))
}

#' @export
print.causalmlr_ate <- function(x, digits = 4, ...) {
  labels <- c(naive = "Naive (difference in means)",
              ipw = "Inverse propensity weighting",
              dr = "Doubly robust (AIPW)",
              dml = "Double machine learning")
  label <- if (x$method %in% names(labels)) labels[[x$method]] else x$method
  cat("ATE estimate -", label, "\n")
  cat("  Estimate:  ", format(x$estimate, digits = digits), "\n", sep = "")
  if (is.finite(x$se)) {
    ci <- confint(x)
    cat("  Std. error:", format(x$se, digits = digits), "\n")
    cat("  95% CI:    [", format(ci[[1L]], digits = digits), ", ",
        format(ci[[2L]], digits = digits), "]\n", sep = "")
  }
  cat("  N:         ", x$n, "\n", sep = "")
  invisible(x)
}

#' @export
coef.causalmlr_ate <- function(object, ...) {
  stats::setNames(object$estimate, "ATE")
}

#' @export
confint.causalmlr_ate <- function(object, parm, level = 0.95, ...) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- object$estimate + c(-1, 1) * z * object$se
  stats::setNames(ci, c("lower", "upper"))
}
