# Meta-learners for conditional average treatment effect (CATE) estimation.

new_cate_learner <- function(subclass, models, outcome, treatment, covariates,
                             refit = NULL) {
  structure(
    list(models = models, outcome = outcome, treatment = treatment,
         covariates = covariates, refit = refit),
    class = c(subclass, "cate_learner")
  )
}

# Build a closure that refits the same learner pipeline on new data. It
# captures only the (small) learner templates and settings -- never the
# training data -- so fitted objects stay lightweight. Used by [cate_ci()] to
# generate bootstrap replicates uniformly across the meta-learners.
make_refit <- function(fn, ...) {
  args <- list(...)
  function(data) do.call(fn, c(list(data), args))
}

# Drop outcome/treatment columns that may or may not be present in newdata.
strip_newdata <- function(object, newdata) {
  missing <- setdiff(object$covariates, names(newdata))
  if (length(missing) > 0L) {
    stop(sprintf("`newdata` is missing covariate column(s): %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  newdata[, object$covariates, drop = FALSE]
}

#' S-Learner for CATE estimation
#'
#' The S-Learner ("single learner") fits one outcome model
#' \eqn{\mu(x, t) = E[Y | X = x, T = t]} that includes the treatment
#' indicator as a regular feature. CATEs are obtained by contrasting the
#' predictions with the treatment switched on and off:
#' \deqn{\hat\tau(x) = \hat\mu(x, 1) - \hat\mu(x, 0)}
#'
#' @param data A `data.frame` with the outcome, treatment and covariates.
#' @param outcome Name of the outcome column. Default `"y"`.
#' @param treatment Name of the binary (0/1) treatment column. Default `"t"`.
#' @param learner An mlr3 regression learner for the outcome model, e.g.
#'   `mlr3::lrn("regr.ranger")`.
#' @param covariates Optional character vector of covariate columns. Defaults
#'   to all columns except the outcome and the treatment.
#'
#' @return A fitted CATE model of class `c("s_learner", "cate_learner")`.
#'   Use [predict()] to obtain CATE estimates for new data and [ate()] for
#'   their average.
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' m <- s_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"))
#' tau_hat <- predict(m, synth_test)
#' pehe(synth_test$tau, tau_hat)
#' }
#'
#' @references
#' Künzel, S. R., Sekhon, J. S., Bickel, P. J., & Yu, B. (2019).
#' Metalearners for estimating heterogeneous treatment effects using machine
#' learning. *Proceedings of the National Academy of Sciences*, 116(10),
#' 4156-4165.
#'
#' @seealso [t_learner()], [x_learner()], [pehe()], [r_loss()]
#' @export
s_learner <- function(data, outcome = "y", treatment = "t", learner,
                      covariates = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  df <- data[, c(covariates, treatment, outcome), drop = FALSE]
  df[[treatment]] <- as.integer(df[[treatment]])
  model <- fit_regr(df, outcome, learner)
  new_cate_learner("s_learner", list(mu = model), outcome, treatment,
                   covariates,
                   refit = make_refit(s_learner, outcome = outcome,
                                      treatment = treatment, learner = learner,
                                      covariates = covariates))
}

#' @describeIn s_learner Predict CATEs for new data. Returns a numeric
#'   vector \eqn{\hat\tau(x)} with one element per row of `newdata`
#'   (only the covariate columns are used).
#' @param object A fitted `s_learner` model.
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @param ... Ignored.
#' @export
predict.s_learner <- function(object, newdata, ...) {
  nd <- strip_newdata(object, newdata)
  nd[[object$treatment]] <- 0L
  y0 <- predict_regr(object$models$mu, nd)
  nd[[object$treatment]] <- 1L
  y1 <- predict_regr(object$models$mu, nd)
  y1 - y0
}

#' T-Learner for CATE estimation
#'
#' The T-Learner ("two learners") fits separate outcome models for the
#' control and treated groups,
#' \eqn{\mu_0(x) = E[Y | X = x, T = 0]} and
#' \eqn{\mu_1(x) = E[Y | X = x, T = 1]}, and estimates the CATE as their
#' difference:
#' \deqn{\hat\tau(x) = \hat\mu_1(x) - \hat\mu_0(x)}
#'
#' @inheritParams s_learner
#' @param learner An mlr3 regression learner used for the control-group
#'   model (and, if `learner1` is `NULL`, also for the treated-group model).
#' @param learner1 Optional separate mlr3 regression learner for the treated
#'   group. Defaults to a clone of `learner`.
#'
#' @return A fitted CATE model of class `c("t_learner", "cate_learner")`;
#'   see [s_learner()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' m <- t_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"))
#' tau_hat <- predict(m, synth_test)
#' pehe(synth_test$tau, tau_hat)
#' }
#'
#' @inherit s_learner references
#' @seealso [s_learner()], [x_learner()]
#' @export
t_learner <- function(data, outcome = "y", treatment = "t", learner,
                      learner1 = NULL, covariates = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  if (is.null(learner1)) learner1 <- learner
  t <- as.integer(data[[treatment]])
  df <- data[, c(covariates, outcome), drop = FALSE]
  m0 <- fit_regr(df[t == 0L, , drop = FALSE], outcome, learner)
  m1 <- fit_regr(df[t == 1L, , drop = FALSE], outcome, learner1)
  new_cate_learner("t_learner", list(mu0 = m0, mu1 = m1), outcome, treatment,
                   covariates,
                   refit = make_refit(t_learner, outcome = outcome,
                                      treatment = treatment, learner = learner,
                                      learner1 = learner1,
                                      covariates = covariates))
}

#' @describeIn t_learner Predict CATEs for new data.
#' @param object A fitted `t_learner` model.
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @param ... Ignored.
#' @export
predict.t_learner <- function(object, newdata, ...) {
  nd <- strip_newdata(object, newdata)
  y0 <- predict_regr(object$models$mu0, nd)
  y1 <- predict_regr(object$models$mu1, nd)
  y1 - y0
}

#' X-Learner for CATE estimation
#'
#' The X-Learner (Künzel et al., 2019) proceeds in three stages:
#' 1. Fit group-specific outcome models \eqn{\hat\mu_0(x)} and
#'    \eqn{\hat\mu_1(x)} as in the T-Learner.
#' 2. Impute individual treatment effects,
#'    \eqn{D_0 = \hat\mu_1(X_0) - Y_0} for controls and
#'    \eqn{D_1 = Y_1 - \hat\mu_0(X_1)} for the treated, and regress them on
#'    the covariates to obtain \eqn{\hat\tau_0(x)} and \eqn{\hat\tau_1(x)}.
#' 3. Combine the two estimates with propensity score weights:
#'    \deqn{\hat\tau(x) = \hat e(x) \hat\tau_0(x) +
#'      (1 - \hat e(x)) \hat\tau_1(x)}
#'
#' The X-Learner is particularly effective when the treated and control
#' groups are of very different sizes.
#'
#' @inheritParams s_learner
#' @param learner An mlr3 regression learner for the stage-1 outcome models.
#' @param ps_learner An mlr3 classification learner for the propensity score
#'   model, e.g. `mlr3::lrn("classif.log_reg")`.
#' @param tau_learner Optional mlr3 regression learner for the stage-2
#'   imputed-effect models. Defaults to a clone of `learner`.
#'
#' @return A fitted CATE model of class `c("x_learner", "cate_learner")`;
#'   see [s_learner()].
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' m <- x_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"),
#'                ps_learner = lrn("classif.rpart"))
#' tau_hat <- predict(m, synth_test)
#' pehe(synth_test$tau, tau_hat)
#' }
#'
#' @inherit s_learner references
#' @seealso [s_learner()], [t_learner()]
#' @export
x_learner <- function(data, outcome = "y", treatment = "t", learner,
                      ps_learner, tau_learner = NULL, covariates = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  if (is.null(tau_learner)) tau_learner <- learner
  t <- as.integer(data[[treatment]])
  df <- data[, c(covariates, outcome), drop = FALSE]
  df0 <- df[t == 0L, , drop = FALSE]
  df1 <- df[t == 1L, , drop = FALSE]

  # Stage 1: group-specific outcome models.
  m0 <- fit_regr(df0, outcome, learner)
  m1 <- fit_regr(df1, outcome, learner)

  # Stage 2: imputed treatment effects regressed on covariates.
  d0 <- df0
  d0[[outcome]] <- predict_regr(m1, df0[, covariates, drop = FALSE]) -
    df0[[outcome]]
  d1 <- df1
  d1[[outcome]] <- df1[[outcome]] -
    predict_regr(m0, df1[, covariates, drop = FALSE])
  tau0 <- fit_regr(d0, outcome, tau_learner)
  tau1 <- fit_regr(d1, outcome, tau_learner)

  # Stage 3: propensity score model for the weighting.
  dfe <- data[, c(covariates, treatment), drop = FALSE]
  e_model <- fit_propensity(dfe, treatment, ps_learner)

  new_cate_learner("x_learner",
                   list(mu0 = m0, mu1 = m1, tau0 = tau0, tau1 = tau1,
                        ps = e_model),
                   outcome, treatment, covariates,
                   refit = make_refit(x_learner, outcome = outcome,
                                      treatment = treatment, learner = learner,
                                      ps_learner = ps_learner,
                                      tau_learner = tau_learner,
                                      covariates = covariates))
}

#' @describeIn x_learner Predict CATEs for new data.
#' @param object A fitted `x_learner` model.
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @param ... Ignored.
#' @export
predict.x_learner <- function(object, newdata, ...) {
  nd <- strip_newdata(object, newdata)
  e <- predict_propensity(object$models$ps, nd)
  tau0 <- predict_regr(object$models$tau0, nd)
  tau1 <- predict_regr(object$models$tau1, nd)
  e * tau0 + (1 - e) * tau1
}

#' DR-Learner for CATE estimation
#'
#' The DR-Learner ("doubly robust learner"; Kennedy, 2023) is a two-stage
#' meta-learner that turns the augmented inverse propensity weighting (AIPW)
#' score used by [ate_dr()] into a model that predicts individual CATEs.
#' 1. Cross-fit the nuisance functions: the propensity score
#'    \eqn{\hat e(x) = P(T = 1 | X = x)} and the potential-outcome models
#'    \eqn{\hat\mu_0(x)}, \eqn{\hat\mu_1(x)} from a single outcome model that
#'    includes the treatment as a feature.
#' 2. Form the doubly robust pseudo-outcome
#'    \deqn{\psi_i = \hat\mu_1(X_i) - \hat\mu_0(X_i) +
#'      \frac{T_i (Y_i - \hat\mu_1(X_i))}{\hat e(X_i)} -
#'      \frac{(1 - T_i)(Y_i - \hat\mu_0(X_i))}{1 - \hat e(X_i)}}
#'    and regress it on the covariates to obtain the final CATE model
#'    \eqn{\hat\tau(x) = E[\psi | X = x]}.
#'
#' Because the pseudo-outcome has conditional mean equal to the true CATE
#' whenever either nuisance is correctly specified, the second-stage
#' regression targets the CATE directly. Cross-fitting the nuisances (the
#' default `folds = 5`) makes the pseudo-outcome orthogonal to nuisance
#' estimation error, so flexible learners can be used without overfitting
#' bias. Unlike [ate_dr()], which averages \eqn{\psi_i} to a scalar ATE, the
#' DR-Learner keeps the second-stage model and can therefore predict CATEs on
#' new test data.
#'
#' @inheritParams ate_dr
#' @param outcome_learner An mlr3 regression learner for the stage-1 outcome
#'   model, e.g. `mlr3::lrn("regr.ranger")`. The treatment indicator is
#'   included as a feature and potential outcomes are predicted by setting it
#'   to 0 and 1.
#' @param tau_learner Optional mlr3 regression learner for the second-stage
#'   pseudo-outcome regression. Defaults to a clone of `outcome_learner`.
#' @param folds Number of cross-fitting folds for the nuisance models.
#'   Default `5`. Cross-fitting is a core ingredient of the DR-Learner;
#'   `folds = 1` (in-sample nuisance estimates) is allowed but not recommended
#'   with flexible learners.
#'
#' @return A fitted CATE model of class `c("dr_learner", "cate_learner")`.
#'   Use [predict()] to obtain CATE estimates for new data and [ate()] for
#'   their average.
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' set.seed(1)
#' m <- dr_learner(synth_train, outcome = "y", treatment = "t",
#'                 outcome_learner = lrn("regr.rpart"),
#'                 ps_learner = lrn("classif.rpart"),
#'                 ps_trim = 0.01)
#' tau_hat <- predict(m, synth_test)
#' pehe(synth_test$tau, tau_hat)
#' }
#'
#' @references
#' Kennedy, E. H. (2023). Towards optimal doubly robust estimation of
#' heterogeneous causal effects. *Electronic Journal of Statistics*, 17(2),
#' 3008-3049.
#'
#' @seealso [s_learner()], [t_learner()], [x_learner()], [r_learner()],
#'   [ate_dr()]
#' @export
dr_learner <- function(data, outcome = "y", treatment = "t", outcome_learner,
                       ps_learner, tau_learner = NULL, covariates = NULL,
                       folds = 5, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  if (is.null(tau_learner)) tau_learner <- outcome_learner
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  fold_id <- make_folds(nrow(data), folds)

  # Stage 1: cross-fitted nuisance functions.
  e <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id, ps_trim)
  po <- crossfit_outcome_po(data, outcome, treatment, covariates,
                            outcome_learner, fold_id)

  # Stage 2: regress the doubly robust pseudo-outcome on the covariates.
  psi <- po$mu1 - po$mu0 +
    (t * (y - po$mu1)) / e -
    ((1L - t) * (y - po$mu0)) / (1 - e)
  d <- data[, covariates, drop = FALSE]
  d[[outcome]] <- psi
  tau <- fit_regr(d, outcome, tau_learner)

  new_cate_learner("dr_learner", list(tau = tau), outcome, treatment,
                   covariates,
                   refit = make_refit(dr_learner, outcome = outcome,
                                      treatment = treatment,
                                      outcome_learner = outcome_learner,
                                      ps_learner = ps_learner,
                                      tau_learner = tau_learner,
                                      covariates = covariates, folds = folds,
                                      ps_trim = ps_trim))
}

#' @describeIn dr_learner Predict CATEs for new data.
#' @param object A fitted `dr_learner` model.
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @param ... Ignored.
#' @export
predict.dr_learner <- function(object, newdata, ...) {
  nd <- strip_newdata(object, newdata)
  predict_regr(object$models$tau, nd)
}

#' R-Learner for CATE estimation
#'
#' The R-Learner (Nie & Wager, 2021) is the heterogeneous-effect
#' generalisation of the partially linear DML estimator ([ate_dml()]). It
#' builds on the Robinson decomposition of the model
#' \eqn{Y = \tau(X) T + g(X) + \varepsilon}:
#' 1. Cross-fit the two nuisance functions, the conditional outcome mean
#'    \eqn{\hat m(x) = E[Y | X = x]} (the treatment is excluded from the
#'    features) and the propensity score \eqn{\hat e(x) = P(T = 1 | X = x)},
#'    and form the residuals \eqn{\tilde Y = Y - \hat m(X)} and
#'    \eqn{\tilde T = T - \hat e(X)}.
#' 2. Estimate the CATE by minimising the R-Loss
#'    \deqn{\hat\tau = \arg\min_\tau \frac{1}{n} \sum_i \left[ \tilde Y_i -
#'      \tilde T_i \, \tau(X_i) \right]^2}
#'    which is fitted as a weighted regression of the pseudo-outcome
#'    \eqn{\tilde Y_i / \tilde T_i} on the covariates with weights
#'    \eqn{\tilde T_i^2}.
#'
#' The objective minimised here is exactly the [r_loss()] used elsewhere for
#' model selection, so the R-Learner is the estimator that directly targets
#' that score. When the treatment effect is constant it reduces to
#' [ate_dml()]; unlike that scalar estimator, the R-Learner keeps the
#' second-stage model and can predict CATEs on new test data.
#'
#' The weighted formulation requires a `tau_learner` that supports
#' observation weights (`"weights" %in% learner$properties`). If it does not,
#' the R-Learner falls back to an unweighted regression of the pseudo-outcome
#' and issues a warning. As with [ate_dml()], small \eqn{\tilde T_i} inflate
#' the pseudo-outcome, so `ps_trim` is recommended with flexible propensity
#' learners.
#'
#' @inheritParams dr_learner
#' @param outcome_learner An mlr3 regression learner for the conditional
#'   outcome mean \eqn{m(x) = E[Y | X]} (the treatment is excluded from the
#'   features), e.g. `mlr3::lrn("regr.ranger")`.
#' @param tau_learner Optional mlr3 regression learner for the second-stage
#'   weighted pseudo-outcome regression. Defaults to a clone of
#'   `outcome_learner`.
#'
#' @return A fitted CATE model of class `c("r_learner", "cate_learner")`.
#'   Use [predict()] to obtain CATE estimates for new data and [ate()] for
#'   their average.
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' set.seed(1)
#' m <- r_learner(synth_train, outcome = "y", treatment = "t",
#'                outcome_learner = lrn("regr.rpart"),
#'                ps_learner = lrn("classif.rpart"),
#'                ps_trim = 0.01)
#' tau_hat <- predict(m, synth_test)
#' pehe(synth_test$tau, tau_hat)
#' }
#'
#' @references
#' Nie, X., & Wager, S. (2021). Quasi-oracle estimation of heterogeneous
#' treatment effects. *Biometrika*, 108(2), 299-319.
#'
#' @seealso [s_learner()], [t_learner()], [x_learner()], [dr_learner()],
#'   [ate_dml()], [r_loss()]
#' @export
r_learner <- function(data, outcome = "y", treatment = "t", outcome_learner,
                      ps_learner, tau_learner = NULL, covariates = NULL,
                      folds = 5, ps_trim = NULL) {
  assert_data(data, outcome, treatment)
  covariates <- resolve_covariates(data, outcome, treatment, covariates)
  if (is.null(tau_learner)) tau_learner <- outcome_learner
  y <- data[[outcome]]
  t <- as.integer(data[[treatment]])
  fold_id <- make_folds(nrow(data), folds)

  # Stage 1: cross-fitted nuisances and Robinson residuals.
  e <- crossfit_ps(data, treatment, covariates, ps_learner, fold_id, ps_trim)
  m <- crossfit_outcome_mean(data, outcome, covariates, outcome_learner,
                             fold_id)
  y_res <- y - m
  t_res <- t - e

  # Stage 2: minimise the R-Loss via a weighted regression of the
  # pseudo-outcome (y_res / t_res) on the covariates with weights t_res^2.
  d <- data[, covariates, drop = FALSE]
  d[[outcome]] <- y_res / t_res
  w <- t_res^2
  if (regr_supports_weights(tau_learner)) {
    tau <- fit_regr_weighted(d, outcome, w, tau_learner)
  } else {
    warning(sprintf(
      paste0("`tau_learner` ('%s') does not support observation weights; ",
             "falling back to an unweighted R-Learner regression."),
      tau_learner$id), call. = FALSE)
    tau <- fit_regr(d, outcome, tau_learner)
  }

  new_cate_learner("r_learner", list(tau = tau), outcome, treatment,
                   covariates,
                   refit = make_refit(r_learner, outcome = outcome,
                                      treatment = treatment,
                                      outcome_learner = outcome_learner,
                                      ps_learner = ps_learner,
                                      tau_learner = tau_learner,
                                      covariates = covariates, folds = folds,
                                      ps_trim = ps_trim))
}

#' @describeIn r_learner Predict CATEs for new data.
#' @param object A fitted `r_learner` model.
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @param ... Ignored.
#' @export
predict.r_learner <- function(object, newdata, ...) {
  nd <- strip_newdata(object, newdata)
  predict_regr(object$models$tau, nd)
}

#' Average treatment effect of a fitted CATE model
#'
#' Averages the CATE predictions of a fitted meta-learner over the rows of
#' `newdata`:
#' \deqn{\widehat{ATE} = \frac{1}{n} \sum_i \hat\tau(x^{(i)})}
#'
#' @param object A fitted CATE model, e.g. from [s_learner()].
#' @param ... Passed on to methods.
#'
#' @return A single numeric value.
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' m <- s_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"))
#' ate(m, synth_train)
#' }
#'
#' @export
ate <- function(object, ...) {
  UseMethod("ate")
}

#' @rdname ate
#' @param newdata A `data.frame` containing at least the covariate columns.
#' @export
ate.cate_learner <- function(object, newdata, ...) {
  mean(predict(object, newdata))
}

#' Bootstrap confidence intervals for predicted CATEs
#'
#' Attaches pointwise confidence intervals to the CATE predictions of any
#' fitted meta-learner ([s_learner()], [t_learner()], [x_learner()],
#' [dr_learner()], [r_learner()]) using the nonparametric bootstrap. For each
#' of `n_boot` replicates the training data is resampled with replacement, the
#' entire learner pipeline (including any nuisance models and cross-fitting) is
#' refitted, and CATEs are predicted for `newdata`. A pointwise interval is
#' then formed from the bootstrap distribution at each row of `newdata`.
#'
#' Because the meta-learners accept arbitrary (and typically non-smooth)
#' machine-learning base learners, no closed-form standard error is available
#' for \eqn{\hat\tau(x)}; the bootstrap is the one mechanism that applies
#' uniformly across all five learners. Two caveats are worth keeping in mind.
#' First, bootstrap coverage for CATEs estimated with flexible learners is only
#' approximate and can be anti-conservative, so the intervals are best read as
#' a measure of estimation stability rather than exact frequentist coverage.
#' Second, refitting the full pipeline `n_boot` times is computationally
#' expensive. Set the seed with [set.seed()] beforehand for reproducibility.
#'
#' @param object A fitted meta-learner of class `"cate_learner"`.
#' @param newdata A `data.frame` of points at which to predict CATEs and their
#'   intervals; must contain the covariate columns.
#' @param train_data The `data.frame` the model was fitted on (with the
#'   outcome, treatment and covariate columns), used to draw the bootstrap
#'   resamples on which the pipeline is refitted.
#' @param n_boot Number of bootstrap replicates. Default `200`.
#' @param level Confidence level. Default `0.95`.
#' @param type Interval type. `"percentile"` (default) uses the empirical
#'   quantiles of the bootstrap CATEs; `"normal"` centres the interval on the
#'   point estimate and uses a normal-quantile multiple of the bootstrap
#'   standard error.
#' @param ... Ignored.
#'
#' @return A `data.frame` with one row per row of `newdata` and columns
#'   `estimate` (the point CATE from `object`), `se` (the bootstrap standard
#'   error) and `lower`/`upper` (the confidence limits).
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' data(synth_train)
#' data(synth_test)
#' set.seed(1)
#' m <- t_learner(synth_train, outcome = "y", treatment = "t",
#'                learner = lrn("regr.rpart"))
#' ci <- cate_ci(m, synth_test, train_data = synth_train, n_boot = 50)
#' head(ci)
#' }
#'
#' @seealso [s_learner()], [t_learner()], [x_learner()], [dr_learner()],
#'   [r_learner()]
#' @export
cate_ci <- function(object, newdata, train_data, n_boot = 200, level = 0.95,
                    type = c("percentile", "normal"), ...) {
  if (!inherits(object, "cate_learner")) {
    stop("`object` must be a fitted meta-learner (class 'cate_learner').",
         call. = FALSE)
  }
  if (!is.function(object$refit)) {
    stop("`object` does not carry a refit method; it may have been created by ",
         "an older version of the package. Refit the model to use `cate_ci()`.",
         call. = FALSE)
  }
  if (!is.data.frame(train_data)) {
    stop("`train_data` must be a data.frame.", call. = FALSE)
  }
  missing <- setdiff(c(object$covariates, object$outcome, object$treatment),
                     names(train_data))
  if (length(missing) > 0L) {
    stop(sprintf("`train_data` is missing column(s): %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 1) {
    stop("`n_boot` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  type <- match.arg(type)
  n_boot <- as.integer(n_boot)

  estimate <- predict(object, newdata)
  n <- nrow(train_data)
  boot <- matrix(NA_real_, nrow = length(estimate), ncol = n_boot)
  n_fail <- 0L
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    pred_b <- tryCatch(
      predict(object$refit(train_data[idx, , drop = FALSE]), newdata),
      error = function(e) NULL)
    # Drop degenerate replicates: a refit that errors, returns the wrong
    # length, or produces non-finite CATEs (e.g. extreme propensities on a
    # resample inflating a pseudo-outcome) would otherwise corrupt the
    # bootstrap distribution.
    if (is.null(pred_b) || length(pred_b) != length(estimate) ||
        !all(is.finite(pred_b))) {
      n_fail <- n_fail + 1L
    } else {
      boot[, b] <- pred_b
    }
  }
  if (n_fail == n_boot) {
    stop(sprintf(paste0("All %d bootstrap replicates failed to fit; cannot ",
                        "form confidence intervals."), n_boot), call. = FALSE)
  }
  if (n_fail > 0L) {
    warning(sprintf(
      "%d of %d bootstrap replicates failed to fit and were dropped.",
      n_fail, n_boot), call. = FALSE)
  }

  se <- apply(boot, 1L, stats::sd, na.rm = TRUE)
  a <- (1 - level) / 2
  if (type == "percentile") {
    q <- apply(boot, 1L, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE,
               names = FALSE)
    lower <- q[1L, ]
    upper <- q[2L, ]
  } else {
    z <- stats::qnorm(1 - a)
    lower <- estimate - z * se
    upper <- estimate + z * se
  }
  data.frame(estimate = estimate, se = se, lower = lower, upper = upper)
}

#' @export
print.cate_learner <- function(x, ...) {
  labels <- c(s_learner = "S-Learner", t_learner = "T-Learner",
              x_learner = "X-Learner", dr_learner = "DR-Learner",
              r_learner = "R-Learner")
  label <- labels[[class(x)[[1L]]]]
  cat(label, "(fitted)\n")
  cat("  Outcome:   ", x$outcome, "\n", sep = "")
  cat("  Treatment: ", x$treatment, "\n", sep = "")
  cat("  Covariates:", length(x$covariates), "\n")
  for (nm in names(x$models)) {
    cat("  Model '", nm, "': ", x$models[[nm]]$id, "\n", sep = "")
  }
  invisible(x)
}
