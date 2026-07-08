# Meta-learners for conditional average treatment effect (CATE) estimation.

new_cate_learner <- function(subclass, models, outcome, treatment, covariates) {
  structure(
    list(models = models, outcome = outcome, treatment = treatment,
         covariates = covariates),
    class = c(subclass, "cate_learner")
  )
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
                   covariates)
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
                   covariates)
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
                   outcome, treatment, covariates)
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

#' @export
print.cate_learner <- function(x, ...) {
  labels <- c(s_learner = "S-Learner", t_learner = "T-Learner",
              x_learner = "X-Learner")
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
