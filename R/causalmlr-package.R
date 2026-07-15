#' causalmlr: Causal Machine Learning with mlr3
#'
#' Estimation and evaluation of average and conditional average treatment
#' effects (ATE/CATE) from observational data using machine learning, built
#' on the mlr3 ecosystem.
#'
#' @section ATE estimators:
#' * [ate_naive()] - difference in means (baseline)
#' * [ate_ipw()] - inverse propensity weighting
#' * [ate_dr()] - doubly robust / AIPW
#' * [ate_dml()] - double machine learning (partially linear model)
#'
#' @section CATE meta-learners:
#' * [s_learner()], [t_learner()], [x_learner()], [dr_learner()] (doubly
#'   robust), [r_learner()] (R-Loss / partially linear DML)
#'
#' @section Evaluation:
#' * [eps_ate()] - absolute ATE error (needs ground truth)
#' * [pehe()] - precision in estimating heterogeneous effects (needs ground truth)
#' * [r_loss()] with [rloss_nuisance()] - observable score for model
#'   selection and hyperparameter tuning on real data
#'
#' All estimators accept arbitrary mlr3 learners as nuisance models, so
#' anything from `mlr3::lrn("regr.lm")` to a tuned `AutoTuner` can be
#' plugged in.
#'
#' @keywords internal
"_PACKAGE"
