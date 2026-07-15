# Simulated datasets with known ground truth used across the tests.

# Confounded data with a constant treatment effect.
sim_confounded <- function(n = 2000, ate = 2, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  e <- stats::plogis(0.8 * x1 - 0.5 * x2)
  t <- rbinom(n, 1, e)
  y <- ate * t + 1.5 * x1 + x2 + rnorm(n, sd = 0.5)
  data.frame(x1 = x1, x2 = x2, t = t, y = y)
}

# Randomised-ish data with a heterogeneous (step function) treatment effect.
sim_hetero <- function(n = 2000, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  tau <- 2 * (x1 > 0)
  t <- rbinom(n, 1, stats::plogis(0.5 * x2))
  y <- tau * t + x1 + 0.5 * x2 + rnorm(n, sd = 0.5)
  data.frame(x1 = x1, x2 = x2, t = t, y = y, tau = tau)
}

skip_if_no_learners <- function() {
  testthat::skip_if_not_installed("mlr3")
  testthat::skip_if_not_installed("rpart")
}

regr_learner <- function() mlr3::lrn("regr.rpart")
classif_learner <- function() mlr3::lrn("classif.rpart")

# A minimal regression learner that does NOT advertise the "weights" property,
# used to exercise the R-Learner's unweighted fallback. Unlike stripping the
# property off a real learner, its `.train` never touches observation weights.
regr_no_weights <- function() {
  cls <- R6::R6Class(
    "LearnerRegrNoWeights",
    inherit = mlr3::LearnerRegr,
    public = list(
      initialize = function() {
        super$initialize(
          id = "regr.no_weights",
          feature_types = c("logical", "integer", "numeric"),
          predict_types = "response",
          properties = character()
        )
      }
    ),
    private = list(
      .train = function(task) {
        x <- as.matrix(task$data(cols = task$feature_names))
        stats::lm.fit(cbind(1, x), task$truth())$coefficients
      },
      .predict = function(task) {
        x <- as.matrix(task$data(cols = task$feature_names))
        list(response = as.numeric(cbind(1, x) %*% self$model))
      }
    )
  )
  cls$new()
}
