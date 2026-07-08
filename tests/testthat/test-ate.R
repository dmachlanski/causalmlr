test_that("ate_naive recovers the ATE under randomisation", {
  set.seed(1)
  n <- 5000
  t <- rbinom(n, 1, 0.5)
  y <- 2 * t + rnorm(n)
  d <- data.frame(x = rnorm(n), t = t, y = y)
  est <- ate_naive(d)
  expect_s3_class(est, "causalmlr_ate")
  expect_equal(est$estimate, 2, tolerance = 0.1)
  expect_gt(est$se, 0)
  ci <- confint(est)
  expect_lt(ci[["lower"]], est$estimate)
  expect_gt(ci[["upper"]], est$estimate)
  expect_equal(unname(coef(est)), est$estimate)
})

test_that("ate_naive is biased under confounding but adjusted estimators are not", {
  skip_if_no_learners()
  d <- sim_confounded(n = 3000, ate = 2)
  naive <- ate_naive(d)
  # x1 raises both treatment probability and outcome -> upward bias
  expect_gt(naive$estimate, 2.3)

  set.seed(7)
  ipw <- ate_ipw(d, ps_learner = classif_learner(), folds = 5,
                 ps_trim = 0.01)
  dr <- ate_dr(d, outcome_learner = regr_learner(),
               ps_learner = classif_learner(), folds = 5, ps_trim = 0.01)
  dml <- ate_dml(d, outcome_learner = regr_learner(),
                 ps_learner = classif_learner(), folds = 5, ps_trim = 0.01)

  expect_lt(eps_ate(2, ipw), eps_ate(2, naive))
  expect_lt(eps_ate(2, dr), eps_ate(2, naive))
  expect_lt(eps_ate(2, dml), eps_ate(2, naive))
  expect_equal(dr$estimate, 2, tolerance = 0.5)
  expect_equal(dml$estimate, 2, tolerance = 0.5)
})

test_that("ATE estimators respect custom column names and covariates", {
  skip_if_no_learners()
  d <- sim_confounded(n = 1000)
  names(d) <- c("age", "income", "treated", "outcome")
  d$noise <- "a"  # non-covariate column excluded via `covariates`
  est <- ate_ipw(d, outcome = "outcome", treatment = "treated",
                 ps_learner = classif_learner(),
                 covariates = c("age", "income"))
  expect_s3_class(est, "causalmlr_ate")
  expect_true(is.finite(est$estimate))
})

test_that("input validation catches bad data", {
  d <- sim_confounded(n = 100)
  expect_error(ate_naive(d, outcome = "nope"), "not found")
  expect_error(ate_naive(d, treatment = "x1"), "binary")
  expect_error(ate_naive("not a df"), "data.frame")
  d_one_group <- d
  d_one_group$t <- 1
  expect_error(ate_naive(d_one_group), "both treated")
})

test_that("learner type mismatches are rejected", {
  skip_if_no_learners()
  d <- sim_confounded(n = 200)
  expect_error(
    ate_ipw(d, ps_learner = regr_learner()),
    "classif"
  )
  expect_error(
    ate_dr(d, outcome_learner = classif_learner(),
           ps_learner = classif_learner()),
    "regr"
  )
})

test_that("print methods run without error", {
  d <- sim_confounded(n = 200)
  expect_output(print(ate_naive(d)), "Naive")
})
