test_that("rloss_nuisance produces residuals and r_loss ranks models sensibly", {
  skip_if_no_learners()
  d <- sim_hetero(n = 2000)
  d$tau <- NULL

  set.seed(11)
  nuis <- rloss_nuisance(d, outcome_learner = regr_learner(),
                         ps_learner = classif_learner(), folds = 5, ps_trim = 0.01)
  expect_s3_class(nuis, "rloss_nuisance")
  expect_length(nuis$y_res, nrow(d))
  expect_length(nuis$t_res, nrow(d))
  expect_true(all(nuis$e_hat >= 0 & nuis$e_hat <= 1))

  # the true CATE should score better than a wildly wrong constant
  tau_true <- 2 * (d$x1 > 0)
  expect_lt(r_loss(nuis, tau_true), r_loss(nuis, rep(-10, nrow(d))))
})

test_that("rloss_nuisance works without cross-fitting", {
  skip_if_no_learners()
  d <- sim_hetero(n = 500)
  d$tau <- NULL
  nuis <- rloss_nuisance(d, outcome_learner = regr_learner(),
                         ps_learner = classif_learner(), folds = 1)
  expect_equal(nuis$folds, 1)
  expect_output(print(nuis), "R-Loss")
})
