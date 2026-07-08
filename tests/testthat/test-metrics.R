test_that("eps_ate computes the absolute error", {
  expect_equal(eps_ate(1.05, 1.05), 0)
  expect_equal(eps_ate(1.05, 2.05), 1)
  expect_equal(eps_ate(2.05, 1.05), 1)
})

test_that("eps_ate unwraps causalmlr_ate objects", {
  est <- causalmlr:::new_ate("naive", 1.5, 0.1, 100)
  expect_equal(eps_ate(1, est), 0.5)
})

test_that("pehe computes the RMSE of CATE predictions", {
  expect_equal(pehe(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_equal(pehe(c(0, 0), c(3, 4)), sqrt(mean(c(9, 16))))
  expect_error(pehe(1:3, 1:2), "same length")
})

test_that("r_loss validates its inputs", {
  expect_error(r_loss(list(), 1:5), "rloss_nuisance")
  nuis <- structure(list(y_res = rnorm(10), t_res = rnorm(10), n = 10),
                    class = "rloss_nuisance")
  expect_error(r_loss(nuis, 1:5), "length")
  expect_true(is.numeric(r_loss(nuis, rep(0, 10))))
})

test_that("r_loss is zero for a perfect fit", {
  tau <- rnorm(10)
  t_res <- rnorm(10)
  nuis <- structure(list(y_res = t_res * tau, t_res = t_res, n = 10),
                    class = "rloss_nuisance")
  expect_equal(r_loss(nuis, tau), 0)
})
