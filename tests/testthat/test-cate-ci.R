test_that("cate_ci returns well-formed intervals for every meta-learner", {
  skip_if_no_learners()
  d <- sim_hetero(n = 800)
  train <- d[1:600, ]
  test <- d[601:800, c("x1", "x2")]
  cov <- c("x1", "x2")

  learners <- list(
    s = s_learner(train, learner = regr_learner(), covariates = cov),
    t = t_learner(train, learner = regr_learner(), covariates = cov),
    x = x_learner(train, learner = regr_learner(),
                  ps_learner = classif_learner(), covariates = cov),
    dr = dr_learner(train, outcome_learner = regr_learner(),
                    ps_learner = classif_learner(), covariates = cov,
                    folds = 3, ps_trim = 0.05),
    r = r_learner(train, outcome_learner = regr_learner(),
                  ps_learner = classif_learner(), covariates = cov,
                  folds = 3, ps_trim = 0.05)
  )

  for (nm in names(learners)) {
    set.seed(1)
    ci <- cate_ci(learners[[nm]], test, train_data = train, n_boot = 10)
    expect_s3_class(ci, "data.frame")
    expect_named(ci, c("estimate", "se", "lower", "upper"))
    expect_equal(nrow(ci), nrow(test))
    expect_true(all(is.finite(ci$estimate)))
    expect_true(all(is.finite(ci$se)) && all(ci$se >= 0))
    expect_true(all(ci$lower <= ci$upper))
    # The point estimate must match a plain predict() call.
    expect_equal(ci$estimate, predict(learners[[nm]], test))
  }
})

test_that("cate_ci normal intervals are centred on the point estimate", {
  skip_if_no_learners()
  d <- sim_hetero(n = 600)
  train <- d[1:400, ]
  test <- d[401:600, c("x1", "x2")]
  m <- t_learner(train, learner = regr_learner(), covariates = c("x1", "x2"))
  set.seed(2)
  ci <- cate_ci(m, test, train_data = train, n_boot = 20, type = "normal")
  expect_true(all(ci$lower <= ci$estimate & ci$estimate <= ci$upper))
  # Symmetry around the point estimate.
  expect_equal(ci$estimate - ci$lower, ci$upper - ci$estimate)
})

test_that("a higher confidence level widens the interval", {
  skip_if_no_learners()
  d <- sim_hetero(n = 600)
  train <- d[1:400, ]
  test <- d[401:600, c("x1", "x2")]
  m <- s_learner(train, learner = regr_learner(), covariates = c("x1", "x2"))
  set.seed(3)
  ci90 <- cate_ci(m, test, train_data = train, n_boot = 30, level = 0.90,
                  type = "normal")
  set.seed(3)
  ci99 <- cate_ci(m, test, train_data = train, n_boot = 30, level = 0.99,
                  type = "normal")
  expect_true(all((ci99$upper - ci99$lower) >= (ci90$upper - ci90$lower)))
})

test_that("cate_ci validates its inputs", {
  skip_if_no_learners()
  d <- sim_hetero(n = 300)
  m <- s_learner(d, learner = regr_learner(), covariates = c("x1", "x2"))
  expect_error(cate_ci(unclass(m), d, train_data = d), "cate_learner")
  expect_error(cate_ci(m, d, train_data = d[, c("x1", "x2")]),
               "missing column")
  expect_error(cate_ci(m, d, train_data = d, level = 1.5), "level")
  expect_error(cate_ci(m, d, train_data = d, n_boot = 0), "n_boot")
})
