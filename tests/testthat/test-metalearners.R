test_that("meta-learners fit, predict and beat a zero baseline", {
  skip_if_no_learners()
  d <- sim_hetero(n = 2000)
  train <- d[1:1500, ]
  test <- d[1501:2000, ]
  tau_true <- test$tau

  models <- list(
    s = s_learner(train, learner = regr_learner(),
                  covariates = c("x1", "x2")),
    t = t_learner(train, learner = regr_learner(),
                  covariates = c("x1", "x2")),
    x = x_learner(train, learner = regr_learner(),
                  ps_learner = classif_learner(),
                  covariates = c("x1", "x2")),
    dr = dr_learner(train, outcome_learner = regr_learner(),
                    ps_learner = classif_learner(),
                    covariates = c("x1", "x2"), folds = 5),
    r = r_learner(train, outcome_learner = regr_learner(),
                  ps_learner = classif_learner(),
                  covariates = c("x1", "x2"), folds = 5, ps_trim = 0.05)
  )

  pehe_zero <- pehe(tau_true, rep(0, nrow(test)))
  for (name in names(models)) {
    m <- models[[name]]
    expect_s3_class(m, "cate_learner")
    tau_hat <- predict(m, test)
    expect_length(tau_hat, nrow(test))
    expect_true(all(is.finite(tau_hat)))
    expect_lt(pehe(tau_true, tau_hat), pehe_zero)
    expect_equal(ate(m, test), mean(tau_true), tolerance = 0.5)
  }
})

test_that("predict works when newdata lacks outcome and treatment columns", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- s_learner(d, learner = regr_learner(), covariates = c("x1", "x2"))
  x_only <- d[, c("x1", "x2")]
  expect_length(predict(m, x_only), nrow(x_only))
  expect_error(predict(m, d[, "x1", drop = FALSE]), "missing covariate")
})

test_that("t_learner accepts a separate treated-group learner", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- t_learner(d, learner = regr_learner(),
                 learner1 = mlr3::lrn("regr.featureless"),
                 covariates = c("x1", "x2"))
  expect_equal(m$models$mu1$id, "regr.featureless")
  expect_length(predict(m, d), nrow(d))
})

test_that("fitting does not mutate the user's learner or data", {
  skip_if_no_learners()
  d <- sim_hetero(n = 500)
  d_copy <- d
  lr <- regr_learner()
  m <- s_learner(d, learner = lr, covariates = c("x1", "x2"))
  expect_identical(d, d_copy)
  expect_null(lr$model)  # user's learner instance stays untrained
})

test_that("dr_learner accepts a separate second-stage learner", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- dr_learner(d, outcome_learner = regr_learner(),
                  ps_learner = classif_learner(),
                  tau_learner = mlr3::lrn("regr.featureless"),
                  covariates = c("x1", "x2"), folds = 2, ps_trim = 0.05)
  expect_equal(m$models$tau$id, "regr.featureless")
  expect_length(predict(m, d), nrow(d))
  # A featureless second stage predicts the mean pseudo-outcome (~ ATE).
  expect_equal(ate(m, d), mean(d$tau), tolerance = 0.5)
})

test_that("dr_learner works without cross-fitting (folds = 1)", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- dr_learner(d, outcome_learner = regr_learner(),
                  ps_learner = classif_learner(),
                  covariates = c("x1", "x2"), folds = 1)
  expect_s3_class(m, "dr_learner")
  expect_true(all(is.finite(predict(m, d))))
})

test_that("r_learner accepts a separate second-stage learner", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- r_learner(d, outcome_learner = regr_learner(),
                 ps_learner = classif_learner(),
                 tau_learner = mlr3::lrn("regr.featureless"),
                 covariates = c("x1", "x2"), folds = 2, ps_trim = 0.05)
  expect_equal(m$models$tau$id, "regr.featureless")
  expect_length(predict(m, d), nrow(d))
})

test_that("r_learner works without cross-fitting (folds = 1)", {
  skip_if_no_learners()
  d <- sim_hetero(n = 1000)
  m <- r_learner(d, outcome_learner = regr_learner(),
                 ps_learner = classif_learner(),
                 covariates = c("x1", "x2"), folds = 1, ps_trim = 0.05)
  expect_s3_class(m, "r_learner")
  expect_true(all(is.finite(predict(m, d))))
})

test_that("r_learner falls back to unweighted fit without weight support", {
  skip_if_no_learners()
  d <- sim_hetero(n = 800)
  # A learner that does not support observation weights.
  no_w <- regr_no_weights()
  expect_warning(
    m <- r_learner(d, outcome_learner = regr_learner(),
                   ps_learner = classif_learner(),
                   tau_learner = no_w,
                   covariates = c("x1", "x2"), folds = 2, ps_trim = 0.05),
    "does not support observation weights"
  )
  expect_s3_class(m, "r_learner")
  expect_true(all(is.finite(predict(m, d))))
})

test_that("print method describes the model", {
  skip_if_no_learners()
  d <- sim_hetero(n = 500)
  m <- x_learner(d, learner = regr_learner(),
                 ps_learner = classif_learner(),
                 covariates = c("x1", "x2"))
  expect_output(print(m), "X-Learner")
})
