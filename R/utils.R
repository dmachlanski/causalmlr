# Internal helpers shared across estimators. None of these are exported.

# Validate the input data frame and the outcome/treatment columns.
assert_data <- function(data, outcome, treatment) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  for (col in c(outcome, treatment)) {
    if (!is.character(col) || length(col) != 1L) {
      stop("`outcome` and `treatment` must be single column names.", call. = FALSE)
    }
    if (!col %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", col), call. = FALSE)
    }
  }
  tvals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(tvals %in% c(0, 1))) {
    stop(sprintf("Treatment column '%s' must be binary (0/1).", treatment),
         call. = FALSE)
  }
  if (length(tvals) < 2L) {
    stop(sprintf("Treatment column '%s' must contain both treated (1) and control (0) units.",
                 treatment), call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve the covariate set: everything except outcome/treatment by default.
resolve_covariates <- function(data, outcome, treatment, covariates = NULL) {
  if (is.null(covariates)) {
    covariates <- setdiff(names(data), c(outcome, treatment))
  } else {
    missing <- setdiff(covariates, names(data))
    if (length(missing) > 0L) {
      stop(sprintf("Covariate column(s) not found in `data`: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    covariates <- setdiff(covariates, c(outcome, treatment))
  }
  if (length(covariates) == 0L) {
    stop("No covariates left after removing outcome and treatment columns.",
         call. = FALSE)
  }
  covariates
}

# Deep-clone a user-supplied mlr3 learner and check its task type.
clone_learner <- function(learner, type = c("regr", "classif"), arg = "learner") {
  type <- match.arg(type)
  if (!inherits(learner, "Learner")) {
    stop(sprintf("`%s` must be an mlr3 Learner, e.g. mlr3::lrn('%s.rpart').",
                 arg, type), call. = FALSE)
  }
  if (learner$task_type != type) {
    stop(sprintf("`%s` must be a %s learner, but '%s' is a %s learner.",
                 arg, type, learner$id, learner$task_type), call. = FALSE)
  }
  lr <- learner$clone(deep = TRUE)
  if (type == "classif" && lr$predict_type != "prob") {
    lr$predict_type <- "prob"
  }
  lr
}

# Fit a regression learner on `df` with the given target column.
fit_regr <- function(df, target, learner) {
  task <- mlr3::as_task_regr(df, target = target, id = "causalmlr")
  lr <- clone_learner(learner, "regr")
  lr$train(task)
  lr
}

predict_regr <- function(model, newdata) {
  model$predict_newdata(newdata = newdata)$response
}

# Does an mlr3 regression learner support observation weights?
regr_supports_weights <- function(learner) {
  "weights" %in% learner$properties
}

# Fit a regression learner with observation weights. The weights are attached
# through the mlr3 "weight" column role (renamed "weights_learner" in newer
# mlr3), which is detected from the task so both APIs are supported.
fit_regr_weighted <- function(df, target, weights, learner) {
  df[["..weights.."]] <- weights
  task <- mlr3::as_task_regr(df, target = target, id = "causalmlr")
  w_role <- if ("weights_learner" %in% names(task$col_roles)) {
    "weights_learner"
  } else {
    "weight"
  }
  task$set_col_roles("..weights..", roles = w_role)
  lr <- clone_learner(learner, "regr")
  lr$train(task)
  lr
}

# Fit a propensity (classification) learner; `df` must contain the treatment
# column, which is converted to a factor with levels c("0", "1").
fit_propensity <- function(df, treatment, learner) {
  df[[treatment]] <- factor(df[[treatment]], levels = c(0, 1))
  task <- mlr3::as_task_classif(df, target = treatment, positive = "1",
                                id = "causalmlr")
  lr <- clone_learner(learner, "classif", arg = "ps_learner")
  lr$train(task)
  lr
}

predict_propensity <- function(model, newdata) {
  model$predict_newdata(newdata = newdata)$prob[, "1"]
}

# Clip propensity scores to [trim, 1 - trim] (or a length-2 range).
trim_ps <- function(e, ps_trim) {
  if (!is.null(ps_trim)) {
    if (length(ps_trim) == 1L) ps_trim <- c(ps_trim, 1 - ps_trim)
    e <- pmin(pmax(e, min(ps_trim)), max(ps_trim))
  }
  if (any(e < 1e-6 | e > 1 - 1e-6)) {
    warning("Extreme propensity scores detected (near 0 or 1); ",
            "estimates may be unstable. Consider setting `ps_trim`.",
            call. = FALSE)
  }
  e
}

# Random fold assignment; returns NULL when folds <= 1 (no cross-fitting).
make_folds <- function(n, folds) {
  if (folds <= 1L) return(NULL)
  if (folds > n) {
    stop("`folds` cannot exceed the number of rows in `data`.", call. = FALSE)
  }
  sample(rep_len(seq_len(folds), n))
}

# Cross-fitted propensity scores e(x) = P(T = 1 | X = x).
# With fold_id = NULL the model is fit and evaluated in-sample.
crossfit_ps <- function(data, treatment, covariates, learner, fold_id,
                        ps_trim = NULL) {
  df <- data[, c(covariates, treatment), drop = FALSE]
  n <- nrow(df)
  e_hat <- numeric(n)
  if (is.null(fold_id)) {
    model <- fit_propensity(df, treatment, learner)
    e_hat <- predict_propensity(model, df[, covariates, drop = FALSE])
  } else {
    for (k in unique(fold_id)) {
      test <- fold_id == k
      model <- fit_propensity(df[!test, , drop = FALSE], treatment, learner)
      e_hat[test] <- predict_propensity(model, df[test, covariates, drop = FALSE])
    }
  }
  trim_ps(e_hat, ps_trim)
}

# Cross-fitted potential-outcome predictions mu0(x), mu1(x) from a single
# outcome model that includes the treatment as a feature (S-Learner style).
crossfit_outcome_po <- function(data, outcome, treatment, covariates, learner,
                                fold_id) {
  df <- data[, c(covariates, treatment, outcome), drop = FALSE]
  df[[treatment]] <- as.integer(df[[treatment]])
  n <- nrow(df)
  mu0 <- numeric(n)
  mu1 <- numeric(n)
  predict_po <- function(model, newdata) {
    nd <- newdata
    nd[[outcome]] <- NULL
    nd[[treatment]] <- 0L
    p0 <- predict_regr(model, nd)
    nd[[treatment]] <- 1L
    p1 <- predict_regr(model, nd)
    list(mu0 = p0, mu1 = p1)
  }
  if (is.null(fold_id)) {
    model <- fit_regr(df, outcome, learner)
    po <- predict_po(model, df)
    mu0 <- po$mu0
    mu1 <- po$mu1
  } else {
    for (k in unique(fold_id)) {
      test <- fold_id == k
      model <- fit_regr(df[!test, , drop = FALSE], outcome, learner)
      po <- predict_po(model, df[test, , drop = FALSE])
      mu0[test] <- po$mu0
      mu1[test] <- po$mu1
    }
  }
  list(mu0 = mu0, mu1 = mu1)
}

# Cross-fitted conditional mean m(x) = E[Y | X = x] (treatment excluded).
crossfit_outcome_mean <- function(data, outcome, covariates, learner, fold_id) {
  df <- data[, c(covariates, outcome), drop = FALSE]
  n <- nrow(df)
  m_hat <- numeric(n)
  if (is.null(fold_id)) {
    model <- fit_regr(df, outcome, learner)
    m_hat <- predict_regr(model, df[, covariates, drop = FALSE])
  } else {
    for (k in unique(fold_id)) {
      test <- fold_id == k
      model <- fit_regr(df[!test, , drop = FALSE], outcome, learner)
      m_hat[test] <- predict_regr(model, df[test, covariates, drop = FALSE])
    }
  }
  m_hat
}
