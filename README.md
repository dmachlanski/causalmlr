# causalmlr

Causal Machine Learning in R, built on the [mlr3](https://mlr3.mlr-org.com/) ecosystem.

`causalmlr` estimates **average treatment effects (ATE)** and **conditional average treatment effects (CATE)** from observational data using machine learning. Any mlr3 regression or classification learner can be plugged in as a nuisance model — from `lrn("regr.lm")` to a fully tuned `AutoTuner`.

The package grew out of the labs of the [Essex Summer School in Social Science Data Analysis](https://essexsummerschool.com/) course *Machine Learning for Estimating Treatment Effects from Observational Data*.

## Features

- **ATE estimators**: naive difference in means (`ate_naive()`), inverse propensity weighting (`ate_ipw()`), doubly robust / AIPW (`ate_dr()`), double machine learning (`ate_dml()`) — all with standard errors, confidence intervals, optional cross-fitting and propensity trimming.
- **CATE meta-learners**: `s_learner()`, `t_learner()`, `x_learner()`, `dr_learner()` (doubly robust), `r_learner()` (R-Loss / partially linear DML), with optional bootstrap confidence intervals on predicted CATEs (`cate_ci()`).
- **Evaluation**: absolute ATE error (`eps_ate()`), PEHE (`pehe()`), and the **R-Loss** (`r_loss()` + `rloss_nuisance()`) for model selection and hyperparameter tuning when no ground truth is available.
- **Benchmark datasets**: `sodium`, `synth_train`/`synth_test`, `ihdp_train`/`ihdp_test`, `jobs_train`/`jobs_test`, plus generic ML practice data (`abalone`, `diabetes`, `housing`, `spirals`).

## Installation

The package is not on CRAN yet. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dmachlanski/causalmlr")
```

## Quick start

### ATE estimation

```r
library(causalmlr)
library(mlr3)

data(sodium)  # true ATE = 1.05, confounded by age

ate_naive(sodium, outcome = "bp", treatment = "sodium")  # biased

est <- ate_dr(sodium, outcome = "bp", treatment = "sodium",
              outcome_learner = lrn("regr.rpart"),
              ps_learner = lrn("classif.rpart"),
              folds = 5)
est
#> ATE estimate - Doubly robust (AIPW)
#>   Estimate:  ~1.05
#>   ...
eps_ate(1.05, est)
```

### CATE estimation

```r
data(synth_train)
data(synth_test)

m <- x_learner(synth_train, outcome = "y", treatment = "t",
               learner = lrn("regr.rpart"),
               ps_learner = lrn("classif.rpart"))

tau_hat <- predict(m, synth_test)   # individual-level effects
ate(m, synth_test)                  # their average
pehe(synth_test$tau, tau_hat)       # evaluation against ground truth

# Pointwise bootstrap confidence intervals for the predicted CATEs
cate_ci(m, synth_test, train_data = synth_train, n_boot = 200)
```

### Causal model selection without ground truth

```r
nuis <- rloss_nuisance(valid_data, outcome = "y", treatment = "t",
                       outcome_learner = lrn("regr.ranger"),
                       ps_learner = lrn("classif.ranger"))

r_loss(nuis, predict(candidate_model, valid_data))  # lower is better
```

See the vignette (`vignette("causalmlr")`) for a complete walkthrough, including both causal hyperparameter tuning strategies (direct nuisance tuning via `mlr3tuning` and indirect tuning via the R-Loss).

## Development setup

After cloning the repository, two one-time steps are needed before building (both require R with `devtools` installed):

```r
# 1. Generate the lazy-loaded datasets in data/ from the CSVs in data-raw/
#    (commit the resulting data/*.rda files so install_github() works)
source("data-raw/prepare_data.R")

# 2. Generate man/ pages from the roxygen comments
devtools::document()

# Then the usual workflow:
devtools::test()
devtools::check()
```

## License

LGPL (>= 2.1). See [LICENSE](LICENSE).
