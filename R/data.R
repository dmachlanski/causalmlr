# Documentation for the datasets shipped with the package.
# The .rda files are generated from the CSVs in data-raw/ by running
# data-raw/prepare_data.R (see that script for details).

#' Sodium intake and blood pressure (synthetic)
#'
#' A simplified simulation of the effect of sodium intake on blood pressure,
#' proposed by Luque-Fernandez et al. (2019). Age confounds both the
#' treatment (high sodium intake) and the outcome (blood pressure). The true
#' ATE used in the data generating process is **1.05**.
#'
#' @format A data frame with 10,000 rows and 3 columns:
#' \describe{
#'   \item{age}{Age in years (confounder).}
#'   \item{sodium}{Binary treatment: high sodium intake (1) or not (0).}
#'   \item{bp}{Systolic blood pressure (outcome).}
#' }
#' @references
#' Luque-Fernandez, M. A., Schomaker, M., Redondo-Sanchez, D., Jose
#' Sanchez Perez, M., Vaidya, A., & Schnitzer, M. E. (2019). Educational
#' Note: Paradoxical collider effect in the analysis of non-communicable
#' disease epidemiological data. *International Journal of Epidemiology*,
#' 48(2), 640-653.
#' @examples
#' data(sodium)
#' ate_naive(sodium, outcome = "bp", treatment = "sodium")
"sodium"

#' Synthetic CATE benchmark (training set)
#'
#' Fully synthetic data generated with the data generating process of
#' Künzel et al. (2019), with five Gaussian covariates, a binary treatment
#' and a continuous outcome. The true CATE is a step function of one of the
#' covariates. See [synth_test] for the accompanying test set with ground
#' truth effects.
#'
#' @format A data frame with 1,000 rows and 7 columns:
#' \describe{
#'   \item{x0, x1, x2, x3, x4}{Continuous covariates.}
#'   \item{t}{Binary treatment indicator (0/1).}
#'   \item{y}{Continuous outcome.}
#' }
#' @references
#' Künzel, S. R., Sekhon, J. S., Bickel, P. J., & Yu, B. (2019).
#' Metalearners for estimating heterogeneous treatment effects using machine
#' learning. *Proceedings of the National Academy of Sciences*, 116(10),
#' 4156-4165.
#' @examples
#' data(synth_train)
#' head(synth_train)
"synth_train"

#' Synthetic CATE benchmark (test set)
#'
#' Test covariates and ground-truth treatment effects accompanying
#' [synth_train]. Contains only the covariates and the true CATE (`tau`),
#' mimicking the deployment setting where a fitted CATE model predicts
#' effects for new individuals.
#'
#' @format A data frame with 250 rows and 6 columns:
#' \describe{
#'   \item{x0, x1, x2, x3, x4}{Continuous covariates.}
#'   \item{tau}{True conditional average treatment effect.}
#' }
#' @inherit synth_train references
#' @examples
#' data(synth_test)
#' head(synth_test)
"synth_test"

#' IHDP semi-synthetic benchmark (training set)
#'
#' The Infant Health and Development Program (IHDP) benchmark introduced by
#' Hill (2011): real covariates describing children and their mothers from a
#' randomised experiment, combined with a simulated outcome, which provides
#' ground-truth treatment effects. The treatment is intensive, high-quality
#' childcare and home visits; the (simulated) outcome mimics future
#' cognitive test scores.
#'
#' @format A data frame with 672 rows and 28 columns:
#' \describe{
#'   \item{x0-x5}{Continuous covariates (standardised).}
#'   \item{x6-x24}{Binary covariates.}
#'   \item{t}{Binary treatment indicator (0/1).}
#'   \item{y}{Continuous (simulated) outcome.}
#'   \item{tau}{True conditional average treatment effect.}
#' }
#' @references
#' Hill, J. L. (2011). Bayesian nonparametric modeling for causal inference.
#' *Journal of Computational and Graphical Statistics*, 20(1), 217-240.
#' @examples
#' data(ihdp_train)
#' head(ihdp_train)
"ihdp_train"

#' IHDP semi-synthetic benchmark (test set)
#'
#' Held-out portion of the IHDP benchmark; see [ihdp_train] for details.
#'
#' @format A data frame with 75 rows and 28 columns; same columns as
#'   [ihdp_train].
#' @inherit ihdp_train references
#' @examples
#' data(ihdp_test)
#' head(ihdp_test)
"ihdp_test"

#' Jobs benchmark (training set)
#'
#' The Jobs benchmark (LaLonde, 1986; composition of Shalit et al., 2017)
#' combines a randomised experiment (the National Supported Work programme)
#' with observational controls. The treatment is participation in a job
#' training programme and the outcome is employment status. The `e` column
#' flags units belonging to the randomised subsample, which enables
#' evaluation on real data via the experimental subset.
#'
#' @format A data frame with 2,891 rows and 20 columns:
#' \describe{
#'   \item{x0-x16}{Covariates (standardised continuous and binary).}
#'   \item{t}{Binary treatment indicator: job training (0/1).}
#'   \item{y}{Binary outcome: employment (0/1).}
#'   \item{e}{Indicator of membership in the randomised experimental
#'     subsample (0/1).}
#' }
#' @references
#' LaLonde, R. J. (1986). Evaluating the econometric evaluations of training
#' programs with experimental data. *The American Economic Review*, 76(4),
#' 604-620.
#'
#' Shalit, U., Johansson, F. D., & Sontag, D. (2017). Estimating individual
#' treatment effect: generalization bounds and algorithms. *Proceedings of
#' the 34th International Conference on Machine Learning*.
#' @examples
#' data(jobs_train)
#' head(jobs_train)
"jobs_train"

#' Jobs benchmark (test set)
#'
#' Held-out portion of the Jobs benchmark; see [jobs_train] for details.
#'
#' @format A data frame with 321 rows and 20 columns; same columns as
#'   [jobs_train].
#' @inherit jobs_train references
#' @examples
#' data(jobs_test)
#' head(jobs_test)
"jobs_test"

#' Abalone (regression practice data)
#'
#' The classic UCI abalone dataset: predict the number of rings (a proxy for
#' age) of abalone from physical measurements. Included as standard
#' supervised-learning practice data; it has no causal structure.
#'
#' @format A data frame with 4,177 rows and 9 columns: `sex` (factor),
#'   seven physical measurements and the target `rings` (integer).
#' @source \url{https://archive.ics.uci.edu/dataset/1/abalone}
#' @examples
#' data(abalone)
#' head(abalone)
"abalone"

#' Pima Indians diabetes (classification practice data)
#'
#' Diagnostic measurements for predicting diabetes onset. Included as
#' standard supervised-learning practice data; it has no causal structure.
#'
#' @format A data frame with 768 rows and 9 columns: eight numeric
#'   diagnostic measurements and the target `diabetes` (factor, 0/1).
#' @examples
#' data(diabetes)
#' head(diabetes)
"diabetes"

#' Boston housing (regression practice data)
#'
#' The Boston housing dataset (Harrison & Rubinfeld, 1978): predict median
#' house value (`MEDV`) from neighbourhood characteristics. Included as
#' standard supervised-learning practice data. Note that this dataset has
#' known ethical concerns (the `B` column encodes a racial statistic) and is
#' provided for teaching purposes only.
#'
#' @format A data frame with 506 rows and 14 columns.
#' @examples
#' data(housing)
#' head(housing)
"housing"

#' Spirals (classification practice data)
#'
#' A synthetic two-dimensional binary classification problem where the two
#' classes form interleaved spirals - a simple example of data that linear
#' models cannot separate. Included as supervised-learning practice data.
#'
#' @format A data frame with 2,000 rows and 3 columns: coordinates `x0`,
#'   `x1` and the class label `y` (factor, 0/1).
#' @examples
#' data(spirals)
#' head(spirals)
"spirals"
