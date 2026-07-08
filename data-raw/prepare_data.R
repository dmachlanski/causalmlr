# Generates the lazy-loaded datasets in data/ from the CSVs in data-raw/.
#
# Run once from the package root (requires only base R):
#   Rscript data-raw/prepare_data.R
#
# Re-run whenever a CSV in data-raw/ changes.

read_raw <- function(name) {
  utils::read.csv(file.path("data-raw", paste0(name, ".csv")))
}

save_data <- function(...) {
  names <- vapply(substitute(list(...))[-1L], deparse, character(1L))
  dir.create("data", showWarnings = FALSE)
  for (i in seq_along(names)) {
    obj <- list(...)[[i]]
    assign(names[i], obj)
    save(list = names[i], file = file.path("data", paste0(names[i], ".rda")),
         compress = "xz", version = 2)
    message("Wrote data/", names[i], ".rda (", nrow(obj), " rows)")
  }
}

# --- Causal datasets ---------------------------------------------------------

sodium <- read_raw("sodium")
sodium$sodium <- as.integer(sodium$sodium)

synth_train <- read_raw("synth_train")
synth_train$t <- as.integer(synth_train$t)
synth_test <- read_raw("synth_test")

ihdp_train <- read_raw("ihdp_train")
ihdp_train$t <- as.integer(ihdp_train$t)
ihdp_test <- read_raw("ihdp_test")
ihdp_test$t <- as.integer(ihdp_test$t)

jobs_train <- read_raw("jobs_train")
jobs_train$t <- as.integer(jobs_train$t)
jobs_train$y <- as.integer(jobs_train$y)
jobs_train$e <- as.integer(jobs_train$e)
jobs_test <- read_raw("jobs_test")
jobs_test$t <- as.integer(jobs_test$t)
jobs_test$y <- as.integer(jobs_test$y)
jobs_test$e <- as.integer(jobs_test$e)

# --- Supervised ML practice datasets -----------------------------------------

abalone <- read_raw("abalone")
abalone$sex <- factor(abalone$sex)

diabetes <- read_raw("diabetes")
diabetes$diabetes <- factor(diabetes$diabetes)

housing <- read_raw("housing")

spirals <- read_raw("spirals")
spirals$y <- factor(as.integer(spirals$y))

# --- Save --------------------------------------------------------------------

save_data(sodium, synth_train, synth_test, ihdp_train, ihdp_test,
          jobs_train, jobs_test, abalone, diabetes, housing, spirals)
