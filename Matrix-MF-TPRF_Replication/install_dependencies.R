# Author: Davide Delfino
# Institution: Aarhus University
# Paper: Matrix Three-Pass Regression Filter: A Nowcasting Framework for Mixed-Frequency Multidimensional Data
# Last updated: 2026-11-24
# Input: installed R environment and CRAN access.
# Output: installs missing packages required by the replication package.

packages <- c(
  "tidyverse", "readxl", "zoo", "lubridate", "abind", "MASS", "tseries",
  "fBasics", "vars", "glmnet", "plsdof", "sandwich", "lmtest", "car",
  "conflicted", "future", "future.apply", "progressr", "forcats",
  "writexl", "xtable", "patchwork", "scales", "Matrix", "png"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) utils::install.packages(missing, repos = "https://cloud.r-project.org")
if (length(packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]))
  stop("Some packages could not be installed.")
message("Dependencies are available.")
