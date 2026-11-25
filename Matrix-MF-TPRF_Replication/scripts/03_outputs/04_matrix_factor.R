# ==============================================================================
# Output main 1 | Matrix static factor, row/column loadings, and loading map
# Run from the project root:
#   source("scripts/03_outputs/04_matrix_factor.R")
# Does not estimate or modify the Matrix MF-TPRF model.
# ==============================================================================

if (!file.exists("R/matrix/matrix.mf.tprf.utils.R"))
  stop("Run this script from the replication package root.")

suppressPackageStartupMessages(library(ggplot2))
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/figures/matrix_static_factor_page.R")

# Keep the existing internal saved-result naming convention (Size / sel).
fit_path <- find_result_file(
  path = "results/empirical/matrix", model = "matrix_scalar",
  stage = "fit", Size = "fixed50", sel = "fixed"
)
res_fit <- readRDS(fit_path)

make_matrix_static_factor_page(
  fit = res_fit,
  figure_dir = "results/paper/appendix/figures",
  table_dir = "results/paper/appendix/tables"
)
message("Source fit: ", fit_path)
