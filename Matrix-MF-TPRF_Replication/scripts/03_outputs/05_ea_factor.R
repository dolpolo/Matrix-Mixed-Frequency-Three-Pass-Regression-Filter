# ==============================================================================
# Output main 2 | EA Matrix vs Vector factor comparison and YoY R2 table
# Run from the project root:
#   source("scripts/03_outputs/05_ea_factor.R")
# Reads pre-existing fit RDS files; does not re-estimate any model.
# ==============================================================================

if (!file.exists("R/matrix/matrix.mf.tprf.utils.R"))
  stop("Run this script from the replication package root.")

suppressPackageStartupMessages(library(ggplot2))
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/figures/ea_factor_comparison_page.R")

# The Size / sel strings below are the internal names of existing RDS outputs.
matrix_path <- find_result_file(
  path = "results/empirical/matrix", model = "matrix_scalar",
  stage = "fit", Size = "fixed50", sel = "fixed"
)
vector_path <- find_result_file(
  path = "results/empirical/vector/EA", model = "vector",
  stage = "fit", Size = "fixed50", sel = "fixed"
)

make_ea_factor_comparison_page(
  matrix_fit = readRDS(matrix_path),
  vector_fit = readRDS(vector_path),
  figure_dir = "results/paper/appendix/figures",
  table_dir = "results/paper/appendix/tables"
)
message("Matrix fit: ", matrix_path)
message("Vector EA fit: ", vector_path)
