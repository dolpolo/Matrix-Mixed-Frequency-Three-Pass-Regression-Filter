# ==============================================================================
# Fixed-50 empirical comparison table
# Author: Davide Delfino
# Institution: Aarhus University
# Paper: Matrix Three-Pass Regression Filter: A Nowcasting Framework for Mixed-Frequency Multidimensional Data
# Last updated: 2026-11-24
# Input: saved Matrix, VEC-C, VEC-P, and DFM national summaries
# Output: results/paper/main/tables/table_2_empirical.tex
# ==============================================================================
# The published table layout and Diebold-Mariano procedure are reused from
# the author's original matrix.mf.tprf.utils.R.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(lubridate)
  library(sandwich); library(zoo)
})
source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
cfg <- replication_config
out_t <- "results/paper/main/tables"
dir.create(out_t, recursive = TRUE, showWarnings = FALSE)

load_summary <- function(path, model, dfm = FALSE) {
  file <- if (dfm) {
    find_dfm_summary_file(path = path, Size = "fixed50", sel = "fixed")
  } else {
    find_result_file(path = path, model = model, stage = "summary",
                     Size = "fixed50", sel = "fixed")
  }
  if (length(file) != 1L || !file.exists(file)) stop("Missing model summary for ", model)
  readRDS(file)
}

mat <- load_summary("results/empirical/matrix", "matrix_scalar")
vecc <- load_summary("results/empirical/vector", "vector")
vecp <- load_summary("results/empirical/pooled", "vectensor")
dfm <- load_summary("results/empirical/dfm/q0", "vector_dfm", dfm = TRUE)
if (!identical(as.character(mat$proxy_mode), "scalar")) stop("Expected scalar EA proxy design.")

# Same pseudo-real-time squared errors and one-sided HAC DM tests as in the paper.
eval_mat <- build_eval_df(mat$df_rt_all, mat$df_yq_eval_all)
eval_vp  <- build_eval_df(vecp$df_rt_all, mat$df_yq_eval_all)
eval_vc  <- build_eval_df(vecc$df_rt_all, mat$df_yq_eval_all)
eval_dfm <- build_eval_df(dfm$df_rt_all, mat$df_yq_eval_all)
dm_vp <- build_dm_wide(eval_mat, eval_vp, "VecP", alternative = "less")
dm_vc <- build_dm_wide(eval_mat, eval_vc, "VecC", alternative = "less")
dm_dfm <- build_dm_wide(eval_mat, eval_dfm, "DFM", alternative = "less")

params <- list(start_eval = cfg$start_eval, end_eval = cfg$end_eval,
               covid_start = cfg$covid_start, covid_end = cfg$covid_end)
tex <- build_final_large_style_latex(
  summary_matrix = mat, summary_vector = vecc,
  summary_vectensor = vecp, summary_dfm = dfm,
  params = params, dm_vecp_wide = dm_vp$wide,
  dm_vecc_wide = dm_vc$wide, dm_dfm_wide = dm_dfm$wide,
  Size = "fixed50", sel = "fixed"
)
writeLines(tex, file.path(out_t, "table_2_empirical.tex"))
dir.create("results/empirical/comparison", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(summaries = list(matrix = mat, vecp = vecp, vecc = vecc, dfm = dfm),
             dm = list(vecp = dm_vp, vecc = dm_vc, dfm = dm_dfm),
             evaluation = list(matrix = eval_mat, vecp = eval_vp,
                               vecc = eval_vc, dfm = eval_dfm)),
        "results/empirical/comparison/empirical_comparison.rds")
cat("Saved the original empirical comparison table.\n")
