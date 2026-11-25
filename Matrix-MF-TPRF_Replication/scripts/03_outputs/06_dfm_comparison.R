# ==============================================================================
# Appendix DFM RMSFE and ratios vs VEC-P and VEC-C | paper output main
# Run from the replication-package root, after empirical models have been run.
# Input: saved q0 DFM, VEC-P and VEC-C summary RDS; no estimation.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(sandwich)
})

if (!file.exists("config/settings.R")) {
  stop("Please run this script from the replication-package root.")
}
source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/tables/dfm_relative_rmsfe_table.R")

cfg <- replication_config
summaries <- load_dfm_benchmark_summaries()
comparison <- prepare_dfm_benchmark_comparison(summaries)

params <- list(
  start_eval = cfg$start_eval,
  end_eval = cfg$end_eval,
  covid_start = cfg$covid_start,
  covid_end = cfg$covid_end
)
outputs <- save_dfm_benchmark_comparison(
  comparison, params,
  out_dir = "results/paper/appendix/tables",
  stem = "table_dfm_relative_rmsfe"
)

message("Saved appendix DFM benchmarking outputs:")
for (p in outputs) message("  ", normalizePath(p, mustWork = TRUE))
message("RMSFEs and one-sided HAC Diebold--Mariano tests come from saved results.")
