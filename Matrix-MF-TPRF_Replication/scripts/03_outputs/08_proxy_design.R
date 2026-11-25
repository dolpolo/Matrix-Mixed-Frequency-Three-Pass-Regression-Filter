# ==============================================================================
# Appendix: unit-target proxy design versus aggregate-target proxy design and DFM
# Run from the replication-package root:
#   source("scripts/03_outputs/08_proxy_design.R")
# Reads already saved current-information-set results; no estimation.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(sandwich)
})

source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/tables/unit_proxy_comparison.R")

cfg <- replication_config
results <- load_unit_comparison_results(size = "fixed50", selection = "fixed")
comparison <- make_unit_comparison_data(
  results = results,
  countries = cfg$countries,
  start_eval = cfg$start_eval,
  end_eval = cfg$end_eval
)

save_unit_comparison(
  result = comparison,
  output_dir = "results/paper/appendix/tables",
  stem = "table_unit_target_comparison"
)

print(comparison$table[, c("country", "period", "type", "n",
  "rmsfe_unit", "ratio_scalar", "ratio_dfm")], row.names = FALSE)
message("Unit-target versus aggregate-target proxy design and DFM comparison completed.")
