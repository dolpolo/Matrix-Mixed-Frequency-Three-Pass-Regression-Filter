# ==============================================================================
# Appendix | Summary measures of relative forecast performance
# Run from the replication-package root:
#   source("scripts/03_outputs/09_forecast_wins.R")
# Uses saved baseline Matrix, VEC-P, VEC-C and DFM results. No re-estimation.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(zoo)
  library(sandwich)
})
source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/tables/forecast_wins.R")

cfg <- replication_config
summaries <- fw_read_summaries(size = "fixed50", selection = "fixed")

common_eval <- fw_prepare_evaluation(
  summaries = summaries,
  countries = cfg$countries,
  start_eval = cfg$start_eval,
  end_eval = cfg$end_eval
)
message("Common evaluable country-quarter-vintage forecasts: ", nrow(common_eval))
message("Before matching (finite observations per model): ",
        paste(names(attr(common_eval, "valid_forecasts_by_model")),
              unlist(attr(common_eval, "valid_forecasts_by_model")),
              sep = "=", collapse = "; "))

win_tables <- fw_save_all(
  evaluation = common_eval,
  output_dir = "results/paper/appendix/tables"
)

print(win_tables$country_rmsfe_share)
print(win_tables$forecast_summary)
message("Forecast win-share appendix tables completed.")
