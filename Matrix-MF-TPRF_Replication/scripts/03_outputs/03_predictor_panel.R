# ==============================================================================
# Output main | Predictor panel: composition and unavailability
# Run from the replication package root:
#   source("scripts/03_outputs/03_predictor_panel.R")
# Called automatically by run_model.R paper / run_all.R outputs in the package.
# Reads the dataset; does not estimate models or select predictors.
# ==============================================================================

if (!file.exists("config/settings.R") ||
    !file.exists("R/tables/predictor_panel_tables.R")) {
  stop("Run from the replication package root and copy the new R/tables module.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(lubridate)
  library(zoo)
  library(abind)
})

source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/matrix/matrix.mf.tprf.prep.R")
source("R/common/fixed_selection.R")
source("R/matrix/matrix.mf.tprf.now.R")
source("R/tables/predictor_panel_tables.R")

# The longtable has the original country columns, description, frequency,
# transformation and release-delay metadata, plus a final count of unavailable
# country-specific series.
composition <- make_information_set_table(
  countries = replication_config$countries,
  csv_path = "data/information_set/selected_variables_50.csv",
  transformed_dir = "data/transformed",
  output_file = "results/paper/appendix/tables/information_set.tex"
)

# Shares of unavailable observations before matrix imputation, averaged across
# M1/M2/M3 within each period. One row per period, three numeric columns.
missingness <- make_predictor_unavailability_table(
  config = replication_config,
  output_file = "results/paper/appendix/tables/realtime_missingness.tex"
)

cat("\nDataset tables created:\n",
    "  ", composition$path, "\n",
    "  ", missingness$path, "\n", sep = "")
cat("\nUnavailable series by country:\n")
print(composition$missing_by_country)
cat("\nAverage unavailable predictor shares (%):\n")
print(missingness$summary, row.names = FALSE)
