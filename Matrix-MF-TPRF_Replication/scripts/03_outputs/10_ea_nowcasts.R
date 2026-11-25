# ==============================================================================
# EA GDP growth nowcasts | Matrix MF-TPRF vs VEC-C
# Run from the replication package root:
#   source("scripts/03_outputs/10_ea_nowcasts.R")
# Uses saved fixed-50 results only; no estimators are rerun.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(scales)
})

source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/figures/ea_nowcasts.R")

saved_ea <- load_ea_nowcast_results(size = "fixed50", selection = "fixed")
data_ea <- prepare_ea_nowcast_data(saved_ea, replication_config)
table_ea <- compute_ea_rmsfe_table(data_ea)
plot_ea <- make_ea_nowcast_plot(data_ea, replication_config)

save_ea_nowcast_outputs(data_ea, table_ea, plot_ea)

cat("\nEA RMSFE comparison (full precision):\n")
print(as.data.frame(table_ea$long), row.names = FALSE)
cat("\nEA evaluation coverage:\n")
print(as.data.frame(data_ea$coverage), row.names = FALSE)
message("EA comparison complete: one M1--M3 figure and one RMSFE table.")
