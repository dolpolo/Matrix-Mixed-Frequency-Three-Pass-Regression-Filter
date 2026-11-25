# ==============================================================================
# Country nowcasts | 8 national figures + unchanged post-Covid figure
# From the replication-package root:
#   source("scripts/03_outputs/07_country_nowcasts.R")
# Uses SAVED summaries; no model estimation.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("config/settings.R")
source("R/matrix/matrix.mf.tprf.utils.R")
source("R/figures/country_nowcasts.R")

cfg <- replication_config
saved <- load_country_nowcasts_summaries(size = "fixed50", selection = "fixed")
out_dir <- "results/paper/appendix/figures"

# 1. Full evaluation: individual M1/M2/M3 figure for EACH country.
for (country in c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")) {
  p <- make_country_nowcast(
    summaries = saved, country = country,
    start_date = cfg$start_eval, end_date = cfg$end_eval,
    covid_start = cfg$covid_start, covid_end = cfg$covid_end,
    model_linewidth = 0.78, gdp_linewidth = 1.10
  )
  save_country_nowcast(
    plot = p, country = country, output_dir = out_dir,
    width = 12.8, height = 4.85, dpi = 450
  )
}

# 2. Post-Covid: monthly M1-M2-M3 paths, one 4x2 plot for 8 countries.
post_start <- seq(as.Date(cfg$covid_end), by = "1 month", length.out = 2L)[2L]
post_end <- as.Date(cfg$end_eval)
post_data <- prepare_postcovid_country_data(
  summaries = saved,
  countries = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  start_date = post_start,
  end_date = post_end
)
p_post <- make_postcovid_country_plot(
  data = post_data,
  start_date = post_start,
  end_date = post_end,
  free_y = TRUE,
  show_monthly_points = TRUE,
  x_labels = "annual"  # Annual labels, vertical; light quarterly gridlines.
  # Change to "quarterly" to label Q1/Q2/Q3/Q4 instead.
)
save_postcovid_country_plot(
  plot = p_post,
  data = post_data,
  output_dir = out_dir,
  stem = "figure_08_postcovid_country_nowcasts_4x2",
  width = 8.1, height = 10.4, dpi = 450
)

message("Generated 8 full-evaluation country figures and 1 post-Covid figure.")
