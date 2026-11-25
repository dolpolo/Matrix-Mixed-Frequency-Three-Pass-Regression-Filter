# One component at a time. All paths relative to the replication root.
# No predictor screening is rerun; each stage is independently executable.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript run_model.R COMPONENT")
component <- args[[1L]]
if (!file.exists("config/settings.R")) stop("Run from the package root")
if (!file.exists("data/information_set/selected_variables_50.csv"))
  stop("Missing selected_variables_50.csv")

components <- list(
  check = "scripts/00_data/01_check_inputs.R",
  matrix = c("scripts/01_empirical/01_matrix.R",
             "scripts/01_empirical/02_matrix_results.R"),
  matrix_unit = c("scripts/01_empirical/01_matrix.R",
                  "scripts/01_empirical/02_matrix_results.R"),
  vecc = "scripts/01_empirical/03_vector_country.R",
  vecp = "scripts/01_empirical/04_vector_pooled.R",
  dfm = "scripts/01_empirical/05_dfm.R",
  ea = "scripts/01_empirical/06_vector_ea.R",
  loading = "scripts/02_monte_carlo/02_loading_maps.R",
  an = "scripts/02_monte_carlo/03_an_pass1.R",
  forecast_mc = "scripts/02_monte_carlo/01_forecasting.R",
  paper_loading = "scripts/03_outputs/02_loading_wald.R",
  paper = c(
    "scripts/03_outputs/01_empirical_comparison.R",
    "scripts/03_outputs/02_loading_wald.R",
    "scripts/03_outputs/03_predictor_panel.R",
    "scripts/03_outputs/04_matrix_factor.R",
    "scripts/03_outputs/05_ea_factor.R",
    "scripts/03_outputs/06_dfm_comparison.R",
    "scripts/03_outputs/07_country_nowcasts.R",
    "scripts/03_outputs/08_proxy_design.R",
    "scripts/03_outputs/09_forecast_wins.R",
    "scripts/03_outputs/10_ea_nowcasts.R"
  ),
  an_diagnostics = "scripts/04_diagnostics/01_an_diagnostics.R",
  finalize = "scripts/00_data/03_finalize_outputs.R",
  verify = "scripts/00_data/02_validate_outputs.R"
)
if (!component %in% names(components))
  stop("Unknown component: ", component, ". Valid: ", paste(names(components), collapse=", "))
# Guarantee the Matrix results script selects the same proxy as the estimator.
if (component %in% c("matrix", "matrix_unit")) {
  Sys.setenv(MF_TPRF_PROXY_MODE=if (component == "matrix") "scalar" else "multivariate")
}
for (script in components[[component]]) {
  if (!file.exists(script)) stop("Missing script: ", script)
  message("Running: ", script)
  # Isolated scope per script, as in the original runner.
  source(script, local = new.env(parent=globalenv()), echo=FALSE)
}
message("Completed: ", component)
