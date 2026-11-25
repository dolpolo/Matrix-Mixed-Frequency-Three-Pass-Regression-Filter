# Matrix MF-TPRF publication replication orchestrator.
# Default: empirical models + loading structure + paper outputs + mandatory AN plots.
# Forecasting Monte Carlo is NEVER executed by default.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0L) "full" else args[[1L]]
valid <- c("full", "check", "empirical", "outputs", "loading", "an", "forecast_mc")
if (length(args) > 1L || !mode %in% valid)
  stop("Usage: Rscript run_all.R [full|check|empirical|outputs|loading|an|forecast_mc]")
if (!file.exists("config/settings.R")) stop("Run from the package root")

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) rscript <- "Rscript"
run <- function(component) {
  message("\n========== ", toupper(component), " ==========")
  result <- system2(rscript, c("--vanilla", "run_model.R", component),
                    stdout = "", stderr = "")
  if (length(result) != 1L || is.na(result) || result != 0L)
    stop("Replication stopped: ", component, " (exit code ",
         paste(result, collapse = ","), ")")
  message("========== DONE: ", component, " ==========")
}
if (mode %in% c("full", "check", "empirical")) run("check")
if (mode %in% c("full", "empirical")) {
  # Both Matrix proxy designs are needed for the appendix robustness comparison.
  for (component in c("matrix", "matrix_unit", "vecc", "vecp", "dfm", "ea"))
    run(component)
}
if (mode == "loading") run("paper_loading")
if (mode == "an") run("an_diagnostics")
if (mode == "forecast_mc") run("forecast_mc")
if (mode %in% c("full", "outputs")) {
  run("paper")
  run("an_diagnostics")
  run("finalize")
  run("verify")
}
message("\nReplication completed: ", mode)
