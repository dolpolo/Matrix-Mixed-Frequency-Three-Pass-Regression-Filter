# ==============================================================================
# MATRIX THREE-PASS REGRESSION FILTER:
# A NOWCASTING FRAMEWORK FOR MIXED-FREQUENCY MULTIDIMENSIONAL DATA
# ==============================================================================
#
# Authors:      Matteo Barigozzi
#               Davide Delfino
#               Massimiliano Marcellino
#
# Institutions: University of Bologna
#               Aarhus University
#               Bocconi University
#
# Contact:      davide.delfino@econ.au.dk
#
# File:         run_all.R
#
# Description:  Main script of the replication package. It coordinates the
#               execution of the empirical models, loading structure analysis,
#               asymptotic normality diagnostics, and publication outputs.
#
# Input:        config/settings.R
#               run_model.R
#               Data and scripts contained in the replication package.
#
# Output:       Empirical nowcasts, estimation results, diagnostic results,
#               and tables and figures required for the paper.
#
# Usage:        Rscript run_all.R [mode]
#
# Execution modes:
#
#   full         Complete paper replication (default).
#   check        Check the replication package requirements.
#   empirical    Estimate all empirical models.
#   outputs      Generate and verify the paper outputs.
#   loading      Perform the loading structure analysis.
#   an           Run the asymptotic normality diagnostics.
#   forecast_mc  Run the forecasting Monte Carlo simulations.
#   benchmark    Measure Matrix, VEC-P and DFM in an isolated workspace.
#
# Note:         Forecasting Monte Carlo simulations are computationally
#               intensive and are excluded from the default replication.
#               They can be executed separately using the forecast_mc mode.
#
# ==============================================================================


# ==============================================================================
# 1. REPLICATION MODE
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Select the execution mode (default: full replication)

if (length(args) == 0L) {
  mode <- "full"
} else {
  mode <- args[[1L]]
}

# Available execution modes

valid_modes <- c(
  "full",
  "check",
  "empirical",
  "outputs",
  "loading",
  "an",
  "forecast_mc",
  "benchmark"
)

# Check the selected execution mode

if (length(args) > 1L || !mode %in% valid_modes) {
  stop(
    "Usage: Rscript run_all.R ",
    "[full|check|empirical|outputs|loading|an|forecast_mc|benchmark]"
  )
}

# Check the working directory

if (!file.exists("config/settings.R")) {
  stop("Run from the package root")
}


# ==============================================================================
# 2. R EXECUTION ENVIRONMENT
# ==============================================================================

# Load the timing utilities (no estimation or data download).
source("R/common/computational_cost.R")

# Locate the Rscript executable according to the operating system

if (.Platform$OS.type == "windows") {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
} else {
  rscript <- file.path(R.home("bin"), "Rscript")
}

# Use the system executable if the specified path is unavailable

if (!file.exists(rscript)) {
  rscript <- "Rscript"
}


# ==============================================================================
# 3. REPLICATION FUNCTION
# ==============================================================================

# Execute each component in a separate R session.
# Stop the replication if a component fails.

# Collect timings only during a complete three-model comparison.
timing_directory <- NULL
timing_workers <- NULL

if (mode %in% c("full", "empirical", "benchmark")) {
  timing_workers <- computational_worker_budget()
  timing_directory <- tempfile("mf_timing_")
  dir.create(timing_directory, recursive = TRUE)
  timing_directory <- normalizePath(timing_directory, winslash = "/")
}

run_component <- function(component) {
  
  message("\n========== ", toupper(component), " ==========")
  
  # Execute the selected component
  
  launch <- function() {
    system2(
      command = rscript,
      args    = c("--vanilla", "run_model.R", component),
      stdout  = "",
      stderr  = ""
    )
  }

  if (!is.null(timing_directory) && component %in% computational_components) {
    exit_code <- with_computational_environment(
      timing_directory, timing_workers, launch()
    )
  } else {
    exit_code <- launch()
  }
  
  # Check whether the execution was successful
  
  if (length(exit_code) != 1L || is.na(exit_code) || exit_code != 0L) {
    stop(
      "Replication stopped: ", component,
      " (exit code ", paste(exit_code, collapse = ","), ")"
    )
  }
  
  message("========== DONE: ", component, " ==========")
  
}


# ==============================================================================
# 4. PRELIMINARY CHECKS
# ==============================================================================

if (mode %in% c("full", "check", "empirical")) {
  
  run_component("check")
  
}


# ==============================================================================
# 5. EMPIRICAL ANALYSIS
# ==============================================================================

if (mode %in% c("full", "empirical")) {
  
  # --------------------------------------------------------------------------
  # 5.1. Matrix MF-TPRF
  # --------------------------------------------------------------------------
  
  run_component("matrix")
  
  # Alternative proxy specification (appendix robustness analysis)
  
  run_component("matrix_unit")
  
  
  # --------------------------------------------------------------------------
  # 5.2. Vector MF-TPRF
  # --------------------------------------------------------------------------
  
  run_component("vecc")
  
  run_component("vecp")
  
  
  # --------------------------------------------------------------------------
  # 5.3. Vector Dynamic Factor Model
  # --------------------------------------------------------------------------
  
  run_component("dfm")

  # All three timings were collected during the normal empirical run.
  collect_computational_records(timing_directory)
  
  
  # --------------------------------------------------------------------------
  # 5.4. Euro Area analysis
  # --------------------------------------------------------------------------
  
  run_component("ea")
  
}


# ==============================================================================
# 6. LOADING STRUCTURE ANALYSIS
# ==============================================================================

if (mode == "loading") {
  
  run_component("paper_loading")
  
}


# ==============================================================================
# 7. ASYMPTOTIC NORMALITY DIAGNOSTICS
# ==============================================================================

if (mode == "an") {
  
  run_component("an_diagnostics")
  
}


# ==============================================================================
# 8. FORECASTING MONTE CARLO SIMULATIONS
# ==============================================================================

# Executed only when explicitly requested.

if (mode == "forecast_mc") {
  
  run_component("forecast_mc")
  
}


# ==============================================================================
# 9. INDEPENDENT COMPUTATIONAL BENCHMARK
# ==============================================================================

# Re-estimate only the three timed models in a temporary copy of the inputs
# and source code. Existing empirical RDS files and figures are not touched.
run_isolated_benchmark <- function() {
  original_root <- getwd()
  workspace <- tempfile("mf_bench_")
  dir.create(workspace, recursive = TRUE)
  on.exit({
    setwd(original_root)
    unlink(workspace, recursive = TRUE)
  }, add = TRUE)

  inputs <- c("R", "scripts", "config", "data", "run_model.R")
  copied <- file.copy(inputs, workspace, recursive = TRUE, copy.mode = TRUE)
  if (!all(copied)) stop("Could not create the isolated benchmark workspace.")

  setwd(workspace)
  for (component in computational_components) {
    run_component(component)
  }
  setwd(original_root)

  collect_computational_records(timing_directory)
  run_component("computational_cost")
}

if (mode == "benchmark") {
  run_isolated_benchmark()
}


# ==============================================================================
# 10. PAPER TABLES AND FIGURES
# ==============================================================================

if (mode %in% c("full", "outputs")) {

  if (!file.exists(computational_records_file)) {
    stop("Computational timings are missing. Run: Rscript run_all.R benchmark")
  }
  
  # Generate the empirical tables and figures
  
  run_component("paper")
  
  # Generate the asymptotic normality diagnostics
  
  run_component("an_diagnostics")
  
  # Organize the final publication outputs
  
  run_component("finalize")
  
  # Verify the completeness of the replication outputs
  
  run_component("verify")
  
}


# ==============================================================================
# 11. END OF REPLICATION
# ==============================================================================

if (!is.null(timing_directory)) {
  unlink(timing_directory, recursive = TRUE)
}

message("\nReplication completed: ", mode)