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
# File:         run_model.R
#
# Description:  Executes individual components of the replication package.
#               Each component consists of one or more scripts executed
#               sequentially in separate local environments.
#
# Input:        config/settings.R
#               data/information_set/selected_variables_50.csv
#               Scripts associated with the selected component.
#
# Output:       Model estimates, nowcasts, Monte Carlo results, diagnostics,
#               and publication tables and figures, depending on the
#               selected component.
#
# Usage:        Rscript run_model.R COMPONENT
#
# Note:         All paths are relative to the replication package root.
#               The information set is fixed throughout the replication.
#               No predictor screening or variable selection is performed.
#
# ==============================================================================


# ==============================================================================
# 1. COMPONENT SELECTION
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Check that exactly one component has been specified

if (length(args) != 1L) {
  stop("Usage: Rscript run_model.R COMPONENT")
}

component <- args[[1L]]


# ==============================================================================
# 2. PRELIMINARY CHECKS
# ==============================================================================

# Check the working directory

if (!file.exists("config/settings.R")) {
  stop("Run from the package root")
}

# Check the availability of the fixed information set

if (!file.exists("data/information_set/selected_variables_50.csv")) {
  stop("Missing selected_variables_50.csv")
}


# ==============================================================================
# 3. GRAPHICAL ENVIRONMENT
# ==============================================================================

# Prevent the automatic creation of Rplots.pdf when running scripts
# in non-interactive R sessions.
#
# Explicitly saved figures are unaffected, while interactive graphical
# output remains available when working in RStudio.

if (!interactive()) {
  options(
    device = function(...) {
      grDevices::pdf(file = NULL, ...)
    }
  )
}


# ==============================================================================
# 4. REPLICATION COMPONENTS
# ==============================================================================

components <- list(
  
  # --------------------------------------------------------------------------
  # 4.1. Preliminary data checks
  # --------------------------------------------------------------------------
  
  check = "scripts/00_data/01_check_inputs.R",
  
  
  # --------------------------------------------------------------------------
  # 4.2. Matrix MF-TPRF
  # --------------------------------------------------------------------------
  
  # Baseline specification: aggregate EA GDP proxy
  
  matrix = c(
    "scripts/01_empirical/01_matrix.R",
    "scripts/01_empirical/02_matrix_results.R"
  ),
  
  # Alternative specification: country-specific GDP proxies
  
  matrix_unit = c(
    "scripts/01_empirical/01_matrix.R",
    "scripts/01_empirical/02_matrix_results.R"
  ),
  
  
  # --------------------------------------------------------------------------
  # 4.3. Vector MF-TPRF
  # --------------------------------------------------------------------------
  
  # Country-specific specification
  
  vecc = "scripts/01_empirical/03_vector_country.R",
  
  # Pooled specification
  
  vecp = "scripts/01_empirical/04_vector_pooled.R",
  
  
  # --------------------------------------------------------------------------
  # 4.4. Vector Dynamic Factor Model
  # --------------------------------------------------------------------------
  
  dfm = "scripts/01_empirical/05_dfm.R",
  
  
  # --------------------------------------------------------------------------
  # 4.5. Euro Area analysis
  # --------------------------------------------------------------------------
  
  ea = "scripts/01_empirical/06_vector_ea.R",
  
  
  # --------------------------------------------------------------------------
  # 4.6. Monte Carlo simulations
  # --------------------------------------------------------------------------
  
  # Forecasting performance
  
  forecast_mc = "scripts/02_monte_carlo/01_forecasting.R",
  
  # Loading structure
  
  loading = "scripts/02_monte_carlo/02_loading_maps.R",
  
  # Asymptotic normality of the Pass 1 estimator
  
  an = "scripts/02_monte_carlo/03_an_pass1.R",
  
  
  # --------------------------------------------------------------------------
  # 4.7. Empirical loading structure and Wald tests
  # --------------------------------------------------------------------------
  
  paper_loading = "scripts/03_outputs/02_loading_wald.R",
  
  
  # --------------------------------------------------------------------------
  # 4.8. Paper tables and figures
  # --------------------------------------------------------------------------
  
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
    "scripts/03_outputs/10_ea_nowcasts.R",
    "scripts/03_outputs/11_computational_cost.R"
  ),
  
  
  # --------------------------------------------------------------------------
  # 4.9. Asymptotic normality diagnostics
  # --------------------------------------------------------------------------
  
  an_diagnostics = "scripts/04_diagnostics/01_an_diagnostics.R",

  # Generate the computational-cost table from saved timings only.
  computational_cost = "scripts/03_outputs/11_computational_cost.R",
  
  
  # --------------------------------------------------------------------------
  # 4.10. Final publication outputs
  # --------------------------------------------------------------------------
  
  finalize = "scripts/00_data/03_finalize_outputs.R",
  
  verify = "scripts/00_data/02_validate_outputs.R"
  
)


# ==============================================================================
# 5. COMPONENT VALIDATION
# ==============================================================================

valid_components <- names(components)

# Check whether the selected component is available

if (!component %in% valid_components) {
  stop(
    "Unknown component: ", component,
    ". Valid components: ",
    paste(valid_components, collapse = ", ")
  )
}


# ==============================================================================
# 6. MATRIX MF-TPRF PROXY SPECIFICATION
# ==============================================================================

# Ensure that the estimation and results scripts use the same
# proxy specification.

if (component == "matrix") {
  
  Sys.setenv(MF_TPRF_PROXY_MODE = "scalar")
  
}

if (component == "matrix_unit") {
  
  Sys.setenv(MF_TPRF_PROXY_MODE = "multivariate")
  
}


# ==============================================================================
# 7. COMPONENT EXECUTION
# ==============================================================================

# Retrieve the scripts associated with the selected component

selected_scripts <- components[[component]]

# Execute each script sequentially in a separate local environment

for (script in selected_scripts) {
  
  # Check the availability of the script
  
  if (!file.exists(script)) {
    stop("Missing script: ", script)
  }
  
  message("Running: ", script)
  
  # Create an isolated environment for the current script
  
  script_environment <- new.env(parent = globalenv())
  
  # Execute the script
  
  source(
    file  = script,
    local = script_environment,
    echo  = FALSE
  )
  
}


# ==============================================================================
# 8. END OF COMPONENT EXECUTION
# ==============================================================================

message("Completed: ", component)