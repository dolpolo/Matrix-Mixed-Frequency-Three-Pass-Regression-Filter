path_main <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Author: Davide Delfino
# Input: data/, config/settings.R, and functions in R/.
# Output: results/; see README.md.
# Original: MonteCarlo_simulation/matrix.loading.main.R

# =============================================================================
# MONTE CARLO -- ONE-DRAW LOADING STRUCTURE COMPARISON (NO NOWCASTING)
# =============================================================================
# Place THIS file in your code project, e.g.
#   <path_main>/MonteCarlo_simulation/matrix.loading.main.R
# Place loading_structure_general.R in:
#   <path_main>/functions/functions_sim/loading_structure_general.R
# Keep your ORIGINAL matrix.dgp.R and matrix.sim.R unchanged.
# =============================================================================

# 1. Project paths
path_main <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

path_func_sim <- file.path(path_main, "R", "simulation")
source(file.path(path_func_sim, "matrix.dgp.R"))
source(file.path(path_func_sim, "loading_structure_general.R"))

# The scientific-brackets figures use only base R grid (no ggplot2/scales).
# MASS remains a dependency of the original DGP.
if (!requireNamespace("MASS", quietly = TRUE)) stop("Install MASS.")

# ============================================================================
# 2. YOUR CHOICES
# ============================================================================

# CHOOSE ONLY ONE:
r_choice <- c(2L, 2L)  # Four predictive factors as in Figure 2
K <- prod(r_choice)

# Original Monte Carlo settings (only change if you want a different design).
Tm_choice   <- 201L     # Design used in the paper
p1_choice   <- 20L
p2_choice   <- 50L
rhoG_choice <- 0.95
rhoE_choice <- 0.90
seed_choice <- 1712L

# First factor always has the original rhoF = 0.90.
# The remaining K-1 factors need different persistences if the ONLY proxy
# is scalar GDP with lags (equal AR dynamics give insufficient proxy rank).
# NULL uses (2,1): (.90,.60), (1,2): (.90,.45),
#           (2,2): (.90,.60,.45,.30).
# These are products of diagonal A_F and B_F, preserving matrix AR dynamics.
rhoF_choice <- NULL

# GDP/proxy include ALL K predictive factors. NULL = equal unit weights.
# Example r=(2,2): c(1,1,1,1). Must have K nonzero values.
target_weights_choice <- NULL

# Number of lags of the SAME simulated quarterly GDP proxy.
# Must be at least K. More lags enlarge the proxy space but are not guaranteed
# to strengthen weak factors; check reported proxy_strength and angles.
Lproxy_choice <- max(4L, K)
projected_updates_choice <- 1L

# ============================================================================
# 3. SAME ORIGINAL DGP PARAMETERS
# ============================================================================
# IMPORTANT: baseline matrix.dgp.R ONLY permits r1=r2=1.
# The new function adds the requested r_choice to its simulated truth.
dgp_params <- make_dgp_params(
  Tm = Tm_choice,
  burn = 400,
  draws = 1,                    # one diagnostic draw, not 500 forecasts
  seed = seed_choice,
  p1 = p1_choice,
  p2 = p2_choice,
  r1 = 1,                       # keep baseline generator restriction
  r2 = 1,
  v1 = 1,
  v2 = 4,
  rhoF = 0.90,
  rhoG = rhoG_choice,
  rhoE = rhoE_choice,
  R2_BI = 0.50,
  sigma_z = 0.25,
  share_monthly = 0.80,
  share_m_flow = 0.50,
  share_q_flow = 0.50,
  delta_r = 0.30,
  delta_c = 0.30,
  quarterly_obs_month = 3,
  ragged_edge = TRUE
)

# ============================================================================
# 4. FIT AND SAVE
# ============================================================================
suffix <- paste0("loading_r", r_choice[1L], "x", r_choice[2L],
                 "_T", Tm_choice, "_p", p1_choice, "x", p2_choice,
                 "_rG", rhoG_choice, "_rE", rhoE_choice)
out_dir <- file.path(path_main,
                     "results", "monte_carlo", suffix)

loading_result <- run_loading_structure_comparison(
  dgp_params = dgp_params,
  r = r_choice,
  out_dir = out_dir,
  rhoF_factors = rhoF_choice,
  target_weights = target_weights_choice,
  Lproxy = Lproxy_choice,
  projected_updates = projected_updates_choice
)

message("Completed loading experiment: ", paste(r_choice, collapse = "x"))
message("Results saved to: ", out_dir)
message("Proxy strength d[K]/d[1] = ",
        signif(loading_result$fit$proxy_strength, 3))
