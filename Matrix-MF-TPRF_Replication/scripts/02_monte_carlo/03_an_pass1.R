# Author: Davide Delfino
# Input: data/, config/settings.R, and functions in R/.
# Output: results/; see README.md.
# Original: MonteCarlo_simulation/mc.asymptotic.normality.R

# ==============================================================================
# MATRIX MF-TPRF
# Monte Carlo: Asymptotic Normality of Pass 1 Row Loadings
# ==============================================================================
# Author      : Davide Delfino
# Institution : Aarhus University
#
# Purpose:
# Verify the asymptotic distribution of the one-step projected row-loading
# estimator derived for Matrix MF-TPRF.
# ==============================================================================


# ==============================================================================
# 0. WORKING DIRECTORY
# ==============================================================================

path_main <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Run from the project root (see README.md).

path_results_root <- file.path(
  path_main,
  "results", "monte_carlo", "an_pass1"
)

if (!dir.exists(path_results_root)) {
  dir.create(path_results_root, recursive = TRUE)
}


# ==============================================================================
# 1. LIBRARIES
# ==============================================================================

library(ggplot2)


# ==============================================================================
# 2. MONTE CARLO PARAMETERS
# ==============================================================================

seed_population <- 1712
seed_sample     <- 1812
seed_mc         <- 1912

B <- 2000

T  <- 200
p1 <- 100
p2 <- 100

r1 <- 1
r2 <- 1
L  <- 1

g1 <- 1
g2 <- 4

i_fixed <- 1
a_fixed <- 1
j_fixed <- 1
b_fixed <- 1

sigma_z <- 0.25

delta_r <- 0.30
delta_c <- 0.30

d <- seq(
  from = 1.5,
  to = 0.5,
  length.out = r1
)

# ------------------------------------------------------------------------------
# Folder names
# ------------------------------------------------------------------------------

scenario_name <- paste0(
  "T", T,
  "_p1", p1,
  "_p2", p2
)

path_scenario <- file.path(
  path_results_root,
  scenario_name
)

path_results_R <- file.path(
  path_scenario,
  "row_loadings"
)

path_results_C <- file.path(
  path_scenario,
  "column_loadings"
)

dir.create(
  path_results_R,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  path_results_C,
  recursive = TRUE,
  showWarnings = FALSE
)
# ==============================================================================
# 3. BASIC UTILITIES AND DIMENSION CHECKS
# ==============================================================================


# ------------------------------------------------------------------------------
# Safe diagonal matrix
# ------------------------------------------------------------------------------

diag_matrix <- function(x) {
  
  x <- as.numeric(x)
  n <- length(x)
  
  diag(
    x,
    nrow = n,
    ncol = n
  )
}


# ------------------------------------------------------------------------------
# Extract one matrix from a three-dimensional array
# ------------------------------------------------------------------------------

matrix_slice <- function(A, t) {
  
  stopifnot(length(dim(A)) == 3)
  
  matrix(
    A[t, , , drop = FALSE],
    nrow = dim(A)[2],
    ncol = dim(A)[3]
  )
}


# ------------------------------------------------------------------------------
# Check admissible dimensions
# ------------------------------------------------------------------------------

check_design_dimensions <- function(
    p1,
    p2,
    r1,
    r2,
    g1,
    g2,
    L,
    d
) {
  
  if (r1 < 1 || r1 > p1) {
    stop("r1 must satisfy 1 <= r1 <= p1.")
  }
  
  if (r2 < 1 || r2 > p2) {
    stop("r2 must satisfy 1 <= r2 <= p2.")
  }
  
  if (g1 < 1 || r1 + g1 > p1) {
    stop("g1 must satisfy r1 + g1 <= p1.")
  }
  
  if (g2 < 1 || r2 + g2 > p2) {
    stop("g2 must satisfy r2 + g2 <= p2.")
  }
  
  if (L < 1) {
    stop("L must be at least one.")
  }
  
  if (length(d) != r1) {
    stop("d must have length r1.")
  }
  
  if (any(d <= 0)) {
    stop("All entries of d must be positive.")
  }
  
  # Necessary conditions for the proxy moments to be able
  # to span both predictive loading dimensions
  if (L * r2 < r1) {
    stop("Need L * r2 >= r1 for row-space identification.")
  }
  
  if (L * r1 < r2) {
    stop("Need L * r1 >= r2 for column-space identification.")
  }
  
  invisible(TRUE)
}


check_design_dimensions(
  p1 = p1,
  p2 = p2,
  r1 = r1,
  r2 = r2,
  g1 = g1,
  g2 = g2,
  L  = L,
  d  = d
)


# ==============================================================================
# 4. FIXED POPULATION DESIGN
# ==============================================================================

set.seed(seed_population)


# ------------------------------------------------------------------------------
# Predictive loading matrices
# ------------------------------------------------------------------------------

make_loading_matrix <- function(n, r) {
  
  if (r > n) {
    stop("The loading rank cannot exceed the number of rows.")
  }
  
  B <- matrix(
    runif(n * r, min = -1, max = 1),
    nrow = n,
    ncol = r
  )
  
  qr_B <- qr(B)
  
  if (qr_B$rank < r) {
    stop("Random loading matrix is rank deficient.")
  }
  
  Q <- qr.Q(qr_B)[
    ,
    seq_len(r),
    drop = FALSE
  ]
  
  B <- sqrt(n) * Q
  
  return(B)
}


# ------------------------------------------------------------------------------
# Orthogonalize non-predictive loadings
# ------------------------------------------------------------------------------

orthogonalize_loading_matrix <- function(
    B_pred,
    B_nonpred
) {
  
  n <- nrow(B_pred)
  g <- ncol(B_nonpred)
  
  if (ncol(B_pred) + g > n) {
    stop("Not enough dimensions for orthogonal non-predictive loadings.")
  }
  
  B_proj <- B_pred %*%
    solve(
      crossprod(B_pred),
      crossprod(B_pred, B_nonpred)
    )
  
  B_perp <- B_nonpred - B_proj
  
  qr_perp <- qr(B_perp)
  
  if (qr_perp$rank < g) {
    stop("Orthogonalized non-predictive loading matrix is rank deficient.")
  }
  
  Q <- qr.Q(qr_perp)[
    ,
    seq_len(g),
    drop = FALSE
  ]
  
  B_perp <- sqrt(n) * Q
  
  return(B_perp)
}


R0 <- make_loading_matrix(
  n = p1,
  r = r1
)

C0 <- make_loading_matrix(
  n = p2,
  r = r2
)

R0_g_raw <- make_loading_matrix(
  n = p1,
  r = g1
)

C0_g_raw <- make_loading_matrix(
  n = p2,
  r = g2
)

R0_g <- orthogonalize_loading_matrix(
  B_pred = R0,
  B_nonpred = R0_g_raw
)

C0_g <- orthogonalize_loading_matrix(
  B_pred = C0,
  B_nonpred = C0_g_raw
)


# ------------------------------------------------------------------------------
# Check loading normalizations and orthogonality
# ------------------------------------------------------------------------------

tol_loading <- 1e-8

stopifnot(
  max(abs(crossprod(R0) / p1 - diag(r1))) < tol_loading,
  max(abs(crossprod(C0) / p2 - diag(r2))) < tol_loading,
  max(abs(crossprod(R0_g) / p1 - diag(g1))) < tol_loading,
  max(abs(crossprod(C0_g) / p2 - diag(g2))) < tol_loading,
  max(abs(crossprod(R0, R0_g) / p1)) < tol_loading,
  max(abs(crossprod(C0, C0_g) / p2)) < tol_loading
)


# ------------------------------------------------------------------------------
# Idiosyncratic covariance matrices
# ------------------------------------------------------------------------------

make_covariance_matrix <- function(n, delta) {
  
  Sigma <- diag(1, n)
  
  if (n > 1) {
    
    for (i in seq_len(n - 1)) {
      
      Sigma[i, i + 1] <- delta
      Sigma[i + 1, i] <- delta
    }
  }
  
  eigenvalues <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (min(eigenvalues) <= 0) {
    stop("The covariance matrix is not positive definite.")
  }
  
  return(Sigma)
}


Sigma_row <- make_covariance_matrix(
  n = p1,
  delta = delta_r
)

Sigma_col <- make_covariance_matrix(
  n = p2,
  delta = delta_c
)


# Cholesky factors are fixed across Monte Carlo replications

L_Xi_row <- t(chol(Sigma_row))
L_Xi_col <- t(chol(Sigma_col))


# ==============================================================================
# 5. POPULATION FACTOR-PROXY DESIGN
# ==============================================================================


# ------------------------------------------------------------------------------
# General proxy design
#
# z_l,tau = <W_l, F_tau>_F + nu_l,tau
#
# Under
#
# F_tau = D^(1/2) U_tau,
#
# we have
#
# A_l,0 = E[z_l,tau F_tau] = D W_l.
# ------------------------------------------------------------------------------

make_proxy_factor_design <- function(
    r1,
    r2,
    L,
    d,
    max_tries = 1000,
    tol_identification = 1e-8,
    tol_eigengap = 1e-8
) {
  
  stopifnot(length(d) == r1)
  stopifnot(all(d > 0))
  
  if (L * r2 < r1) {
    stop("Need L * r2 >= r1 for row-space identification.")
  }
  
  if (L * r1 < r2) {
    stop("Need L * r1 >= r2 for column-space identification.")
  }
  
  D <- diag_matrix(d)
  
  identified <- FALSE
  
  for (attempt in seq_len(max_tries)) {
    
    W_list  <- vector("list", L)
    A0_list <- vector("list", L)
    
    
    # --------------------------------------------------------------------------
    # Construct proxy weights
    # --------------------------------------------------------------------------
    
    for (l in seq_len(L)) {
      
      W_l <- matrix(
        rnorm(r1 * r2),
        nrow = r1,
        ncol = r2
      )
      
      # Variance of the proxy signal <W_l, F_tau>
      signal_var <- sum(
        d * rowSums(W_l^2)
      )
      
      if (signal_var <= 0) {
        stop("Non-positive proxy signal variance.")
      }
      
      # Normalize each proxy signal to variance one
      W_l <- W_l / sqrt(signal_var)
      
      # Fix an arbitrary sign convention
      idx_sign <- which.max(abs(W_l))
      
      if (W_l[idx_sign] < 0) {
        W_l <- -W_l
      }
      
      W_list[[l]] <- W_l
      
      # Population moment A_l,0 = E[z_l F] = D W_l
      A0_list[[l]] <- D %*% W_l
    }
    
    
    # --------------------------------------------------------------------------
    # Population row and column identification matrices
    # --------------------------------------------------------------------------
    
    B_R <- matrix(
      0,
      nrow = r1,
      ncol = r1
    )
    
    B_C <- matrix(
      0,
      nrow = r2,
      ncol = r2
    )
    
    for (l in seq_len(L)) {
      
      A_l <- A0_list[[l]]
      
      B_R <- B_R +
        A_l %*%
        t(A_l)
      
      B_C <- B_C +
        t(A_l) %*%
        A_l
    }
    
    
    eig_R <- eigen(
      B_R,
      symmetric = TRUE
    )
    
    eig_C <- eigen(
      B_C,
      symmetric = TRUE
    )
    
    
    positive_R <-
      min(eig_R$values) > tol_identification
    
    positive_C <-
      min(eig_C$values) > tol_identification
    
    
    distinct_R <- if (r1 == 1) {
      TRUE
    } else {
      min(abs(diff(eig_R$values))) > tol_eigengap
    }
    
    distinct_C <- if (r2 == 1) {
      TRUE
    } else {
      min(abs(diff(eig_C$values))) > tol_eigengap
    }
    
    
    if (
      positive_R &&
      positive_C &&
      distinct_R &&
      distinct_C
    ) {
      
      identified <- TRUE
      break
    }
  }
  
  
  if (!identified) {
    stop("Could not generate an identified proxy design.")
  }
  
  
  # ---------------------------------------------------------------------------
  # Population covariance of proxy signal components
  # ---------------------------------------------------------------------------
  
  Sigma_signal <- matrix(
    0,
    nrow = L,
    ncol = L
  )
  
  for (l in seq_len(L)) {
    
    for (m in seq_len(L)) {
      
      Sigma_signal[l, m] <- sum(
        d *
          rowSums(
            W_list[[l]] *
              W_list[[m]]
          )
      )
    }
  }
  
  
  Gamma_R <- eig_R$vectors
  Lambda_R <- diag_matrix(eig_R$values)
  
  Gamma_C <- eig_C$vectors
  Lambda_C <- diag_matrix(eig_C$values)
  
  
  return(
    list(
      W_list = W_list,
      A0_list = A0_list,
      Sigma_signal = Sigma_signal,
      B_R = B_R,
      B_C = B_C,
      Gamma_R = Gamma_R,
      Lambda_R = Lambda_R,
      Gamma_C = Gamma_C,
      Lambda_C = Lambda_C
    )
  )
}


proxy_factor_population <- make_proxy_factor_design(
  r1 = r1,
  r2 = r2,
  L = L,
  d = d
)


# ==============================================================================
# 6. DATA-GENERATING PROCESS
# ==============================================================================


# ------------------------------------------------------------------------------
# Predictive factors
# ------------------------------------------------------------------------------

simulate_predictive_factors <- function(
    T,
    r1,
    r2,
    d
) {
  
  stopifnot(length(d) == r1)
  stopifnot(all(d > 0))
  
  D_sqrt <- diag_matrix(
    sqrt(d)
  )
  
  F0 <- array(
    0,
    dim = c(T, r1, r2)
  )
  
  for (tau in seq_len(T)) {
    
    U_tau <- matrix(
      rnorm(r1 * r2),
      nrow = r1,
      ncol = r2
    )
    
    F0[tau, , ] <-
      D_sqrt %*% U_tau
  }
  
  return(F0)
}


# ------------------------------------------------------------------------------
# Proxy variables
# ------------------------------------------------------------------------------

simulate_proxies <- function(
    F0,
    W_list,
    sigma_z
) {
  
  T  <- dim(F0)[1]
  r1 <- dim(F0)[2]
  r2 <- dim(F0)[3]
  
  L <- length(W_list)
  
  z <- matrix(
    0,
    nrow = T,
    ncol = L
  )
  
  for (l in seq_len(L)) {
    
    if (!all(dim(W_list[[l]]) == c(r1, r2))) {
      stop("A proxy weight matrix W_l has incompatible dimensions.")
    }
  }
  
  
  for (tau in seq_len(T)) {
    
    F_tau <- matrix_slice(
      F0,
      tau
    )
    
    nu_tau <- rnorm(
      L,
      mean = 0,
      sd = sigma_z
    )
    
    for (l in seq_len(L)) {
      
      proxy_signal <- sum(
        W_list[[l]] *
          F_tau
      )
      
      z[tau, l] <-
        proxy_signal +
        nu_tau[l]
    }
  }
  
  return(z)
}


# ------------------------------------------------------------------------------
# Non-predictive factors
# ------------------------------------------------------------------------------

simulate_nonpredictive_factors <- function(
    T,
    g1,
    g2
) {
  
  G0 <- array(
    0,
    dim = c(T, g1, g2)
  )
  
  for (tau in seq_len(T)) {
    
    G0[tau, , ] <- matrix(
      rnorm(g1 * g2),
      nrow = g1,
      ncol = g2
    )
  }
  
  return(G0)
}


# ------------------------------------------------------------------------------
# Matrix-normal idiosyncratic error
# ------------------------------------------------------------------------------

draw_matrix_normal <- function(
    L_row,
    L_col
) {
  
  p1 <- nrow(L_row)
  p2 <- nrow(L_col)
  
  Z <- matrix(
    rnorm(p1 * p2),
    nrow = p1,
    ncol = p2
  )
  
  Xi <-
    L_row %*%
    Z %*%
    t(L_col)
  
  return(Xi)
}


simulate_idiosyncratic_errors <- function(
    T,
    L_row,
    L_col
) {
  
  p1 <- nrow(L_row)
  p2 <- nrow(L_col)
  
  Xi0 <- array(
    0,
    dim = c(T, p1, p2)
  )
  
  for (tau in seq_len(T)) {
    
    Xi0[tau, , ] <- draw_matrix_normal(
      L_row = L_row,
      L_col = L_col
    )
  }
  
  return(Xi0)
}


# ------------------------------------------------------------------------------
# Observed predictor matrices
# ------------------------------------------------------------------------------

construct_predictor_matrices <- function(
    F0,
    G0,
    Xi0,
    R0,
    C0,
    R0_g,
    C0_g
) {
  
  T  <- dim(F0)[1]
  p1 <- nrow(R0)
  p2 <- nrow(C0)
  
  X <- array(
    0,
    dim = c(T, p1, p2)
  )
  
  for (tau in seq_len(T)) {
    
    F_tau <- matrix_slice(
      F0,
      tau
    )
    
    G_tau <- matrix_slice(
      G0,
      tau
    )
    
    Xi_tau <- matrix_slice(
      Xi0,
      tau
    )
    
    
    X_pred <-
      R0 %*%
      F_tau %*%
      t(C0)
    
    
    X_nonpred <-
      R0_g %*%
      G_tau %*%
      t(C0_g)
    
    
    X[tau, , ] <-
      X_pred +
      X_nonpred +
      Xi_tau
  }
  
  return(X)
}


# ==============================================================================
# 7. THEORETICAL ASYMPTOTIC COVARIANCE
# ==============================================================================

compute_theoretical_covariance <- function(
    C0,
    Sigma_row,
    Sigma_col,
    A0_list,
    Sigma_signal,
    Gamma_R,
    Lambda_R,
    sigma_z,
    i
) {
  
  p2 <- nrow(C0)
  r1 <- nrow(Lambda_R)
  L  <- length(A0_list)
  
  if (i < 1 || i > nrow(Sigma_row)) {
    stop("Invalid fixed row i.")
  }
  
  
  # ---------------------------------------------------------------------------
  # Column-side covariance after projection by C0
  # ---------------------------------------------------------------------------
  
  Q_C <-
    crossprod(
      C0,
      Sigma_col %*% C0
    ) / p2
  
  
  # ---------------------------------------------------------------------------
  # Full population covariance matrix of the proxies
  # ---------------------------------------------------------------------------
  
  Sigma_z <-
    Sigma_signal +
    sigma_z^2 * diag(L)
  
  
  # ---------------------------------------------------------------------------
  # Variance of the limiting score for row i
  # ---------------------------------------------------------------------------
  
  sigma_row_i <-
    Sigma_row[i, i]
  
  
  V_R_i <- matrix(
    0,
    nrow = r1,
    ncol = r1
  )
  
  
  for (l in seq_len(L)) {
    
    for (m in seq_len(L)) {
      
      V_R_i <- V_R_i +
        Sigma_z[l, m] *
        A0_list[[l]] %*%
        Q_C %*%
        t(A0_list[[m]])
    }
  }
  
  
  V_R_i <-
    sigma_row_i * V_R_i
  
  V_R_i <-
    (V_R_i + t(V_R_i)) / 2
  
  
  # ---------------------------------------------------------------------------
  # Asymptotic covariance of the row-loading estimator
  # ---------------------------------------------------------------------------
  
  Lambda_R_inv <-
    solve(Lambda_R)
  
  
  Omega_R_i <-
    Lambda_R_inv %*%
    t(Gamma_R) %*%
    V_R_i %*%
    Gamma_R %*%
    Lambda_R_inv
  
  
  Omega_R_i <-
    (Omega_R_i + t(Omega_R_i)) / 2
  
  
  return(
    list(
      Q_C = Q_C,
      Sigma_z = Sigma_z,
      V_R_i = V_R_i,
      Omega_R_i = Omega_R_i
    )
  )
}


theoretical_covariance_R <- compute_theoretical_covariance(
  C0 = C0,
  Sigma_row = Sigma_row,
  Sigma_col = Sigma_col,
  A0_list = proxy_factor_population$A0_list,
  Sigma_signal = proxy_factor_population$Sigma_signal,
  Gamma_R = proxy_factor_population$Gamma_R,
  Lambda_R = proxy_factor_population$Lambda_R,
  sigma_z = sigma_z,
  i = i_fixed
)



compute_theoretical_covariance_C <- function(
    R0,
    Sigma_row,
    Sigma_col,
    A0_list,
    Sigma_signal,
    Gamma_C,
    Lambda_C,
    sigma_z,
    j
) {
  
  p1 <- nrow(R0)
  r2 <- nrow(Lambda_C)
  L  <- length(A0_list)
  
  
  if (j < 1 || j > nrow(Sigma_col)) {
    stop("Invalid fixed column j.")
  }
  
  
  # ---------------------------------------------------------------------------
  # Row-side covariance after projection by R0
  # ---------------------------------------------------------------------------
  
  Q_R <-
    crossprod(
      R0,
      Sigma_row %*% R0
    ) / p1
  
  
  # ---------------------------------------------------------------------------
  # Full proxy covariance
  # ---------------------------------------------------------------------------
  
  Sigma_z <-
    Sigma_signal +
    sigma_z^2 * diag(L)
  
  
  # ---------------------------------------------------------------------------
  # Limiting score variance for column j
  # ---------------------------------------------------------------------------
  
  sigma_col_j <-
    Sigma_col[j, j]
  
  
  V_C_j <- matrix(
    0,
    nrow = r2,
    ncol = r2
  )
  
  
  for (l in seq_len(L)) {
    
    for (m in seq_len(L)) {
      
      V_C_j <- V_C_j +
        Sigma_z[l, m] *
        t(A0_list[[l]]) %*%
        Q_R %*%
        A0_list[[m]]
    }
  }
  
  
  V_C_j <-
    sigma_col_j * V_C_j
  
  V_C_j <-
    (V_C_j + t(V_C_j)) / 2
  
  
  # ---------------------------------------------------------------------------
  # Asymptotic covariance
  # ---------------------------------------------------------------------------
  
  Lambda_C_inv <-
    solve(Lambda_C)
  
  
  Omega_C_j <-
    Lambda_C_inv %*%
    t(Gamma_C) %*%
    V_C_j %*%
    Gamma_C %*%
    Lambda_C_inv
  
  
  Omega_C_j <-
    (Omega_C_j + t(Omega_C_j)) / 2
  
  
  return(
    list(
      Q_R = Q_R,
      Sigma_z = Sigma_z,
      V_C_j = V_C_j,
      Omega_C_j = Omega_C_j
    )
  )
}


theoretical_covariance_C <- compute_theoretical_covariance_C(
  R0 = R0,
  Sigma_row = Sigma_row,
  Sigma_col = Sigma_col,
  A0_list = proxy_factor_population$A0_list,
  Sigma_signal = proxy_factor_population$Sigma_signal,
  Gamma_C = proxy_factor_population$Gamma_C,
  Lambda_C = proxy_factor_population$Lambda_C,
  sigma_z = sigma_z,
  j = j_fixed
)

# ==============================================================================
# 8. PASS 1: ONE-STEP PROJECTED ROW ESTIMATOR
# ==============================================================================

mf_tprf_pass1_one_step <- function(
    X,
    z,
    r1,
    r2
) {
  
  stopifnot(length(dim(X)) == 3)
  
  T  <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(X)[3]
  
  z <- as.matrix(z)
  
  if (nrow(z) != T) {
    stop("nrow(z) must equal dim(X)[1].")
  }
  
  if (r1 < 1 || r1 > p1) {
    stop("Invalid r1.")
  }
  
  if (r2 < 1 || r2 > p2) {
    stop("Invalid r2.")
  }
  
  L <- ncol(z)
  
  
  # ---------------------------------------------------------------------------
  # 1. Proxy-weighted moments
  #
  # M_l = (1/T) sum_tau z_tau,l X_tau
  # ---------------------------------------------------------------------------
  
  M_list <- vector(
    "list",
    L
  )
  
  
  for (l in seq_len(L)) {
    
    M_l <- matrix(
      0,
      nrow = p1,
      ncol = p2
    )
    
    
    for (tau in seq_len(T)) {
      
      X_tau <- matrix_slice(
        X,
        tau
      )
      
      M_l <- M_l +
        z[tau, l] *
        X_tau
    }
    
    
    M_list[[l]] <-
      M_l / T
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Spectral initialization
  # ---------------------------------------------------------------------------
  
  S_R_0 <- matrix(
    0,
    nrow = p1,
    ncol = p1
  )
  
  S_C_0 <- matrix(
    0,
    nrow = p2,
    ncol = p2
  )
  
  
  for (l in seq_len(L)) {
    
    M_l <- M_list[[l]]
    
    S_R_0 <- S_R_0 +
      M_l %*%
      t(M_l)
    
    S_C_0 <- S_C_0 +
      t(M_l) %*%
      M_l
  }
  
  
  eig_R_0 <- eigen(
    S_R_0,
    symmetric = TRUE
  )
  
  eig_C_0 <- eigen(
    S_C_0,
    symmetric = TRUE
  )
  
  
  R_hat_0 <-
    sqrt(p1) *
    eig_R_0$vectors[
      ,
      seq_len(r1),
      drop = FALSE
    ]
  
  
  C_hat_0 <-
    sqrt(p2) *
    eig_C_0$vectors[
      ,
      seq_len(r2),
      drop = FALSE
    ]
  
  
  # ---------------------------------------------------------------------------
  # 3. One projected row update
  # ---------------------------------------------------------------------------
  
  P_C_0 <-
    C_hat_0 %*%
    t(C_hat_0) / p2
  
  
  S_R_1 <- matrix(
    0,
    nrow = p1,
    ncol = p1
  )
  
  
  for (l in seq_len(L)) {
    
    M_l <- M_list[[l]]
    
    S_R_1 <- S_R_1 +
      M_l %*%
      P_C_0 %*%
      t(M_l)
  }
  
  
  eig_R_1 <- eigen(
    S_R_1,
    symmetric = TRUE
  )
  
  
  leading_values <-
    eig_R_1$values[
      seq_len(r1)
    ]
  
  
  if (min(leading_values) <= 0) {
    stop("The leading projected row eigenvalues must be positive.")
  }
  
  
  R_hat_1 <-
    sqrt(p1) *
    eig_R_1$vectors[
      ,
      seq_len(r1),
      drop = FALSE
    ]
  
  
  D_hat_R_1 <-
    diag_matrix(
      leading_values
    )
  
  # ---------------------------------------------------------------------------
  # 4. One projected column update
  # ---------------------------------------------------------------------------
  
  P_R_1 <-
    R_hat_1 %*%
    t(R_hat_1) / p1
  
  
  S_C_1 <- matrix(
    0,
    nrow = p2,
    ncol = p2
  )
  
  
  for (l in seq_len(L)) {
    
    M_l <- M_list[[l]]
    
    S_C_1 <- S_C_1 +
      t(M_l) %*%
      P_R_1 %*%
      M_l
  }
  
  
  eig_C_1 <- eigen(
    S_C_1,
    symmetric = TRUE
  )
  
  
  leading_values_C <-
    eig_C_1$values[
      seq_len(r2)
    ]
  
  
  if (min(leading_values_C) <= 0) {
    stop("The leading projected column eigenvalues must be positive.")
  }
  
  
  C_hat_1 <-
    sqrt(p2) *
    eig_C_1$vectors[
      ,
      seq_len(r2),
      drop = FALSE
    ]
  
  
  D_hat_C_1 <-
    diag_matrix(
      leading_values_C
    )
  
  return(
    list(
      M_list = M_list,
      S_R_0 = S_R_0,
      S_C_0 = S_C_0,
      R_hat_0 = R_hat_0,
      C_hat_0 = C_hat_0,
      P_C_0 = P_C_0,
      S_R_1 = S_R_1,
      R_hat_1 = R_hat_1,
      D_hat_R_1 = D_hat_R_1,
      P_R_1 = P_R_1,
      S_C_1 = S_C_1,
      C_hat_1 = C_hat_1,
      D_hat_C_1 = D_hat_C_1
    )
  )
}


# ==============================================================================
# 9. IDENTIFICATION / ROTATION MATRIX
# ==============================================================================


# ------------------------------------------------------------------------------
# Theorem-defined H_R^(1)
# ------------------------------------------------------------------------------

compute_H_R_theory <- function(
    F0,
    z,
    C0,
    C_hat_0,
    R0,
    R_hat_1,
    D_hat_R_1
) {
  
  T  <- dim(F0)[1]
  r1 <- dim(F0)[2]
  r2 <- dim(F0)[3]
  
  p1 <- nrow(R0)
  p2 <- nrow(C0)
  
  z <- as.matrix(z)
  L <- ncol(z)
  
  
  if (nrow(z) != T) {
    stop("nrow(z) must equal the number of factor observations.")
  }
  
  if (!all(dim(C_hat_0) == c(p2, r2))) {
    stop("C_hat_0 has incompatible dimensions.")
  }
  
  if (!all(dim(R_hat_1) == c(p1, r1))) {
    stop("R_hat_1 has incompatible dimensions.")
  }
  
  if (!all(dim(D_hat_R_1) == c(r1, r1))) {
    stop("D_hat_R_1 has incompatible dimensions.")
  }
  
  
  # ---------------------------------------------------------------------------
  # 1. Sample proxy-factor moments
  #
  # A_hat_l = (1/T) sum_tau z_tau,l F_tau
  # ---------------------------------------------------------------------------
  
  A_hat_list <- vector(
    "list",
    L
  )
  
  
  for (l in seq_len(L)) {
    
    A_hat_l <- matrix(
      0,
      nrow = r1,
      ncol = r2
    )
    
    
    for (tau in seq_len(T)) {
      
      F_tau <- matrix_slice(
        F0,
        tau
      )
      
      A_hat_l <- A_hat_l +
        z[tau, l] *
        F_tau
    }
    
    
    A_hat_list[[l]] <-
      A_hat_l / T
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Initial column-space overlap
  #
  # H_C,0 = C0' C_hat_0 / p2
  # ---------------------------------------------------------------------------
  
  H_C_0 <-
    crossprod(
      C0,
      C_hat_0
    ) / p2
  
  
  K_C <-
    H_C_0 %*%
    t(H_C_0)
  
  
  # ---------------------------------------------------------------------------
  # 3. Predictive sample matrix entering H_R^(1)
  # ---------------------------------------------------------------------------
  
  B_hat_R <- matrix(
    0,
    nrow = r1,
    ncol = r1
  )
  
  
  for (l in seq_len(L)) {
    
    B_hat_R <- B_hat_R +
      A_hat_list[[l]] %*%
      K_C %*%
      t(A_hat_list[[l]])
  }
  
  
  # ---------------------------------------------------------------------------
  # 4. True-estimated row overlap
  # ---------------------------------------------------------------------------
  
  H_overlap <-
    crossprod(
      R0,
      R_hat_1
    ) / p1
  
  
  # ---------------------------------------------------------------------------
  # 5. Scaled eigenvalue matrix
  # ---------------------------------------------------------------------------
  
  D_scaled <-
    D_hat_R_1 /
    (p1 * p2)
  
  
  if (min(diag(D_scaled)) <= 0) {
    stop("D_scaled must have positive diagonal elements.")
  }
  
  
  # ---------------------------------------------------------------------------
  # 6. Theorem-defined row rotation
  # ---------------------------------------------------------------------------
  
  H_R_1 <-
    B_hat_R %*%
    H_overlap %*%
    solve(D_scaled)
  
  
  return(
    list(
      H_R_1 = H_R_1,
      A_hat_list = A_hat_list,
      H_C_0 = H_C_0,
      K_C = K_C,
      B_hat_R = B_hat_R,
      H_overlap = H_overlap,
      D_scaled = D_scaled
    )
  )
}

compute_H_C_theory <- function(
    F0,
    z,
    R0,
    R_hat_1,
    C0,
    C_hat_1,
    D_hat_C_1
) {
  
  T  <- dim(F0)[1]
  r1 <- dim(F0)[2]
  r2 <- dim(F0)[3]
  
  p1 <- nrow(R0)
  p2 <- nrow(C0)
  
  z <- as.matrix(z)
  L <- ncol(z)
  
  
  # ---------------------------------------------------------------------------
  # 1. Sample proxy-factor moments
  # ---------------------------------------------------------------------------
  
  A_hat_list <- vector(
    "list",
    L
  )
  
  
  for (l in seq_len(L)) {
    
    A_hat_l <- matrix(
      0,
      nrow = r1,
      ncol = r2
    )
    
    for (tau in seq_len(T)) {
      
      F_tau <- matrix_slice(
        F0,
        tau
      )
      
      A_hat_l <- A_hat_l +
        z[tau, l] *
        F_tau
    }
    
    A_hat_list[[l]] <-
      A_hat_l / T
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Refined row-space overlap
  # ---------------------------------------------------------------------------
  
  H_R_1 <-
    crossprod(
      R0,
      R_hat_1
    ) / p1
  
  
  K_R <-
    H_R_1 %*%
    t(H_R_1)
  
  
  # ---------------------------------------------------------------------------
  # 3. Predictive sample matrix
  # ---------------------------------------------------------------------------
  
  B_hat_C <- matrix(
    0,
    nrow = r2,
    ncol = r2
  )
  
  
  for (l in seq_len(L)) {
    
    B_hat_C <- B_hat_C +
      t(A_hat_list[[l]]) %*%
      K_R %*%
      A_hat_list[[l]]
  }
  
  
  # ---------------------------------------------------------------------------
  # 4. True-estimated column overlap
  # ---------------------------------------------------------------------------
  
  H_overlap <-
    crossprod(
      C0,
      C_hat_1
    ) / p2
  
  
  # ---------------------------------------------------------------------------
  # 5. Scaled eigenvalues
  # ---------------------------------------------------------------------------
  
  D_scaled <-
    D_hat_C_1 /
    (p1 * p2)
  
  
  # ---------------------------------------------------------------------------
  # 6. Theorem-defined column rotation
  # ---------------------------------------------------------------------------
  
  H_C_1 <-
    B_hat_C %*%
    H_overlap %*%
    solve(D_scaled)
  
  
  return(
    list(
      H_C_1 = H_C_1,
      A_hat_list = A_hat_list,
      H_R_1 = H_R_1,
      K_R = K_R,
      B_hat_C = B_hat_C,
      H_overlap = H_overlap,
      D_scaled = D_scaled
    )
  )
}

compute_H_C_LS <- function(
    C0,
    C_hat_1
) {
  
  p2 <- nrow(C0)
  
  H_C_LS <-
    crossprod(
      C0,
      C_hat_1
    ) / p2
  
  return(H_C_LS)
}

# ------------------------------------------------------------------------------
# Least-squares rotation used only as a diagnostic
# ------------------------------------------------------------------------------

compute_H_R_LS <- function(
    R0,
    R_hat_1
) {
  
  p1 <- nrow(R0)
  
  H_R_LS <-
    crossprod(
      R0,
      R_hat_1
    ) / p1
  
  return(H_R_LS)
}


# ==============================================================================
# 10. CLT STANDARDIZATION
# ==============================================================================


# ------------------------------------------------------------------------------
# Inverse square root of a positive definite matrix
# ------------------------------------------------------------------------------

matrix_inverse_sqrt <- function(Sigma) {
  
  Sigma <- as.matrix(Sigma)
  
  if (nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be square.")
  }
  
  Sigma <-
    (Sigma + t(Sigma)) / 2
  
  
  eig <- eigen(
    Sigma,
    symmetric = TRUE
  )
  
  
  if (min(eig$values) <= 0) {
    stop("Sigma must be positive definite.")
  }
  
  
  D_inv_sqrt <-
    diag_matrix(
      1 / sqrt(eig$values)
    )
  
  
  Sigma_inv_sqrt <-
    eig$vectors %*%
    D_inv_sqrt %*%
    t(eig$vectors)
  
  
  return(Sigma_inv_sqrt)
}


# ------------------------------------------------------------------------------
# Fixed-row CLT statistic
# ------------------------------------------------------------------------------

compute_row_clt_statistic <- function(
    R_hat_1,
    R0,
    H_R_1,
    Omega_R_i,
    T,
    p2,
    i,
    a
) {
  
  r1 <- ncol(R0)
  
  if (i < 1 || i > nrow(R0)) {
    stop("Invalid row index i.")
  }
  
  if (a < 1 || a > r1) {
    stop("Invalid coordinate index a.")
  }
  
  
  # True loading row expressed in estimated coordinates
  
  R0_rotated <-
    R0[
      i,
      ,
      drop = FALSE
    ] %*%
    H_R_1
  
  
  # Row-loading estimation error
  
  row_error <-
    R_hat_1[
      i,
      ,
      drop = FALSE
    ] -
    R0_rotated
  
  
  # CLT scaling
  
  scaled_error <-
    sqrt(T * p2) *
    t(row_error)
  
  
  # Theoretical covariance standardization
  
  Omega_inv_sqrt <-
    matrix_inverse_sqrt(
      Omega_R_i
    )
  
  
  Z <-
    Omega_inv_sqrt %*%
    scaled_error
  
  
  return(
    list(
      row_error = as.numeric(row_error),
      scaled_error = as.numeric(scaled_error),
      Z = as.numeric(Z),
      Z_fixed = as.numeric(Z[a])
    )
  )
}

compute_column_clt_statistic <- function(
    C_hat_1,
    C0,
    H_C_1,
    Omega_C_j,
    T,
    p1,
    j,
    b
) {
  
  r2 <- ncol(C0)
  
  
  if (j < 1 || j > nrow(C0)) {
    stop("Invalid column-loading index j.")
  }
  
  if (b < 1 || b > r2) {
    stop("Invalid coordinate index b.")
  }
  
  
  C0_rotated <-
    C0[
      j,
      ,
      drop = FALSE
    ] %*%
    H_C_1
  
  
  column_error <-
    C_hat_1[
      j,
      ,
      drop = FALSE
    ] -
    C0_rotated
  
  
  scaled_error <-
    sqrt(T * p1) *
    t(column_error)
  
  
  Omega_inv_sqrt <-
    matrix_inverse_sqrt(
      Omega_C_j
    )
  
  
  Z <-
    Omega_inv_sqrt %*%
    scaled_error
  
  
  return(
    list(
      column_error = as.numeric(column_error),
      scaled_error = as.numeric(scaled_error),
      Z = as.numeric(Z),
      Z_fixed = as.numeric(Z[b])
    )
  )
}

# ==============================================================================
# 11. GENERATE AND INSPECT ONE SAMPLE
# ==============================================================================

set.seed(seed_sample)


# ------------------------------------------------------------------------------
# Predictive factors
# ------------------------------------------------------------------------------

F0 <- simulate_predictive_factors(
  T = T,
  r1 = r1,
  r2 = r2,
  d = d
)


# ------------------------------------------------------------------------------
# Proxies
# ------------------------------------------------------------------------------

z <- simulate_proxies(
  F0 = F0,
  W_list = proxy_factor_population$W_list,
  sigma_z = sigma_z
)


# ------------------------------------------------------------------------------
# Non-predictive factors
# ------------------------------------------------------------------------------

G0 <- simulate_nonpredictive_factors(
  T = T,
  g1 = g1,
  g2 = g2
)


# ------------------------------------------------------------------------------
# Idiosyncratic errors
# ------------------------------------------------------------------------------

Xi0 <- simulate_idiosyncratic_errors(
  T = T,
  L_row = L_Xi_row,
  L_col = L_Xi_col
)


# ------------------------------------------------------------------------------
# Observed predictor tensor
# ------------------------------------------------------------------------------

X <- construct_predictor_matrices(
  F0 = F0,
  G0 = G0,
  Xi0 = Xi0,
  R0 = R0,
  C0 = C0,
  R0_g = R0_g,
  C0_g = C0_g
)


# ------------------------------------------------------------------------------
# Pass 1
# ------------------------------------------------------------------------------

pass1 <- mf_tprf_pass1_one_step(
  X = X,
  z = z,
  r1 = r1,
  r2 = r2
)


# ==============================================================================
# 11.1 ROW LOADINGS
# ==============================================================================


# ------------------------------------------------------------------------------
# Theorem-defined row rotation
# ------------------------------------------------------------------------------

H_R_out <- compute_H_R_theory(
  F0 = F0,
  z = z,
  C0 = C0,
  C_hat_0 = pass1$C_hat_0,
  R0 = R0,
  R_hat_1 = pass1$R_hat_1,
  D_hat_R_1 = pass1$D_hat_R_1
)

H_R_1 <- H_R_out$H_R_1


# ------------------------------------------------------------------------------
# Least-squares row rotation
# ------------------------------------------------------------------------------

H_R_LS <- compute_H_R_LS(
  R0 = R0,
  R_hat_1 = pass1$R_hat_1
)


# ------------------------------------------------------------------------------
# Row CLT statistic
# ------------------------------------------------------------------------------

clt_sample_R <- compute_row_clt_statistic(
  R_hat_1 = pass1$R_hat_1,
  R0 = R0,
  H_R_1 = H_R_1,
  Omega_R_i = theoretical_covariance_R$Omega_R_i,
  T = T,
  p2 = p2,
  i = i_fixed,
  a = a_fixed
)


# ==============================================================================
# 11.2 COLUMN LOADINGS
# ==============================================================================


# ------------------------------------------------------------------------------
# Theorem-defined column rotation
# ------------------------------------------------------------------------------

H_C_out <- compute_H_C_theory(
  F0 = F0,
  z = z,
  R0 = R0,
  R_hat_1 = pass1$R_hat_1,
  C0 = C0,
  C_hat_1 = pass1$C_hat_1,
  D_hat_C_1 = pass1$D_hat_C_1
)

H_C_1 <- H_C_out$H_C_1


# ------------------------------------------------------------------------------
# Least-squares column rotation
# ------------------------------------------------------------------------------

H_C_LS <- compute_H_C_LS(
  C0 = C0,
  C_hat_1 = pass1$C_hat_1
)


# ------------------------------------------------------------------------------
# Column CLT statistic
# ------------------------------------------------------------------------------

clt_sample_C <- compute_column_clt_statistic(
  C_hat_1 = pass1$C_hat_1,
  C0 = C0,
  H_C_1 = H_C_1,
  Omega_C_j = theoretical_covariance_C$Omega_C_j,
  T = T,
  p1 = p1,
  j = j_fixed,
  b = b_fixed
)


# ==============================================================================
# 11.3 ONE-SAMPLE DIAGNOSTICS
# ==============================================================================


# ------------------------------------------------------------------------------
# Row-loading diagnostics
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ONE-SAMPLE DIAGNOSTICS: ROW LOADINGS\n")
cat("============================================================\n")

cat("\nTheorem-defined H_R:\n")
print(H_R_1)

cat("\nLeast-squares H_R:\n")
print(H_R_LS)

cat("\nFrobenius distance between row rotations:\n")
print(
  norm(
    H_R_1 - H_R_LS,
    type = "F"
  )
)

cat("\nScaled row-loading error:\n")
print(
  clt_sample_R$scaled_error
)

cat("\nStandardized row-loading CLT statistic:\n")
print(
  clt_sample_R$Z
)

cat("\nFixed row-loading coordinate:\n")
print(
  clt_sample_R$Z_fixed
)


# ------------------------------------------------------------------------------
# Column-loading diagnostics
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("ONE-SAMPLE DIAGNOSTICS: COLUMN LOADINGS\n")
cat("============================================================\n")

cat("\nTheorem-defined H_C:\n")
print(H_C_1)

cat("\nLeast-squares H_C:\n")
print(H_C_LS)

cat("\nFrobenius distance between column rotations:\n")
print(
  norm(
    H_C_1 - H_C_LS,
    type = "F"
  )
)

cat("\nScaled column-loading error:\n")
print(
  clt_sample_C$scaled_error
)

cat("\nStandardized column-loading CLT statistic:\n")
print(
  clt_sample_C$Z
)

cat("\nFixed column-loading coordinate:\n")
print(
  clt_sample_C$Z_fixed
)
# ==============================================================================
# 12. ONE MONTE CARLO REPLICATION
# ==============================================================================

run_one_replication <- function(
    T,
    d,
    sigma_z,
    W_list,
    R0,
    C0,
    R0_g,
    C0_g,
    L_Xi_row,
    L_Xi_col,
    Omega_R_i,
    i_fixed,
    a_fixed,
    Omega_C_j,
    j_fixed,
    b_fixed
) {
  
  r1 <- ncol(R0)
  r2 <- ncol(C0)
  
  g1 <- ncol(R0_g)
  g2 <- ncol(C0_g)
  
  p2 <- nrow(C0)
  
  
  # ---------------------------------------------------------------------------
  # 1. Predictive factors
  # ---------------------------------------------------------------------------
  
  F0 <- simulate_predictive_factors(
    T = T,
    r1 = r1,
    r2 = r2,
    d = d
  )
  
  
  # ---------------------------------------------------------------------------
  # 2. Proxies
  # ---------------------------------------------------------------------------
  
  z <- simulate_proxies(
    F0 = F0,
    W_list = W_list,
    sigma_z = sigma_z
  )
  
  
  # ---------------------------------------------------------------------------
  # 3. Non-predictive factors
  # ---------------------------------------------------------------------------
  
  G0 <- simulate_nonpredictive_factors(
    T = T,
    g1 = g1,
    g2 = g2
  )
  
  
  # ---------------------------------------------------------------------------
  # 4. Idiosyncratic errors
  # ---------------------------------------------------------------------------
  
  Xi0 <- simulate_idiosyncratic_errors(
    T = T,
    L_row = L_Xi_row,
    L_col = L_Xi_col
  )
  
  
  # ---------------------------------------------------------------------------
  # 5. Predictor tensor
  # ---------------------------------------------------------------------------
  
  X <- construct_predictor_matrices(
    F0 = F0,
    G0 = G0,
    Xi0 = Xi0,
    R0 = R0,
    C0 = C0,
    R0_g = R0_g,
    C0_g = C0_g
  )
  
  
  # ---------------------------------------------------------------------------
  # 6. One-step projected estimator
  # ---------------------------------------------------------------------------
  
  pass1 <- mf_tprf_pass1_one_step(
    X = X,
    z = z,
    r1 = r1,
    r2 = r2
  )
  
  
  # ---------------------------------------------------------------------------
  # 7. Theorem-defined rotation
  # ---------------------------------------------------------------------------
  
  H_out <- compute_H_R_theory(
    F0 = F0,
    z = z,
    C0 = C0,
    C_hat_0 = pass1$C_hat_0,
    R0 = R0,
    R_hat_1 = pass1$R_hat_1,
    D_hat_R_1 = pass1$D_hat_R_1
  )
  
  
  H_R_1 <- H_out$H_R_1
  
  # ---------------------------------------------------------------------------
  # Column rotation
  # ---------------------------------------------------------------------------
  
  H_C_out <- compute_H_C_theory(
    F0 = F0,
    z = z,
    R0 = R0,
    R_hat_1 = pass1$R_hat_1,
    C0 = C0,
    C_hat_1 = pass1$C_hat_1,
    D_hat_C_1 = pass1$D_hat_C_1
  )
  
  H_C_1 <- H_C_out$H_C_1
  
  
  H_C_LS <- compute_H_C_LS(
    C0 = C0,
    C_hat_1 = pass1$C_hat_1
  )
  
  
  clt_C <- compute_column_clt_statistic(
    C_hat_1 = pass1$C_hat_1,
    C0 = C0,
    H_C_1 = H_C_1,
    Omega_C_j = Omega_C_j,
    T = T,
    p1 = nrow(R0),
    j = j_fixed,
    b = b_fixed
  )
  
  
  H_C_gap <-
    norm(
      H_C_1 - H_C_LS,
      type = "F"
    )
  
  
  H_C_orthogonality <-
    norm(
      crossprod(H_C_1) -
        diag(r2),
      type = "F"
    )
  # ---------------------------------------------------------------------------
  # 8. Least-squares rotation
  # ---------------------------------------------------------------------------
  
  H_R_LS <- compute_H_R_LS(
    R0 = R0,
    R_hat_1 = pass1$R_hat_1
  )
  
  
  # ---------------------------------------------------------------------------
  # 9. CLT statistic
  # ---------------------------------------------------------------------------
  
  clt_R <- compute_row_clt_statistic(
    R_hat_1 = pass1$R_hat_1,
    R0 = R0,
    H_R_1 = H_R_1,
    Omega_R_i = Omega_R_i,
    T = T,
    p2 = p2,
    i = i_fixed,
    a = a_fixed
  )
  
  
  # ---------------------------------------------------------------------------
  # 10. Diagnostics
  # ---------------------------------------------------------------------------
  
  H_gap <-
    norm(
      H_R_1 - H_R_LS,
      type = "F"
    )
  
  
  H_orthogonality <-
    norm(
      crossprod(H_R_1) -
        diag(r1),
      type = "F"
    )
  
  
  return(
    list(
      
      Z_R = clt_R$Z,
      Z_R_fixed = clt_R$Z_fixed,
      scaled_error_R = clt_R$scaled_error,
      
      H_R_gap = H_gap,
      H_R_orthogonality = H_orthogonality,
      
      Z_C = clt_C$Z,
      Z_C_fixed = clt_C$Z_fixed,
      scaled_error_C = clt_C$scaled_error,
      
      H_C_gap = H_C_gap,
      H_C_orthogonality = H_C_orthogonality
    )
  )
}


# ==============================================================================
# 13. MONTE CARLO SIMULATION
# ==============================================================================

set.seed(seed_mc)


Z_R_mc <- matrix(
  NA_real_,
  nrow = B,
  ncol = r1
)

scaled_error_R_mc <- matrix(
  NA_real_,
  nrow = B,
  ncol = r1
)

H_R_gap_mc <- numeric(B)
H_R_orthogonality_mc <- numeric(B)


Z_C_mc <- matrix(
  NA_real_,
  nrow = B,
  ncol = r2
)

scaled_error_C_mc <- matrix(
  NA_real_,
  nrow = B,
  ncol = r2
)

H_C_gap_mc <- numeric(B)
H_C_orthogonality_mc <- numeric(B)


for (b in seq_len(B)) {
  
  out_b <- run_one_replication(
    T = T,
    d = d,
    sigma_z = sigma_z,
    W_list = proxy_factor_population$W_list,
    R0 = R0,
    C0 = C0,
    R0_g = R0_g,
    C0_g = C0_g,
    L_Xi_row = L_Xi_row,
    L_Xi_col = L_Xi_col,
    Omega_R_i = theoretical_covariance_R$Omega_R_i,
    Omega_C_j = theoretical_covariance_C$Omega_C_j,
    i_fixed = i_fixed,
    a_fixed = a_fixed,
    j_fixed = j_fixed,
    b_fixed = b_fixed
  )
  
  Z_R_mc[b, ] <- out_b$Z_R
  scaled_error_R_mc[b, ] <- out_b$scaled_error_R
  
  H_R_gap_mc[b] <- out_b$H_R_gap
  H_R_orthogonality_mc[b] <- out_b$H_R_orthogonality
  
  
  Z_C_mc[b, ] <- out_b$Z_C
  scaled_error_C_mc[b, ] <- out_b$scaled_error_C
  
  H_C_gap_mc[b] <- out_b$H_C_gap
  H_C_orthogonality_mc[b] <- out_b$H_C_orthogonality
  
  
  if (
    b == 1 ||
    b %% 100 == 0 ||
    b == B
  ) {
    
    cat(
      "Replication",
      b,
      "of",
      B,
      "\n"
    )
  }
}


colnames(Z_R_mc) <-
  paste0(
    "coordinate_",
    seq_len(r1)
  )

colnames(scaled_error_R_mc) <-
  paste0(
    "coordinate_",
    seq_len(r1)
  )


colnames(Z_C_mc) <-
  paste0(
    "coordinate_",
    seq_len(r2)
  )

colnames(scaled_error_C_mc) <-
  paste0(
    "coordinate_",
    seq_len(r2)
  )


Z_R_fixed_mc <-
  Z_R_mc[, a_fixed]

Z_C_fixed_mc <-
  Z_C_mc[, b_fixed]


# ==============================================================================
# 14. MONTE CARLO RESULTS
# ==============================================================================


process_loading_results <- function(
    Z_mc,
    scaled_error_mc,
    theoretical_covariance_matrix,
    H_gap_mc,
    H_orthogonality_mc,
    path_results,
    loading_name,
    index_fixed,
    coordinate_fixed,
    parameters
) {
  
  Z_mc <- as.matrix(Z_mc)
  scaled_error_mc <- as.matrix(scaled_error_mc)
  theoretical_covariance_matrix <-
    as.matrix(theoretical_covariance_matrix)
  
  r <- ncol(Z_mc)
  
  
  # ============================================================================
  # 14.1 NORMALITY AND COVERAGE
  # ============================================================================
  
  diagnostics_table <- data.frame(
    coordinate = seq_len(r),
    mean_Z      = NA_real_,
    sd_Z        = NA_real_,
    q025_Z      = NA_real_,
    q975_Z      = NA_real_,
    coverage_90 = NA_real_,
    coverage_95 = NA_real_,
    coverage_99 = NA_real_
  )
  
  
  for (a in seq_len(r)) {
    
    Z_a <- Z_mc[, a]
    
    diagnostics_table$mean_Z[a] <-
      mean(Z_a)
    
    diagnostics_table$sd_Z[a] <-
      sd(Z_a)
    
    diagnostics_table$q025_Z[a] <-
      unname(
        quantile(Z_a, 0.025)
      )
    
    diagnostics_table$q975_Z[a] <-
      unname(
        quantile(Z_a, 0.975)
      )
    
    diagnostics_table$coverage_90[a] <-
      mean(
        abs(Z_a) <= qnorm(0.95)
      )
    
    diagnostics_table$coverage_95[a] <-
      mean(
        abs(Z_a) <= qnorm(0.975)
      )
    
    diagnostics_table$coverage_99[a] <-
      mean(
        abs(Z_a) <= qnorm(0.995)
      )
  }
  
  
  cat("\n============================================================\n")
  cat(
    toupper(loading_name),
    "LOADING: NORMALITY AND COVERAGE\n"
  )
  cat("============================================================\n\n")
  
  print(
    round(
      diagnostics_table,
      4
    )
  )
  
  
  # ============================================================================
  # 14.2 EMPIRICAL VERSUS THEORETICAL COVARIANCE
  # ============================================================================
  
  empirical_covariance <-
    cov(
      scaled_error_mc
    )
  
  
  covariance_gap <-
    empirical_covariance -
    theoretical_covariance_matrix
  
  
  covariance_gap_F <-
    norm(
      covariance_gap,
      type = "F"
    )
  
  
  relative_covariance_gap_F <-
    covariance_gap_F /
    norm(
      theoretical_covariance_matrix,
      type = "F"
    )
  
  
  cat("\n============================================================\n")
  cat(
    toupper(loading_name),
    "LOADING: EMPIRICAL COVARIANCE\n"
  )
  cat("============================================================\n\n")
  
  print(
    round(
      empirical_covariance,
      5
    )
  )
  
  
  cat("\nTheoretical covariance:\n\n")
  
  print(
    round(
      theoretical_covariance_matrix,
      5
    )
  )
  
  
  cat(
    "\nFrobenius gap          :",
    round(covariance_gap_F, 5),
    "\n"
  )
  
  cat(
    "Relative Frobenius gap :",
    round(relative_covariance_gap_F, 5),
    "\n"
  )
  
  
  # ============================================================================
  # 14.3 COVARIANCE OF STANDARDIZED STATISTICS
  # ============================================================================
  
  standardized_covariance <-
    cov(
      Z_mc
    )
  
  
  standardized_covariance_gap_F <-
    norm(
      standardized_covariance -
        diag(r),
      type = "F"
    )
  
  
  cat("\n============================================================\n")
  cat(
    toupper(loading_name),
    "LOADING: COVARIANCE OF STANDARDIZED ERRORS\n"
  )
  cat("Target: I\n")
  cat("============================================================\n\n")
  
  print(
    round(
      standardized_covariance,
      5
    )
  )
  
  cat(
    "\nDistance from identity:",
    round(
      standardized_covariance_gap_F,
      5
    ),
    "\n"
  )
  
  
  # ============================================================================
  # 14.4 ROTATION DIAGNOSTICS
  # ============================================================================
  
  rotation_summary <- data.frame(
    mean_H_gap =
      mean(H_gap_mc),
    
    mean_H_orthogonality =
      mean(H_orthogonality_mc)
  )
  
  
  cat("\n============================================================\n")
  cat(
    toupper(loading_name),
    "LOADING: ROTATION DIAGNOSTICS\n"
  )
  cat("============================================================\n\n")
  
  print(
    round(
      rotation_summary,
      5
    )
  )
  
  
  # ============================================================================
  # 15. SAVE NUMERICAL RESULTS
  # ============================================================================
  
  write.csv(
    diagnostics_table,
    file = file.path(
      path_results,
      "MC_diagnostics.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    empirical_covariance,
    file = file.path(
      path_results,
      "empirical_covariance.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    theoretical_covariance_matrix,
    file = file.path(
      path_results,
      "theoretical_covariance.csv"
    ),
    row.names = FALSE
  )
  
  
  write.csv(
    standardized_covariance,
    file = file.path(
      path_results,
      "standardized_covariance.csv"
    ),
    row.names = FALSE
  )
  
  
  saveRDS(
    list(
      
      parameters =
        parameters,
      
      diagnostics =
        diagnostics_table,
      
      empirical_covariance =
        empirical_covariance,
      
      theoretical_covariance =
        theoretical_covariance_matrix,
      
      standardized_covariance =
        standardized_covariance,
      
      covariance_gap_F =
        covariance_gap_F,
      
      relative_covariance_gap_F =
        relative_covariance_gap_F,
      
      standardized_covariance_gap_F =
        standardized_covariance_gap_F,
      
      H_gap_mc =
        H_gap_mc,
      
      H_orthogonality_mc =
        H_orthogonality_mc,
      
      Z_mc =
        Z_mc,
      
      scaled_error_mc =
        scaled_error_mc
      
    ),
    
    file = file.path(
      path_results,
      "MC_asymptotic_normality.rds"
    )
  )
  
  
  # ============================================================================
  # 16. FIGURES
  # ============================================================================
  
  Z_fixed_mc <-
    Z_mc[, coordinate_fixed]
  
  
  mean_Z <-
    mean(Z_fixed_mc)
  
  sd_Z <-
    sd(Z_fixed_mc)
  
  
  coverage_95 <-
    mean(
      abs(Z_fixed_mc) <=
        qnorm(0.975)
    )
  
  
  plot_data <- data.frame(
    Z = Z_fixed_mc
  )
  
  
  theme_mc <- theme_minimal(
    base_size = 13
  ) +
    theme(
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15,
        margin = margin(b = 7)
      ),
      
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11.5,
        margin = margin(b = 14)
      ),
      
      axis.title = element_text(
        size = 12
      ),
      
      axis.text = element_text(
        size = 10.5
      ),
      
      panel.grid.minor =
        element_blank(),
      
      panel.grid.major =
        element_line(
          color = "grey90",
          linewidth = 0.35
        ),
      
      plot.margin =
        margin(
          12, 16, 12, 12
        )
    )
  
  
  # ---------------------------------------------------------------------------
  # Histogram
  # ---------------------------------------------------------------------------
  
  p_hist <- ggplot(
    plot_data,
    aes(x = Z)
  ) +
    
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      bins = 100,
      fill = "#B8C4CE",
      color = "white",
      linewidth = 0.35
    ) +
    
    stat_function(
      fun = dnorm,
      args = list(
        mean = 0,
        sd = 1
      ),
      color = "#17365D",
      linewidth = 1.15
    ) +
    
    geom_vline(
      xintercept = 0,
      color = "grey45",
      linewidth = 0.45
    ) +
    
    geom_vline(
      xintercept = c(
        qnorm(0.025),
        qnorm(0.975)
      ),
      color = "#A63D40",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    
    labs(
      title = paste0(
        "Monte Carlo Distribution of the Standardized ",
        loading_name,
        "-Loading Error"
      ),
      
      subtitle = paste0(
        "Mean = ",
        sprintf("%.3f", mean_Z),
        "    |    SD = ",
        sprintf("%.3f", sd_Z),
        "    |    95% coverage = ",
        sprintf(
          "%.1f",
          100 * coverage_95
        ),
        "%"
      ),
      
      x = "Standardized estimation error",
      y = "Density"
    ) +
    
    theme_mc + ggplot2::labs(title = NULL, subtitle = NULL)
  
  
  print(p_hist)
  
  
  ggsave(
    filename = file.path(
      path_results,
      "histogram_standardized_error.png"
    ),
    plot = p_hist,
    width = 8,
    height = 5.5,
    dpi = 320,
    bg = "white"
  )
  
  
  ggsave(
    filename = file.path(
      path_results,
      "histogram_standardized_error.pdf"
    ),
    plot = p_hist,
    width = 8,
    height = 5.5,
    bg = "white"
  )
  
  
  # ---------------------------------------------------------------------------
  # Q-Q plot
  # ---------------------------------------------------------------------------
  
  p_qq <- ggplot(
    plot_data,
    aes(sample = Z)
  ) +
    
    stat_qq(
      color = "#394B59",
      size = 1.5,
      alpha = 0.70
    ) +
    
    stat_qq_line(
      color = "#A63D40",
      linewidth = 0.95
    ) +
    
    labs(
      title = paste0(
        "Normal Q-Q Plot of the Standardized ",
        loading_name,
        "-Loading Error"
      ),
      
      subtitle = paste0(
        loading_name,
        " ", index_fixed,
        ", coordinate ",
        coordinate_fixed,
        "    |    ",
        nrow(Z_mc),
        " Monte Carlo replications"
      ),
      
      x = "Theoretical normal quantiles",
      y = "Monte Carlo quantiles"
    ) +
    
    theme_mc + ggplot2::labs(title = NULL, subtitle = NULL)
  
  
  print(p_qq)
  
  
  ggsave(
    filename = file.path(
      path_results,
      "qqplot_standardized_error.png"
    ),
    plot = p_qq,
    width = 8,
    height = 5.5,
    dpi = 320,
    bg = "white"
  )
  
  
  ggsave(
    filename = file.path(
      path_results,
      "qqplot_standardized_error.pdf"
    ),
    plot = p_qq,
    width = 8,
    height = 5.5,
    bg = "white"
  )
  
  
  return(
    list(
      diagnostics =
        diagnostics_table,
      
      empirical_covariance =
        empirical_covariance,
      
      theoretical_covariance =
        theoretical_covariance_matrix,
      
      standardized_covariance =
        standardized_covariance,
      
      covariance_gap_F =
        covariance_gap_F,
      
      relative_covariance_gap_F =
        relative_covariance_gap_F,
      
      rotation_summary =
        rotation_summary
    )
  )
}

mc_parameters <- list(
  B = B,
  T = T,
  p1 = p1,
  p2 = p2,
  r1 = r1,
  r2 = r2,
  L = L,
  g1 = g1,
  g2 = g2,
  i_fixed = i_fixed,
  a_fixed = a_fixed,
  j_fixed = j_fixed,
  b_fixed = b_fixed,
  sigma_z = sigma_z,
  delta_r = delta_r,
  delta_c = delta_c,
  d = d
)

results_R <- process_loading_results(
  
  Z_mc =
    Z_R_mc,
  
  scaled_error_mc =
    scaled_error_R_mc,
  
  theoretical_covariance_matrix =
    theoretical_covariance_R$Omega_R_i,
  
  H_gap_mc =
    H_R_gap_mc,
  
  H_orthogonality_mc =
    H_R_orthogonality_mc,
  
  path_results =
    path_results_R,
  
  loading_name =
    "Row",
  
  index_fixed =
    i_fixed,
  
  coordinate_fixed =
    a_fixed,
  
  parameters =
    mc_parameters
)

results_C <- process_loading_results(
  
  Z_mc =
    Z_C_mc,
  
  scaled_error_mc =
    scaled_error_C_mc,
  
  theoretical_covariance_matrix =
    theoretical_covariance_C$Omega_C_j,
  
  H_gap_mc =
    H_C_gap_mc,
  
  H_orthogonality_mc =
    H_C_orthogonality_mc,
  
  path_results =
    path_results_C,
  
  loading_name =
    "Column",
  
  index_fixed =
    j_fixed,
  
  coordinate_fixed =
    b_fixed,
  
  parameters =
    mc_parameters
)
