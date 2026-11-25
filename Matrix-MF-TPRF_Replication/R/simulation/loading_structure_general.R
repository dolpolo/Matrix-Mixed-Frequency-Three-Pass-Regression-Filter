# Author: Davide Delfino
# Input: see function arguments; source from a model main.
# Output: function return values; no files written unless documented.
# Source: original functions_sim/loading_structure_general.R; numerical implementation retained.

# =============================================================================
# MATRIX MF-TPRF VS VECTOR MF-TPRF: LOADING-SPACE DIAGNOSTIC
# =============================================================================
# Supported predictive factor dimensions: 2x1, 1x2, 2x2.
# Source the ORIGINAL functions/functions_sim/matrix.dgp.R before this file.
#
# This is a matched Pass-1 experiment (NO nowcasting, NO imputation):
#   Matrix: targeted row/column spectral moments + bilateral projections.
#   Vector: unrestricted rank-K PCA of proxy-projected vec(X).
# Both methods use exactly the same panel, proxy basis and estimation dates.
#
# The ORIGINAL DGP is not modified. We draw it with r1=r2=1 and extend the
# existing predictive part, keeping its G, Xi, release specifications and F11.
# =============================================================================

loading_factor_labels <- function(r = c(2L, 1L)) {
  r <- as.integer(r)
  if (length(r) != 2L || anyNA(r) || !all(r %in% c(1L, 2L)) ||
      prod(r) < 2L) {
    stop("Allowed dimensions are r = c(2,1), c(1,2), or c(2,2).")
  }
  # Countries (row factors) vary fastest under column-major vectorization.
  unlist(lapply(seq_len(r[2L]), function(b) {
    paste0("f", seq_len(r[1L]), b)
  }), use.names = FALSE)
}

loading_default_rho <- function(r, rho_f1) {
  # Separable matrix AR: A_F diagonal and B_F diagonal.
  # For (2,2), the vectorized order is (.90,.60,.45,.30) when rho_f1=.90.
  row_ratios <- c(1, 2 / 3)[seq_len(r[1L])]
  col_ratios <- c(1, 1 / 2)[seq_len(r[2L])]
  as.vector(rho_f1 * outer(row_ratios, col_ratios))
}


# Extend the baseline 1x1 simulation to an r1 x r2 predictive factor.
simulate_loading_dgp <- function(params, r = c(2L, 1L),
                                 rhoF_factors = NULL,
                                 target_weights = NULL) {
  labels <- loading_factor_labels(r)
  r <- as.integer(r)
  K <- prod(r)
  needed <- c("simulate_matrix_dgp_truth", "build_random_loading_block",
              "mariano_murasawa_agg", "aggregate_lf_panel")
  missing <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop("Source the original matrix.dgp.R first. Missing: ",
         paste(missing, collapse = ", "))
  }
  if (params$r1 != 1L || params$r2 != 1L) {
    stop("Supply baseline params from make_dgp_params(r1=1,r2=1). ",
         "The requested dimensions are set ONLY by argument r.")
  }
  
  if (is.null(rhoF_factors)) {
    rhoF_factors <- loading_default_rho(r, params$rhoF)
  }
  if (is.null(target_weights)) target_weights <- rep(1, K)
  if (length(rhoF_factors) != K || any(!is.finite(rhoF_factors)) ||
      any(rhoF_factors <= 0) || any(rhoF_factors >= 1)) {
    stop("rhoF_factors must have K=", K,
         " distinct positive stationary AR(1) coefficients.")
  }
  if (abs(rhoF_factors[1L] - params$rhoF) > 1e-12) {
    stop("The first predictive persistence must equal baseline params$rhoF.")
  }
  rho_distances <- abs(outer(rhoF_factors, rhoF_factors, "-"))
  diag(rho_distances) <- Inf
  if (any(rho_distances < 1e-8)) {
    stop("The K factor persistences must be distinct when using GDP lags.")
  }
  if (all(r == c(2L, 2L)) &&
      abs(rhoF_factors[1L] * rhoF_factors[4L] -
          rhoF_factors[2L] * rhoF_factors[3L]) > 1e-8) {
    stop("For a matrix AR with diagonal A_F and B_F, ",
         "rho11*rho22 must equal rho21*rho12.")
  }
  if (length(target_weights) != K ||
      any(!is.finite(target_weights)) || any(target_weights == 0)) {
    stop("target_weights must contain K=", K, " nonzero finite values.")
  }
  
  sim <- simulate_matrix_dgp_truth(params)
  Tm <- dim(sim$truth$X_monthly_raw)[1L]
  p1 <- params$p1
  p2 <- params$p2
  R1 <- sim$loadings$R_signal
  C1 <- sim$loadings$C_signal
  Rg <- sim$loadings$R_noise
  Cg <- sim$loadings$C_noise
  if (r[1L] == 2L && p1 <= ncol(Rg) + 1L) {
    stop("Need p1 > 1 + v1 to construct two orthogonal predictive rows.")
  }
  if (r[2L] == 2L && p2 <= ncol(Cg) + 1L) {
    stop("Need p2 > 1 + v2 to construct two orthogonal predictive columns.")
  }
  
  # Keep the original first loading; add new orthogonal predictive columns.
  R <- R1
  if (r[1L] == 2L) {
    R <- cbind(R1, build_random_loading_block(
      n_rows = p1, n_cols = 1L, scale = p1,
      B_ref = cbind(R1, Rg)))
  }
  C <- C1
  if (r[2L] == 2L) {
    C <- cbind(C1, build_random_loading_block(
      n_rows = p2, n_cols = 1L, scale = p2,
      B_ref = cbind(C1, Cg)))
  }
  if (!all(dim(R) == c(p1, r[1L])) ||
      !all(dim(C) == c(p2, r[2L]))) {
    stop("Unexpected predictive loading dimensions.")
  }
  
  # F11 is exactly the original simulated factor; the other K-1 processes
  # have independent N(0,1) innovations and factor-specific AR coefficients.
  factor_time <- matrix(NA_real_, nrow = Tm, ncol = K)
  factor_time[, 1L] <- as.numeric(sim$factors$F_signal[, 1L, 1L])
  TT <- Tm + params$burn
  for (k in seq.int(2L, K)) {
    path <- numeric(TT)
    for (tt in 2:TT) {
      path[tt] <- rhoF_factors[k] * path[tt - 1L] + stats::rnorm(1L)
    }
    factor_time[, k] <- tail(path, Tm)
  }
  F <- array(NA_real_, dim = c(Tm, r[1L], r[2L]))
  for (k in seq_len(K)) {
    a <- (k - 1L) %% r[1L] + 1L
    b <- (k - 1L) %/% r[1L] + 1L
    F[, a, b] <- factor_time[, k]
  }
  
  # Original G and Xi remain unchanged: add only the extra predictive maps.
  X <- sim$truth$X_monthly_raw
  for (k in seq.int(2L, K)) {
    a <- (k - 1L) %% r[1L] + 1L
    b <- (k - 1L) %/% r[1L] + 1L
    map <- tcrossprod(R[, a], C[, b])
    for (tt in seq_len(Tm)) {
      X[tt, , ] <- X[tt, , ] + factor_time[tt, k] * map
    }
  }
  
  # The target and EA GDP proxy depend on ALL predictive factors.
  # The original country-specific target noise is preserved up to a common
  # variance rescaling to maintain R2_BI = params$R2_BI.
  common <- as.numeric(factor_time %*% target_weights)
  sigma_y <- sqrt(stats::var(common) *
                    (1 - params$R2_BI) / params$R2_BI)
  baseline_noise <- sim$truth$y_monthly - sim$truth$y_monthly_signal
  if (!is.finite(sim$params$sigma_y) || sim$params$sigma_y <= 0) {
    stop("Invalid baseline sigma_y.")
  }
  y_signal <- matrix(common, nrow = Tm, ncol = p1)
  y_monthly <- y_signal +
    baseline_noise * (sigma_y / sim$params$sigma_y)
  y_q <- mariano_murasawa_agg(y_monthly)
  y_q_signal <- mariano_murasawa_agg(y_signal)
  w <- params$w_ea
  if (is.null(w)) w <- rep(1 / p1, p1)
  z_q_signal <- as.numeric(y_q_signal %*% w)
  # Preserve the actual N(0,sigma_z^2) proxy errors of the original draw.
  z_noise <- sim$truth$z_quarterly - sim$truth$z_quarterly_signal
  z_q <- z_q_signal + z_noise
  
  # Recompute quarterly series AFTER burn-in so X and target share q indexes.
  sim$loadings$R_signal <- R
  sim$loadings$C_signal <- C
  sim$factors$F_signal <- F
  sim$factors$f_true <- common
  sim$truth$X_monthly_raw <- X
  sim$truth$X_quarterly_lf_raw <- aggregate_lf_panel(X, params$var_type)
  sim$truth$y_monthly <- y_monthly
  sim$truth$y_monthly_signal <- y_signal
  sim$truth$y_quarterly <- y_q
  sim$truth$y_quarterly_signal <- y_q_signal
  sim$truth$z_quarterly <- z_q
  sim$truth$z_quarterly_signal <- z_q_signal
  sim$params$r1 <- r[1L]
  sim$params$r2 <- r[2L]
  sim$params$rhoF_factors <- rhoF_factors
  A_diag <- rhoF_factors[seq_len(r[1L])] /
    sqrt(rhoF_factors[1L])
  B_diag <- rhoF_factors[seq.int(1L, by = r[1L], length.out = r[2L])] /
    sqrt(rhoF_factors[1L])
  sim$factor_dynamics <- list(A_F = diag(A_diag, r[1L], r[1L]),
                              B_F = diag(B_diag, r[2L], r[2L]))
  sim$params$target_weights <- target_weights
  sim$params$sigma_y <- sigma_y
  sim$factor_labels <- labels
  sim
}


# Rotate an orthonormal estimated space towards another orthonormal space.
# This is only for interpretable plots; principal angles do not use rotation.
align_loading_basis <- function(Q_est, Q_ref) {
  S <- svd(crossprod(Q_est, Q_ref))
  Q_est %*% (S$u %*% t(S$v))
}

loading_subspace_diagnostics <- function(Q_est, Q_ref) {
  sv <- svd(crossprod(Q_est, Q_ref), nu = 0L, nv = 0L)$d
  sv <- pmax(0, pmin(1, sv))
  list(angles = acos(sv) * 180 / pi,
       mean_cos2 = mean(sv^2),
       projection_distance = sqrt(max(0, 2 * length(sv) -
                                        2 * sum(sv^2))))
}


# Targeted Pass-1 loading extraction from the original mixed quarterly panel.
# Flow variables use quarterly SUMS; stock variables use quarterly MEANS.
# Thus the true quarterly column loading C_q equals D*C, with D_jj=3 for
# flows and D_jj=1 for stocks. This preserves the low-rank matrix structure.
estimate_loading_structure <- function(sim, Lproxy = NULL,
                                       projected_updates = 1L) {
  r <- as.integer(c(sim$params$r1, sim$params$r2))
  K <- prod(r)
  if (K < 2L || K > 4L || any(r > 2L)) {
    stop("Supported predictive ranks: 2x1, 1x2, 2x2.")
  }
  if (is.null(Lproxy)) Lproxy <- max(4L, K)
  if (length(Lproxy) != 1L || !is.finite(Lproxy) ||
      Lproxy < K || Lproxy != as.integer(Lproxy)) {
    stop("Lproxy must be an integer >= K = ", K)
  }
  Lproxy <- as.integer(Lproxy)
  if (length(projected_updates) != 1L || !is.finite(projected_updates) ||
      projected_updates < 0 || projected_updates != as.integer(projected_updates)) {
    stop("projected_updates must be a nonnegative integer.")
  }
  projected_updates <- as.integer(projected_updates)
  
  Xq <- sim$truth$X_quarterly_lf_raw
  d <- dim(Xq)
  if (length(d) != 3L) stop("Missing complete quarterly predictor tensor.")
  Tq <- d[1L]
  p1 <- d[2L]
  p2 <- d[3L]
  if (any(!is.finite(Xq))) stop("This diagnostic requires true completed Xq.")
  z <- sim$truth$z_quarterly
  if (length(z) != Tq || Lproxy >= Tq - 10L) {
    stop("GDP proxy and Xq are not aligned, or T is too short for Lproxy.")
  }
  
  # R vec(X) convention: countries vary fastest; C %x% R.
  X <- t(matrix(aperm(Xq, c(2L, 3L, 1L)),
                nrow = p1 * p2, ncol = Tq))
  Zraw <- sapply(0:(Lproxy - 1L), function(h) {
    if (h == 0L) z else c(rep(NA_real_, h), head(z, -h))
  })
  Zraw <- matrix(Zraw, nrow = Tq, ncol = Lproxy)
  keep <- stats::complete.cases(Zraw)
  if (sum(keep) <= Lproxy + 10L) {
    stop("Insufficient complete proxy observations.")
  }
  X <- X[keep, , drop = FALSE]
  Z <- Zraw[keep, , drop = FALSE]
  X <- sweep(X, 2L, colMeans(X), "-")
  Z <- scale(Z, center = TRUE, scale = TRUE)
  if (any(!is.finite(Z))) stop("The GDP proxy has a constant lag column.")
  zqr <- qr(Z)
  if (zqr$rank < Lproxy) stop("GDP proxy lags are linearly dependent.")
  Qz <- qr.Q(zqr, complete = FALSE)[, seq_len(Lproxy), drop = FALSE]
  
  # Same proxy projection for both estimators.
  XP <- Qz %*% crossprod(Qz, X)
  sv <- svd(XP, nu = 0L, nv = K)
  if (sv$d[K] <= 1e-8 * sv$d[1L]) {
    stop("The proxy-projected panel has rank less than K=", K)
  }
  strength <- sv$d[K] / sv$d[1L]
  if (strength < 0.05) {
    warning("Weak K-th proxy direction: d[K]/d[1] = ",
            signif(strength, 3),
            ". Consider larger T or a richer proxy system.")
  }
  Qvec <- sv$v[, seq_len(K), drop = FALSE]
  
  # Matrix spectral moments, initialization and bilateral refinement.
  M <- lapply(seq_len(Lproxy), function(ell) {
    matrix(as.numeric(crossprod(Qz[, ell], X)), nrow = p1, ncol = p2)
  })
  eig_top <- function(S, kk) {
    S <- (S + t(S)) / 2
    eigen(S, symmetric = TRUE)$vectors[, seq_len(kk), drop = FALSE]
  }
  SR <- matrix(0, p1, p1)
  SC <- matrix(0, p2, p2)
  for (Mm in M) {
    SR <- SR + tcrossprod(Mm)
    SC <- SC + crossprod(Mm)
  }
  Rhat <- eig_top(SR, r[1L])
  Chat <- eig_top(SC, r[2L])
  if (projected_updates > 0L) {
    for (iter in seq_len(projected_updates)) {
      PC <- tcrossprod(Chat)
      SR <- matrix(0, p1, p1)
      for (Mm in M) SR <- SR + Mm %*% PC %*% t(Mm)
      Rhat <- eig_top(SR, r[1L])
      
      PR <- tcrossprod(Rhat)
      SC <- matrix(0, p2, p2)
      for (Mm in M) SC <- SC + t(Mm) %*% PR %*% Mm
      Chat <- eig_top(SC, r[2L])
    }
  }
  Qmat <- kronecker(Chat, Rhat)
  
  # Population loading space on THE SAME mixed-frequency quarterly scale.
  Rtrue <- sim$loadings$R_signal
  Ctrue <- sim$loadings$C_signal
  flow <- grepl("_flow$", sim$params$var_type)
  Ctrue[flow, ] <- 3 * Ctrue[flow, , drop = FALSE]
  Qrtrue <- qr.Q(qr(Rtrue))[, seq_len(r[1L]), drop = FALSE]
  Qctrue <- qr.Q(qr(Ctrue))[, seq_len(r[2L]), drop = FALSE]
  Qtrue <- kronecker(Qctrue, Qrtrue)
  
  # CRUCIAL for r=(2,2): align row and column spaces SEPARATELY.
  # A generic KxK Procrustes rotation would destroy rank-one Kronecker maps.
  Raligned <- align_loading_basis(Rhat, Qrtrue)
  Caligned <- align_loading_basis(Chat, Qctrue)
  Qmat_aligned <- kronecker(Caligned, Raligned)
  Qvec_aligned <- align_loading_basis(Qvec, Qtrue)
  
  dm <- loading_subspace_diagnostics(Qmat, Qtrue)
  dv <- loading_subspace_diagnostics(Qvec, Qtrue)
  dmv <- loading_subspace_diagnostics(Qmat, Qvec)
  metrics <- data.frame(
    comparison = c("Matrix vs true", "Vector vs true", "Matrix vs Vector"),
    mean_cos2 = c(dm$mean_cos2, dv$mean_cos2, dmv$mean_cos2),
    projection_distance = c(dm$projection_distance,
                            dv$projection_distance,
                            dmv$projection_distance),
    max_principal_angle_deg = c(max(dm$angles), max(dv$angles),
                                max(dmv$angles))
  )
  all_angles <- rbind(dm$angles, dv$angles, dmv$angles)
  colnames(all_angles) <- paste0("angle", seq_len(K), "_degrees")
  angles <- cbind(comparison = metrics$comparison,
                  as.data.frame(all_angles, check.names = FALSE))
  
  # Basis-dependent descriptive separability for the displayed maps.
  rank_one_fraction <- function(Q, k) {
    A <- matrix(Q[, k], nrow = p1, ncol = p2)
    d <- svd(A, nu = 0L, nv = 0L)$d
    d[1L]^2 / sum(d^2)
  }
  separability <- data.frame(
    component = loading_factor_labels(r),
    true_rank1_share = vapply(seq_len(K), function(k)
      rank_one_fraction(Qtrue, k), numeric(1L)),
    matrix_rank1_share = vapply(seq_len(K), function(k)
      rank_one_fraction(Qmat_aligned, k), numeric(1L)),
    vector_rank1_share = vapply(seq_len(K), function(k)
      rank_one_fraction(Qvec_aligned, k), numeric(1L))
  )
  list(Q_true = Qtrue, Q_matrix = Qmat, Q_vector = Qvec,
       Q_matrix_aligned = Qmat_aligned, Q_vector_aligned = Qvec_aligned,
       R_hat = Rhat, C_hat = Chat, metrics = metrics, angles = angles,
       separability = separability, labels = loading_factor_labels(r),
       p1 = p1, p2 = p2, r = r, K = K,
       n_monthly_vars = sim$params$Nm, n_used = sum(keep),
       proxy_strength = strength)
}


plot_loading_structure <- function(fit, filename = NULL,
                                   show_truth = FALSE,
                                   show_indicators = FALSE) {
  # Scientific-brackets layout: no overall title or country-row labels.
  # Factor headings sit inside full-panel-width grey rectangles.
  # The PDF has vector text/brackets; only the heatmaps are rasterized.
  p1 <- as.integer(fit$p1)
  p2 <- as.integer(fit$p2)
  K <- as.integer(fit$K)
  Nm <- as.integer(fit$n_monthly_vars)
  if (length(p1) != 1L || length(p2) != 1L ||
      length(K) != 1L || anyNA(c(p1, p2, K)) ||
      p1 < 1L || p2 < 2L || !(K %in% c(2L, 4L))) {
    stop("Invalid panel dimensions or predictive factor count in fit.")
  }
  if (length(Nm) != 1L || is.na(Nm) || Nm < 1L || Nm >= p2) {
    stop("fit$n_monthly_vars must lie strictly between 0 and p2.")
  }
  if (length(fit$labels) != K ||
      any(!grepl("^f[0-9]+$", fit$labels))) {
    stop("Expected factor labels such as f11, f21, f12, f22.")
  }
  
  groups <- list()
  if (show_truth) groups[["True loadings"]] <- fit$Q_true
  groups[["Matrix MF-TPRF"]] <- fit$Q_matrix_aligned
  groups[["Vector MF-TPRF"]] <- fit$Q_vector_aligned
  n_models <- length(groups)
  
  # Restore the p1 x p2 loading maps and normalize each factor separately.
  # Crucially, both estimators retain their existing alignment to the truth.
  maps <- lapply(groups, function(Q) {
    Q <- as.matrix(Q)
    if (!identical(dim(Q), c(p1 * p2, K)) || any(!is.finite(Q))) {
      stop("A loading matrix has unexpected dimensions or non-finite values.")
    }
    lapply(seq_len(K), function(k) {
      M <- matrix(Q[, k], nrow = p1, ncol = p2)
      rms <- sqrt(mean(M^2))
      if (!is.finite(rms) || rms <= 0) stop("A loading map is degenerate.")
      M / rms
    })
  })
  all_values <- unlist(lapply(maps, unlist), use.names = FALSE)
  lim <- as.numeric(stats::quantile(abs(all_values), 0.98, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  
  # One symmetric colour scale for every factor and both estimators.
  palette <- grDevices::colorRampPalette(c(
    "#306C9E", "#9DB8D2", "#F9F8F6", "#F0A29A", "#C95055"
  ))(257L)
  as_colour_raster <- function(M) {
    z <- pmax(-1, pmin(1, M / lim))
    idx <- 1L + as.integer(round((z + 1) * 128))
    as.raster(matrix(palette[idx], nrow = nrow(M), ncol = ncol(M)))
  }
  
  ink <- "#1D2E3E"
  faint <- "#DCE2E5"
  bracket_colour <- "#7C8D98"
  header_fill <- "#EEEEEE"
  panel_left <- 0.027
  panel_right <- 0.940
  panel_gap <- if (K == 4L) 0.010 else 0.015
  panel_width <- (panel_right - panel_left - (K - 1L) * panel_gap) / K
  panel_x <- panel_left + (seq_len(K) - 1L) * (panel_width + panel_gap)
  
  if (n_models == 3L) {
    panel_height <- 0.204
    panel_bottom <- c(0.608, 0.379, 0.150)
    legend_y <- 0.103
    tick_y <- 0.074
  } else {
    panel_height <- 0.295
    panel_bottom <- c(0.518, 0.195)
    legend_y <- 0.132
    tick_y <- 0.099
  }
  
  pieces <- list()
  add <- function(grob) pieces[[length(pieces) + 1L]] <<- grob
  add_line <- function(x0, y0, x1, y1, colour = bracket_colour,
                       linewidth = 0.75) {
    add(grid::segmentsGrob(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                           default.units = "npc",
                           gp = grid::gpar(col = colour, lwd = linewidth)))
  }
  
  for (k in seq_len(K)) {
    left <- panel_x[k]
    right <- left + panel_width
    divider <- left + panel_width * Nm / p2
    
    # Grey heading exactly as wide as the heatmap beneath it.
    header_y <- 0.948
    header_height <- 0.062
    add(grid::rectGrob(x = (left + right) / 2, y = header_y,
                       width = panel_width, height = header_height,
                       default.units = "npc",
                       gp = grid::gpar(fill = header_fill, col = NA)))
    
    # Plotmath indices: f_{11}, f_{21}, f_{12}, f_{22}.
    subscript <- substring(fit$labels[k], 2L)
    factor_label <- parse(text = paste0("italic(f)[", subscript, "]"))
    add(grid::textGrob(factor_label, x = (left + right) / 2, y = header_y,
                       default.units = "npc",
                       gp = grid::gpar(col = ink, fontsize = 16,
                                       fontface = "bold")))
    
    # Two scientific brackets placed ABOVE the first row of heatmaps.
    for (ends in list(c(left, divider), c(divider, right))) {
      a <- ends[1L] + 0.0012
      b <- ends[2L] - 0.0012
      add_line(a, 0.834, b, 0.834, linewidth = 0.85)
      add_line(a, 0.834, a, 0.823, linewidth = 0.85)
      add_line(b, 0.834, b, 0.823, linewidth = 0.85)
    }
    add(grid::textGrob("MONTHLY", x = (left + divider) / 2, y = 0.866,
                       default.units = "npc",
                       gp = grid::gpar(col = ink, fontsize = 9, fontface = "bold")))
    add(grid::textGrob("QUARTERLY", x = (divider + right) / 2, y = 0.866,
                       default.units = "npc",
                       gp = grid::gpar(col = ink, fontsize = 8, fontface = "bold")))
    
    for (j in seq_len(n_models)) {
      bottom <- panel_bottom[j]
      M <- maps[[j]][[k]]
      add(grid::rasterGrob(
        image = as_colour_raster(M),
        x = left + panel_width / 2, y = bottom + panel_height / 2,
        width = panel_width, height = panel_height,
        default.units = "npc", interpolate = FALSE))
      
      # Subtle country-row divisions: no row names or y-axis ticks.
      if (p1 > 1L) {
        row_levels <- bottom + panel_height * (seq_len(p1 - 1L) / p1)
        add(grid::segmentsGrob(
          x0 = rep(left, length(row_levels)),
          y0 = row_levels,
          x1 = rep(right, length(row_levels)),
          y1 = row_levels,
          default.units = "npc",
          gp = grid::gpar(col = grDevices::adjustcolor("white", 0.35),
                          lwd = 0.35)))
      }
      add_line(divider, bottom, divider, bottom + panel_height,
               colour = "#46545F", linewidth = 0.85)
      
      if (show_indicators) {
        step <- max(1L, ceiling(p2 / 10L))
        ticks <- unique(c(1L, seq.int(step, p2, by = step), p2))
        for (i in ticks) {
          add(grid::textGrob(as.character(i),
                             x = left + panel_width * (i - 0.5) / p2,
                             y = bottom - 0.012, default.units = "npc",
                             gp = grid::gpar(col = "#6A7884", fontsize = 6.2)))
        }
      }
    }
  }
  
  # Model names vertically on the RIGHT, as in the original loading figure.
  for (j in seq_len(n_models)) {
    mid <- panel_bottom[j] + panel_height / 2
    add(grid::rectGrob(x = 0.968, y = mid, width = 0.037,
                       height = panel_height, default.units = "npc",
                       gp = grid::gpar(col = NA, fill = "#F0F2F4")))
    add(grid::textGrob(names(groups)[j], x = 0.968, y = mid,
                       rot = 270, default.units = "npc",
                       gp = grid::gpar(col = ink, fontsize = 10.5,
                                       fontface = "bold")))
  }
  
  # Compact, hairline-height horizontal legend directly below the maps.
  legend_width <- if (K == 4L) 0.205 else 0.260
  legend_left <- (1 - legend_width) / 2
  add(grid::rasterGrob(as.raster(matrix(palette, nrow = 1L)),
                       x = 0.5, y = legend_y,
                       width = legend_width, height = 0.014,
                       default.units = "npc", interpolate = TRUE))
  add(grid::rectGrob(x = 0.5, y = legend_y,
                     width = legend_width, height = 0.014,
                     default.units = "npc",
                     gp = grid::gpar(col = faint, fill = NA, lwd = 0.55)))
  legend_breaks <- pretty(c(-lim, lim), n = 5L)
  legend_breaks <- legend_breaks[abs(legend_breaks) <= lim + 1e-10]
  if (length(legend_breaks) < 3L) legend_breaks <- c(-lim, 0, lim)
  for (value in legend_breaks) {
    xpos <- legend_left + legend_width * (value / lim + 1) / 2
    add(grid::textGrob(formatC(value, digits = 2, format = "fg"),
                       x = xpos, y = tick_y,
                       default.units = "npc",
                       gp = grid::gpar(col = ink, fontsize = 9.5)))
  }
  
  figure <- grid::gTree(children = do.call(grid::gList, pieces),
                        name = "scientific_brackets_loading_maps")
  
  if (is.null(filename)) {
    grid::grid.newpage()
    grid::grid.draw(figure)
  } else {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    width_in <- 3.7 * K + 1.2
    height_in <- if (show_truth) 7.5 else 5.8
    ext <- tolower(tools::file_ext(filename))
    if (ext == "png") {
      grDevices::png(filename, width = width_in, height = height_in,
                     units = "in", res = 300, bg = "white")
    } else if (ext == "pdf") {
      grDevices::pdf(filename, width = width_in, height = height_in,
                     useDingbats = FALSE, onefile = FALSE)
    } else {
      stop("Use a .png or .pdf filename for loading figures.")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(figure)
  }
  invisible(figure)
}


# Main callable function: ONE DGP draw, TWO fits, diagnostic figures and CSVs.
run_loading_structure_comparison <- function(
    dgp_params, r = c(2L, 1L),
    out_dir = "MonteCarlo_simulation/results/loading_structure",
    rhoF_factors = NULL, target_weights = NULL,
    Lproxy = NULL, projected_updates = 1L) {
  r <- as.integer(r)
  loading_factor_labels(r) # validate early
  sim <- simulate_loading_dgp(
    dgp_params, r = r, rhoF_factors = rhoF_factors,
    target_weights = target_weights)
  fit <- estimate_loading_structure(
    sim, Lproxy = Lproxy, projected_updates = projected_updates)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- paste0(r[1L], "x", r[2L])
  # Both raster (slides) and vector (LaTeX) versions use the same grobs.
  pair_base <- file.path(out_dir, paste0("01_matrix_vs_vector_", suffix))
  truth_base <- file.path(out_dir, paste0("02_truth_matrix_vector_", suffix))
  plot_pair <- plot_loading_structure(
    fit, paste0(pair_base, ".png"))
  plot_loading_structure(fit, paste0(pair_base, ".pdf"))
  plot_truth <- plot_loading_structure(
    fit, paste0(truth_base, ".png"),
    show_truth = TRUE)
  plot_loading_structure(fit, paste0(truth_base, ".pdf"),
                         show_truth = TRUE)
  utils::write.csv(fit$metrics,
                   file.path(out_dir, "03_subspace_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(fit$angles,
                   file.path(out_dir, "04_principal_angles.csv"),
                   row.names = FALSE)
  utils::write.csv(fit$separability,
                   file.path(out_dir, "05_rank_one_shares.csv"),
                   row.names = FALSE)
  saveRDS(list(sim = sim, fit = fit),
          file.path(out_dir, "06_sim_and_fit.rds"))
  print(fit$metrics)
  print(fit$separability)
  invisible(list(sim = sim, fit = fit,
                 plot_pair = plot_pair, plot_truth = plot_truth))
}
