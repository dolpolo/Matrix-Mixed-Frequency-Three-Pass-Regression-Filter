# Author: Davide Delfino
# Input: see function arguments; source from a model main.
# Output: function return values; no files written unless documented.
# Source: original functions_mat/matrix.mf.tprf.wald.R; numerical implementation retained.

# ==============================================================================
# MATRIX MF-TPRF: WALD TESTS FOR THE EIGHT ROW LOADINGS
# ==============================================================================
# Tests H0: R_{0,i} = 1, separately for DE, FR, IT, ES, NL, BE, AT, PT.
#
# Uses the FIRST projected update R^(1), as in Pass 1, Theorem 3.
# Estimates the score covariance by Newey-West (Bartlett weights).
#
# IMPORTANT: exploratory inference. Consistency of the projected-residual
# score and validity of the chi-square approximation for p1 = 8 fixed
# remain to be established.
# ==============================================================================

matrix_mf_tprf_hac <- function(scores, lags) {
  # scores: quarters x countries
  scores <- scale(scores, center = TRUE, scale = FALSE)
  Tq <- nrow(scores)
  p1 <- ncol(scores)

  V <- crossprod(scores) / Tq

  if (lags > 0) {
    for (h in seq_len(lags)) {
      cov_h <- crossprod(
        scores[(h + 1):Tq, , drop = FALSE],
        scores[1:(Tq - h), , drop = FALSE]
      ) / Tq

      weight <- 1 - h / (lags + 1)
      V <- V + weight * (cov_h + t(cov_h))
    }
  }

  dimnames(V) <- list(colnames(scores), colnames(scores))
  V
}

matrix_mf_tprf_row_wald <- function(res_fit, hac_lags = 4L) {

  # 1. DATA AND CHECKS -------------------------------------------------------
  countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  X <- res_fit$preprocessing$X_cl_q
  fit <- res_fit$fit

  if (length(dim(X)) != 3L) stop("Missing quarterly tensor X_cl_q.")
  if (is.null(fit$Z_q)) stop("Missing proxy matrix fit$Z_q.")
  if (!identical(as.integer(fit$r_selected), c(1L, 1L))) {
    stop("The Wald function requires r1 = r2 = 1.")
  }
  if (!setequal(dimnames(X)[[2]], countries)) {
    stop("Quarterly tensor must contain exactly the eight countries.")
  }

  # The multivariate proxy uses a shorter, aligned quarterly window.
  rows <- fit$target_proxy_rows
  if (is.null(rows)) rows <- seq_len(dim(X)[1])
  if (anyNA(rows) || any(rows < 1L | rows > dim(X)[1])) {
    stop("Invalid target_proxy_rows.")
  }
  X <- X[rows, countries, , drop = FALSE]
  Z <- as.matrix(fit$Z_q)

  Tq <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(X)[3]
  L  <- ncol(Z)

  if (nrow(Z) != Tq) stop("Quarterly tensor and proxies are not aligned.")
  if (Tq < 8L) stop("Too few quarters for HAC inference.")
  if (any(!is.finite(X)) || any(!is.finite(Z))) {
    stop("X or Z contains NA, NaN or Inf.")
  }
  if (length(hac_lags) != 1L || !is.finite(hac_lags) ||
      hac_lags < 0 || hac_lags != as.integer(hac_lags) || hac_lags >= Tq) {
    stop("hac_lags must be an integer between 0 and Tq - 1.")
  }
  hac_lags <- as.integer(hac_lags)

  # Scalar proxies are QR-orthogonalized in the original Pass 1.
  # Multiplying by sqrt(Tq) restores a non-vanishing sample second moment.
  # This common rescaling does not change the selected eigenvectors.
  if (identical(as.character(res_fit$proxy_mode), "scalar")) {
    Z <- scale(Z, center = TRUE, scale = FALSE)
    if (qr(Z)$rank < L) stop("The proxy columns are linearly dependent.")
    Z <- sqrt(Tq) * qr.Q(qr(Z))
  }

  # 2. PROXY-WEIGHTED MOMENTS M_l -------------------------------------------
  M <- lapply(seq_len(L), function(l) matrix(0, p1, p2))

  for (tau in seq_len(Tq)) {
    X_tau <- X[tau, , ]
    for (l in seq_len(L)) {
      M[[l]] <- M[[l]] + Z[tau, l] * X_tau / Tq
    }
  }

  # 3. FIRST PROJECTED UPDATE: C^(0), R^(1), C^(1) ----------------------------
  S_C0 <- matrix(0, p2, p2)
  for (l in seq_len(L)) {
    S_C0 <- S_C0 + crossprod(M[[l]])
  }
  C0 <- sqrt(p2) * eigen(S_C0, symmetric = TRUE)$vectors[, 1]

  S_R1 <- matrix(0, p1, p1)
  for (l in seq_len(L)) {
    S_R1 <- S_R1 +
      M[[l]] %*% tcrossprod(C0) %*% t(M[[l]]) / p2
  }
  eig_R <- eigen(S_R1, symmetric = TRUE)
  R1 <- sqrt(p1) * eig_R$vectors[, 1]

  # The test of R_i = +1 needs a fixed sign convention.
  if (R1[1] < 0) R1 <- -R1
  names(R1) <- countries

  S_C1 <- matrix(0, p2, p2)
  for (l in seq_len(L)) {
    S_C1 <- S_C1 +
      t(M[[l]]) %*% tcrossprod(R1) %*% M[[l]] / p1
  }
  C1 <- sqrt(p2) * eigen(S_C1, symmetric = TRUE)$vectors[, 1]

  # Theorem 3 uses the eigenvalue divided by p1*p2.
  lambda_R <- eig_R$values[1] / (p1 * p2)
  if (!is.finite(lambda_R) || lambda_R <= 0) {
    stop("The leading normalized row eigenvalue is not positive.")
  }

  # 4. QUARTERLY FACTORS AND PROJECTED RESIDUALS -----------------------------
  F_q <- numeric(Tq)
  projected_residuals <- matrix(NA_real_, Tq, p1)

  for (tau in seq_len(Tq)) {
    X_tau <- X[tau, , ]
    F_q[tau] <- as.numeric(crossprod(R1, X_tau %*% C1)) / (p1 * p2)

    # X_tau - R^(1) F_tau C^(1)' ; project along C^(0).
    X_res <- X_tau - F_q[tau] * tcrossprod(R1, C1)
    projected_residuals[tau, ] <- as.vector(X_res %*% C0) / sqrt(p2)
  }

  # 5. SCORE IN ASSUMPTION F AND HAC VARIANCE --------------------------------
  # Sample counterpart of E[z^(l) F_0,tau].
  EzF <- as.vector(crossprod(Z, F_q) / Tq)

  # A row of scores contains the eight country-specific contributions.
  proxy_weight <- as.vector(Z %*% EzF)
  scores <- sweep(projected_residuals, 1, proxy_weight, "*")
  colnames(scores) <- countries

  V_R <- matrix_mf_tprf_hac(scores, lags = hac_lags)
  omega_R <- diag(V_R) / lambda_R^2
  var_R <- omega_R / (Tq * p2)

  if (any(!is.finite(var_R)) || any(var_R <= 0)) {
    stop("Non-positive score variance. Inspect the projected residuals.")
  }

  # 6. EIGHT INDIVIDUAL WALD TESTS AND BONFERRONI ----------------------------
  wald <- (R1 - 1)^2 / var_R
  p_value <- pchisq(wald, df = 1, lower.tail = FALSE)
  p_bonf <- p.adjust(p_value, method = "bonferroni")

  results <- data.frame(
    country = countries,
    loading = as.numeric(R1),
    se = sqrt(var_R),
    wald = as.numeric(wald),
    p_value = as.numeric(p_value),
    p_bonferroni = as.numeric(p_bonf),
    reject_5pct = as.logical(p_bonf < 0.05),
    row.names = NULL
  )

  list(
    results = results,
    Tq = Tq, p1 = p1, p2 = p2, L = L,
    hac_lags = hac_lags,
    lambda_R = lambda_R,
    R_first = R1,
    R_ils = fit$R,
    C_initial = C0,
    C_first = C1,
    V_R = V_R,
    omega_R = omega_R,
    scores = scores,
    note = "Exploratory: fixed-p1 and feasible-score inference not yet justified."
  )
}


# ==============================================================================
# MATRIX MF-TPRF: LATEX TABLE FOR EIGHT ROW-LOADING WALD TESTS
# ============================================================================== 
# Requires: a result from matrix_mf_tprf_row_wald().
# Requires: \usepackage{booktabs} in the LaTeX preamble.

matrix_mf_tprf_wald_latex <- function(wald_result, file) {
  tab <- wald_result$results
  
  # Display very small p-values without printing misleading zeros.
  show_p <- function(p) {
    if (p < 0.0001) "$<10^{-4}$" else sprintf("%.4f", p)
  }
  
  # Two LaTeX backslashes terminate a table row.
  row_end <- " \\\\" 
  
  latex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Wald tests for predictive row loadings}",
    "\\label{tab:matrix_mf_tprf_row_wald}",
    "\\small",
    "\\begin{tabular}{lrrrrr}",
    "\\toprule",
    paste0(
      "Country & $\\widehat R_i^{(1)}$ & Std. error & Wald ",
      "& $p$-value & Bonferroni $p$-value", row_end
    ),
    "\\midrule"
  )
  
  for (i in seq_len(nrow(tab))) {
    line <- paste(
      tab$country[i],
      sprintf("%.3f", tab$loading[i]),
      sprintf("%.3f", tab$se[i]),
      sprintf("%.2f", tab$wald[i]),
      show_p(tab$p_value[i]),
      show_p(tab$p_bonferroni[i]),
      sep = " & "
    )
    latex <- c(latex, paste0(line, row_end))
  }
  
  latex <- c(
    latex,
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.12cm}",
    "\\begin{minipage}{0.96\\textwidth}",
    "\\footnotesize \\emph{Notes:} Each country is tested separately ",
    "against $H_{0,i}:R_{0,i}=1$. Standard errors use a Bartlett HAC ",
    paste0("estimator with ", wald_result$hac_lags, " quarterly lags. "),
    "Bonferroni-adjusted $p$-values control eight comparisons. ",
    "The inference is exploratory pending a fixed-$p_1$ justification.",
    "\\end{minipage}",
    "\\end{table}"
  )
  
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(latex, con = file, useBytes = TRUE)
  invisible(latex)
}

# ==============================================================================
# MATRIX MF-TPRF: JOINT WALD TEST OF EQUAL COUNTRY ROW LOADINGS
# ==============================================================================
# H0: R_DE = R_FR = R_IT = R_ES = R_NL = R_BE = R_AT = R_PT
#
# Uses the first projected row loadings and the FULL HAC covariance already
# computed by matrix_mf_tprf_row_wald(). No new model estimation is needed.
#
# IMPORTANT: Exploratory inference. Fixed-p1 validity and consistency of the
# estimated projected-residual score still need a theoretical justification.
# ==============================================================================

matrix_mf_tprf_joint_wald <- function(wald_result) {
  
  # 1. Inputs ---------------------------------------------------------------
  countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  
  R <- wald_result$R_first[countries]
  V_R <- wald_result$V_R[countries, countries, drop = FALSE]
  lambda_R <- wald_result$lambda_R
  Tq <- wald_result$Tq
  p2 <- wald_result$p2
  
  if (length(R) != 8L || any(!is.finite(R))) {
    stop("Missing row loadings. Run matrix_mf_tprf_row_wald() first.")
  }
  if (!identical(dim(V_R), c(8L, 8L)) || any(!is.finite(V_R))) {
    stop("Missing 8 x 8 HAC covariance matrix V_R.")
  }
  if (length(lambda_R) != 1L || !is.finite(lambda_R) || lambda_R <= 0) {
    stop("The normalized leading eigenvalue must be positive.")
  }
  if (!is.finite(Tq) || !is.finite(p2) || Tq <= 0 || p2 <= 0) {
    stop("Invalid sample dimensions Tq or p2.")
  }
  
  # 2. Seven restrictions: each country equals Germany ---------------------
  # Each row computes R_country - R_DE.
  Q <- cbind(-rep(1, 7), diag(7))
  rownames(Q) <- paste0(countries[-1], " - DE")
  colnames(Q) <- countries
  
  differences <- as.vector(Q %*% R)
  
  # 3. Covariance matrix of the seven estimated differences ----------------
  # V_R is the HAC covariance of the score (Assumption F).
  # Dividing by lambda_R^2 and Tq*p2 yields the estimated covariance
  # of the row loadings (under the proposed plug-in approximation).
  covariance_R <- V_R / (lambda_R^2 * Tq * p2)
  covariance_diff <- Q %*% covariance_R %*% t(Q)
  covariance_diff <- (covariance_diff + t(covariance_diff)) / 2
  
  # The seven restriction directions must have a nonsingular covariance.
  eigenvalues <- eigen(covariance_diff, symmetric = TRUE,
                       only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) ||
      min(eigenvalues) <= max(eigenvalues) * 1e-9) {
    stop("Covariance of the seven restrictions is singular or nearly singular.")
  }
  
  # 4. Joint Wald statistic ------------------------------------------------
  wald <- as.numeric(t(differences) %*%
                       solve(covariance_diff, differences))
  df <- 7L
  p_value <- pchisq(wald, df = df, lower.tail = FALSE)
  
  results <- data.frame(
    test = "Equal row loadings across eight countries",
    wald = wald,
    df = df,
    p_value = p_value,
    reject_5pct = p_value < 0.05
  )
  
  contrasts <- data.frame(
    comparison = paste0(countries[-1], " - DE"),
    difference = differences,
    se = sqrt(diag(covariance_diff))
  )
  
  list(
    results = results,
    contrasts = contrasts,
    restrictions = Q,
    covariance_diff = covariance_diff,
    hac_lags = wald_result$hac_lags,
    note = "Exploratory: validity for fixed p1=8 remains to be established."
  )
}


# ==============================================================================
# LATEX TABLE FOR THE JOINT WALD TEST
# ==============================================================================
# Requires \usepackage{booktabs} in the LaTeX preamble.

matrix_mf_tprf_joint_wald_latex <- function(joint_result, file) {
  
  tab <- joint_result$contrasts
  test <- joint_result$results
  row_end <- " \\\\"
  
  # Print small p-values in scientific notation rather than as 0.0000.
  show_p <- function(p) {
    if (p == 0) return("$<10^{-300}$")
    if (p < 0.0001) {
      exponent <- floor(log10(p))
      coefficient <- p / 10^exponent
      sprintf("$%.2f\\times 10^{%d}$", coefficient, exponent)
    } else {
      sprintf("%.4f", p)
    }
  }
  
  latex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Joint Wald test for equal predictive row loadings}",
    "\\label{tab:matrix_mf_tprf_joint_wald}",
    "\\small",
    "\\begin{tabular}{lrr}",
    "\\toprule",
    paste0("Restriction & Estimated difference & Std. error", row_end),
    "\\midrule"
  )
  
  for (i in seq_len(nrow(tab))) {
    latex <- c(latex, paste0(
      tab$comparison[i], " & ",
      sprintf("%.3f", tab$difference[i]), " & ",
      sprintf("%.3f", tab$se[i]), row_end
    ))
  }
  
  latex <- c(
    latex,
    "\\midrule",
    paste0("Joint Wald statistic & \\multicolumn{2}{r}{",
           sprintf("%.3f", test$wald), "}", row_end),
    paste0("Degrees of freedom & \\multicolumn{2}{r}{",
           test$df, "}", row_end),
    paste0("$p$-value & \\multicolumn{2}{r}{",
           show_p(test$p_value), "}", row_end),
    "\\bottomrule",
    "\\end{tabular}",
    "\\vspace{0.12cm}",
    "\\begin{minipage}{0.96\\textwidth}",
    paste0("\\footnotesize \\emph{Notes:} The null is equality of all eight ",
           "predictive row loadings. Germany is the reference country. ",
           "The joint test imposes seven restrictions and uses a Bartlett ",
           "HAC covariance matrix with ", joint_result$hac_lags,
           " quarterly lags. Inference is exploratory pending a ",
           "fixed-$p_1$ justification."),
    "\\end{minipage}",
    "\\end{table}"
  )
  
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(latex, file, useBytes = TRUE)
  invisible(latex)
}


