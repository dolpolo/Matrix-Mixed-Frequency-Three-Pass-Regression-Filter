# ==============================================================================
# Appendix | Matrix versus Vector factor comparison: output functions only
# Requires: ggplot2
# Inputs: saved scalar-proxy Matrix fit and Euro Area vector (VEC-C) fit RDSs
# Outputs: EA GDP YoY factor figure, 2x2 R2/adj-R2 LaTeX table, diagnostics CSV.
# ==============================================================================

ea_page_save_plot <- function(plot, path_stem, width, height) {
  dir.create(dirname(path_stem), recursive = TRUE, showWarnings = FALSE)
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(
      filename = paste0(path_stem, ".", ext), plot = plot,
      width = width, height = height, units = "in", dpi = 500,
      bg = "white", limitsize = FALSE
    )
  }
}

ea_page_factor_table <- function(fit, vector = FALSE) {
  if (vector) {
    factors <- fit$fit
    quarterly_dates <- as.Date(fit$dates_q)
  } else {
    factors <- fit$fit$factors
    quarterly_dates <- as.Date(factors$dates_q)
  }
  get_month <- function(key) {
    x <- as.matrix(factors[[key]])
    if (is.null(x) || ncol(x) < 1L || nrow(x) != length(quarterly_dates))
      stop("Invalid quarterly factor: ", key)
    as.numeric(x[, 1L])
  }
  if (anyNA(quarterly_dates) || anyDuplicated(quarterly_dates))
    stop("Missing or duplicated quarterly factor dates.")
  suffix <- if (vector) "Vector" else "Matrix"
  d <- data.frame(date = quarterly_dates,
                  m1 = get_month("F1"),
                  m2 = get_month("F2"),
                  m3 = get_month("F3"))
  names(d)[-1L] <- paste0(suffix, c("_M1", "_M2", "_M3"))
  d
}

ea_page_yoy_from_qoq <- function(matrix_fit) {
  # Use the identical four-quarter sum convention as the earlier EA figure.
  dates_q <- as.Date(matrix_fit$dates_q)
  y <- as.matrix(matrix_fit$target$y_q_all)
  if (!"EA" %in% colnames(y) || nrow(y) != length(dates_q))
    stop("EA GDP quarterly series is absent or misaligned.")
  if (length(dates_q) < 5L || anyNA(dates_q) || anyDuplicated(dates_q))
    stop("Quarterly target dates are missing, duplicated, or too short.")
  qoq <- as.numeric(y[, "EA"])
  lag_values <- function(x, k) c(rep(NA_real_, k), head(x, -k))
  yoy <- qoq + lag_values(qoq, 1L) + lag_values(qoq, 2L) +
    lag_values(qoq, 3L)
  data.frame(date = dates_q, GDP_YoY = yoy)
}

ea_page_r2_fit <- function(df, suffix, model_name) {
  avg <- rowMeans(df[, paste0(suffix, c("_M1", "_M2", "_M3"))])
  dat <- data.frame(
    GDP_YoY = df$GDP_YoY, F_quarter = avg,
    M1 = df[[paste0(suffix, "_M1")]],
    M2 = df[[paste0(suffix, "_M2")]],
    M3 = df[[paste0(suffix, "_M3")]]
  )
  fit_avg <- stats::lm(GDP_YoY ~ F_quarter, data = dat)
  fit_month <- stats::lm(GDP_YoY ~ M1 + M2 + M3, data = dat)
  r_avg <- summary(fit_avg)
  r_month <- summary(fit_month)
  if (r_month$r.squared + 1e-9 < r_avg$r.squared)
    stop("Unexpected R2 ordering: check that the samples are identical.")
  data.frame(
    Model = model_name,
    Specification = c("Quarterly average", "M1 + M2 + M3"),
    N = c(stats::nobs(fit_avg), stats::nobs(fit_month)),
    R2 = c(r_avg$r.squared, r_month$r.squared),
    Adj_R2 = c(r_avg$adj.r.squared, r_month$adj.r.squared),
    stringsAsFactors = FALSE
  )
}

ea_page_write_latex <- function(results, file) {
  stopifnot(nrow(results) == 4L)
  get_row <- function(model) {
    q <- results[results$Model == model &
                   results$Specification == "Quarterly average", ]
    m <- results[results$Model == model &
                   results$Specification == "M1 + M2 + M3", ]
    if (nrow(q) != 1L || nrow(m) != 1L)
      stop("Missing/duplicated factor regression for model: ", model)
    sprintf("%.1f", 100 * c(q$R2, q$Adj_R2, m$R2, m$Adj_R2))
  }
  eol <- paste0(" ", strrep("\\", 2L))
  model_row <- function(name, label) {
    paste0(label, " & ", paste(get_row(name), collapse = " & "), eol)
  }
  tex <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\renewcommand{\\arraystretch}{1.15}",
    "\\begin{tabular}{lcccc}",
    "\\toprule",
    paste0(" & \\multicolumn{2}{c}{Quarterly average} & ",
           "\\multicolumn{2}{c}{Monthly factors ($M1$--$M3$)}", eol),
    "\\cmidrule(lr){2-3}",
    "\\cmidrule(lr){4-5}",
    paste0("Model & $R^2$ & Adj. $R^2$ & $R^2$ & Adj. $R^2$", eol),
    "\\midrule",
    model_row("Matrix MF-TPRF", "Matrix MF--TPRF"),
    model_row("Vector MF-TPRF", "Vector MF--TPRF"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\caption{In-sample explanatory power for euro-area GDP growth (YoY).}",
    "\\label{tab:ea_factor_fit}",
    "\\vspace{0.1cm}",
    "\\begin{minipage}{0.85\\textwidth}",
    "\\footnotesize",
    "\\textit{Notes:} The table reports $R^2$ and adjusted $R^2$ in percentages. The quarterly-average specification uses the mean of the three monthly factors, while the alternative includes $M1$, $M2$, and $M3$ separately.",
    "\\end{minipage}",
    "\\end{table}"
  )
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(tex, file, useBytes = TRUE)
}

make_ea_factor_comparison_page <- function(
    matrix_fit, vector_fit,
    figure_dir = "results/paper/appendix/figures",
    table_dir  = "results/paper/appendix/tables") {

  if (!identical(as.character(matrix_fit$proxy_mode), "scalar"))
    stop("Expected scalar-proxy Matrix MF-TPRF baseline fit.")

  gdp <- ea_page_yoy_from_qoq(matrix_fit)
  M <- ea_page_factor_table(matrix_fit, vector = FALSE)
  V <- ea_page_factor_table(vector_fit, vector = TRUE)
  if (anyDuplicated(gdp$date)) stop("Duplicate target quarters.")

  # All four regressions use EXACTLY the same quarterly observations.
  common <- merge(merge(gdp, M, by = "date"), V, by = "date")
  numeric_cols <- setdiff(names(common), "date")
  common <- common[stats::complete.cases(common[, numeric_cols, drop = FALSE]), ]
  common <- common[order(common$date), , drop = FALSE]
  if (nrow(common) < 8L) stop("Too few common quarters for factor comparison.")

  matrix_fit_stats <- ea_page_r2_fit(common, "Matrix", "Matrix MF-TPRF")
  vector_fit_stats <- ea_page_r2_fit(common, "Vector", "Vector MF-TPRF")
  fit_results <- rbind(matrix_fit_stats, vector_fit_stats)
  if (length(unique(fit_results$N)) != 1L)
    stop("Regressions are not using the same quarterly sample.")

  # Orient factor signs for the graph only. OLS fits are unchanged by sign.
  standardize <- function(x) {
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) stop("Cannot standardize a constant series.")
    as.numeric(scale(x))
  }
  quarter_mean <- function(prefix) rowMeans(common[, paste0(prefix, c("_M1", "_M2", "_M3"))])
  matrix_f <- quarter_mean("Matrix")
  vector_f <- quarter_mean("Vector")
  sign_m <- sign(stats::cor(common$GDP_YoY, matrix_f))
  sign_v <- sign(stats::cor(common$GDP_YoY, vector_f))
  if (!is.finite(sign_m) || !is.finite(sign_v) || sign_m == 0 || sign_v == 0)
    stop("Factor cannot be oriented relative to EA YoY GDP.")

  chart_data <- rbind(
    data.frame(date = common$date, series = "EA GDP YoY", value = standardize(common$GDP_YoY)),
    data.frame(date = common$date, series = "Matrix MF-TPRF", value = standardize(sign_m * matrix_f)),
    data.frame(date = common$date, series = "Vector MF-TPRF", value = standardize(sign_v * vector_f))
  )
  chart_data$series <- factor(
    chart_data$series,
    levels = c("EA GDP YoY", "Matrix MF-TPRF", "Vector MF-TPRF")
  )

  # Original paper color and linetype scheme. Preserve the title because the
  # author's LaTeX includegraphics instruction trims the upper figure margin.
  p <- ggplot2::ggplot(chart_data,
                       ggplot2::aes(x = date, y = value,
                                    colour = series, linetype = series)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::scale_colour_manual(
      values = c("EA GDP YoY" = "#C00000", "Matrix MF-TPRF" = "#228B22",
                 "Vector MF-TPRF" = "#1F77B4"), name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("EA GDP YoY" = "solid", "Matrix MF-TPRF" = "solid",
                 "Vector MF-TPRF" = "33"), name = NULL
    ) +
    ggplot2::scale_x_date(
      date_breaks = "3 years", date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(title = "EA year-on-year GDP growth and estimated factors",
                  x = NULL, y = "Standardised units") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  ea_page_save_plot(p, file.path(figure_dir, "figure_06_ea_yoy_factor_comparison"),
                    width = 11.5, height = 5.0)
  ea_page_write_latex(fit_results, file.path(table_dir, "ea_factor_fit_yoy.tex"))
  utils::write.csv(fit_results, file.path(table_dir, "ea_factor_fit_yoy.csv"), row.names = FALSE)
  message("EA factor-comparison page exported: ", figure_dir)
  print(fit_results, row.names = FALSE)
  invisible(list(plot = p, regressions = fit_results, sample = common))
}
