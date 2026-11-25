# ==============================================================================
# Appendix | Matrix static factor: output functions only
# Requires: ggplot2
# Input: full-sample scalar-proxy Matrix MF-TPRF fit RDS (already estimated)
# Output: factor time series, column loadings, R C' loading map,
#         and a compact country loading tabular for use inside a minipage.
# ==============================================================================

matrix_page_save_plot <- function(plot, path_stem, width, height) {
  dir.create(dirname(path_stem), recursive = TRUE, showWarnings = FALSE)
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(
      filename = paste0(path_stem, ".", ext), plot = plot,
      width = width, height = height, units = "in", dpi = 500,
      bg = "white", limitsize = FALSE
    )
  }
}

matrix_page_loading_labels <- function(fit, R, C) {
  # Tensor dimensions in this project: time x country x indicator.
  names_x <- dimnames(fit$preprocessing$X_cl)
  if (is.null(rownames(R)) && length(names_x) >= 3L)
    rownames(R) <- names_x[[2L]]
  if (is.null(rownames(C)) && length(names_x) >= 3L)
    rownames(C) <- names_x[[3L]]
  list(R = R, C = C)
}

make_matrix_static_factor_page <- function(
    fit,
    figure_dir = "results/paper/appendix/figures",
    table_dir  = "results/paper/appendix/tables") {

  if (!identical(as.character(fit$proxy_mode), "scalar"))
    stop("Expected the baseline scalar-proxy Matrix fit.")
  if (is.null(fit$fit$R) || is.null(fit$fit$C))
    stop("Saved Matrix fit has no R and C loading matrices.")

  R <- as.matrix(fit$fit$R)
  C <- as.matrix(fit$fit$C)
  if (!identical(dim(R), c(8L, 1L)) || !identical(dim(C), c(50L, 1L)))
    stop("Expected row/column loading dimensions 8x1 and 50x1.")
  if (any(!is.finite(R)) || any(!is.finite(C)))
    stop("Non-finite row or column loadings.")

  named <- matrix_page_loading_labels(fit, R, C)
  R <- named$R
  C <- named$C
  countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  if (is.null(rownames(R)) || !setequal(rownames(R), countries))
    stop("Country names are absent/inconsistent in the saved fit.")
  if (is.null(rownames(C)) || anyNA(rownames(C)) || anyDuplicated(rownames(C)))
    stop("Indicator names are absent/inconsistent in the saved fit.")

  # Eigenvector signs are not identified: change R and C together, if needed.
  # The implied loading map R %*% t(C) does not change.
  if (mean(R[, 1L]) < 0) {
    R <- -R
    C <- -C
  }
  R <- R[countries, , drop = FALSE]

  row_values <- data.frame(
    Country = countries, Loading = as.numeric(R[, 1L]),
    stringsAsFactors = FALSE
  )
  col_values <- data.frame(
    Indicator = rownames(C), Loading = as.numeric(C[, 1L]),
    stringsAsFactors = FALSE
  )
  # Signed descending order; same order in column plot and heatmap.
  row_order <- row_values$Country[order(-row_values$Loading, row_values$Country)]
  col_values <- col_values[order(-col_values$Loading, col_values$Indicator), , drop = FALSE]
  col_order <- col_values$Indicator

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------------------------------------------------------
  # Figure 1: original black-line monthly static-factor time series.
  # Keep the original ggplot title, since the manuscript trims its top margin.
  # ---------------------------------------------------------------------------
  factors <- fit$fit$factors
  F <- factors$F_hf
  dates_m <- as.Date(factors$dates_m)
  if (length(dim(F)) != 3L || dim(F)[2L] != 1L || dim(F)[3L] != 1L ||
      dim(F)[1L] != length(dates_m))
    stop("Expected an aligned monthly factor array T x 1 x 1.")
  factor_df <- data.frame(date = dates_m, value = as.numeric(F[, 1L, 1L]),
                          row_factor = "Row 1", col_factor = "Col 1")
  if (anyNA(factor_df$date)) stop("Missing factor dates.")
  covid_start <- as.Date(fit$params$covid_start)
  covid_end <- as.Date(fit$params$covid_end)
  if (length(covid_start) != 1L || length(covid_end) != 1L ||
      is.na(covid_start) || is.na(covid_end))
    stop("Covid shading dates are missing from fit$params.")

  p_factor <- ggplot2::ggplot(factor_df, ggplot2::aes(x = date, y = value)) +
    ggplot2::annotate("rect", xmin = covid_start, xmax = covid_end,
                      ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.16) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_line(linewidth = 0.8, colour = "black") +
    ggplot2::facet_grid(row_factor ~ col_factor, scales = "free_y") +
    ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(title = "Estimated latent factor", x = "Date", y = "Factor value") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  matrix_page_save_plot(
    p_factor, file.path(figure_dir, "figure_04_baseline_matrix_factor"),
    width = 8.0, height = 4.5
  )

  # ---------------------------------------------------------------------------
  # Figure 2: 50 column loadings, ranked and split into 25 + 25 columns.
  # No subplot titles, as requested; both panels share the loading x-scale.
  # ---------------------------------------------------------------------------
  col_values$Panel <- factor(rep(c("A", "B"), each = 25L), levels = c("A", "B"))
  col_values$Indicator <- factor(col_values$Indicator, levels = rev(col_order))
  xmin <- min(0, col_values$Loading)
  xmax <- max(0, col_values$Loading)
  pad <- max(0.04 * (xmax - xmin), 0.02)

  p_col <- ggplot2::ggplot(
    col_values,
    ggplot2::aes(x = Loading, y = Indicator, fill = Loading)
  ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.32, colour = "grey50") +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~ Panel, nrow = 1, scales = "free_y") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, guide = "none"
    ) +
    ggplot2::coord_cartesian(xlim = c(xmin - pad, xmax + pad)) +
    ggplot2::labs(x = "Loading", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey91", linewidth = 0.25),
      panel.spacing.x = grid::unit(0.18, "in"),
      axis.text.y = ggplot2::element_text(size = 9, colour = "grey20"),
      axis.text.x = ggplot2::element_text(size = 9, colour = "grey25"),
      legend.position = "none",
      plot.margin = ggplot2::margin(2, 3, 2, 2)
    )
  matrix_page_save_plot(
    p_col, file.path(figure_dir, "column_loadings"),
    width = 7.2, height = 4.1
  )

  # ---------------------------------------------------------------------------
  # Figure 3: the country x indicator version of C %x% R.
  # Matrix entries are R_i * C_j, with 8 country rows and 50 indicator columns.
  # ---------------------------------------------------------------------------
  map_matrix <- R[row_order, 1L, drop = FALSE] %*%
    t(C[col_order, 1L, drop = FALSE])
  dimnames(map_matrix) <- list(row_order, col_order)
  map_df <- as.data.frame(as.table(map_matrix), stringsAsFactors = FALSE)
  names(map_df) <- c("Country", "Indicator", "Loading")
  map_df$Country <- factor(map_df$Country, levels = rev(row_order))
  map_df$Indicator <- factor(map_df$Indicator, levels = col_order)
  color_limit <- max(abs(map_df$Loading))
  if (color_limit == 0) color_limit <- 1

  p_map <- ggplot2::ggplot(
    map_df, ggplot2::aes(x = Indicator, y = Country, fill = Loading)
  ) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-color_limit, color_limit), name = "Loading"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5, size = 8.2, colour = "grey20"
      ),
      axis.text.y = ggplot2::element_text(size = 10.5, colour = "grey20"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  matrix_page_save_plot(
    p_map, file.path(figure_dir, "loading_map"),
    width = 10.2, height = 4.1
  )

  # ---------------------------------------------------------------------------
  # Country loading table: ONLY the tabular environment. The manuscript's
  # captionof{table} and labels stay outside, in the minipage.
  # ---------------------------------------------------------------------------
  tex_eol <- paste0(" ", strrep("\\", 2L))
  tab_lines <- c(
    "\\begin{tabular}{lr}",
    "\\toprule",
    paste0("Country & Loading", tex_eol),
    "\\midrule",
    sprintf("%s & %.3f%s", row_values$Country, row_values$Loading, tex_eol),
    "\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(tab_lines, file.path(table_dir, "row_loadings_vertical.tex"), useBytes = TRUE)
  utils::write.csv(row_values, file.path(table_dir, "row_loadings.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(Indicator = col_order,
               Loading = as.numeric(C[col_order, 1L])),
    file.path(table_dir, "column_loadings.csv"), row.names = FALSE
  )

  message("Matrix static-factor page exported: ", figure_dir)
  invisible(list(
    factor = p_factor, column = p_col, loading_map = p_map,
    rows = row_values, columns = col_values,
    matrix = map_matrix
  ))
}
