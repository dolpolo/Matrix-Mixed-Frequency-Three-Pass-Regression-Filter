# ==============================================================================
# Euro Area GDP growth nowcasts | Matrix MF-TPRF versus VEC-C
# Functions only: read saved results, compute matched RMSFEs, make original plot
# Input:  results/empirical/matrix and results/empirical/vector/EA (fixed50)
# No model estimation or predictor selection is performed here.
# ==============================================================================

# Requires ggplot2, dplyr, tidyr, scales and the project's
# R/matrix/matrix.mf.tprf.utils.R, sourced by the main script.

load_ea_nowcast_results <- function(size = "fixed50", selection = "fixed") {
  required_helpers <- c("find_result_file", "extract_matrix_ea_nowcasts",
                        "extract_vector_ea_nowcasts", "extract_matrix_ea_gdp",
                        "extract_vector_ea_gdp", "make_ea_quarter_id",
                        "theme_country_compare", "get_ea_gdp_ylim")
  absent <- required_helpers[!vapply(required_helpers, exists, logical(1),
                                     mode = "function")]
  if (length(absent)) {
    stop("First source R/matrix/matrix.mf.tprf.utils.R. Missing: ",
         paste(absent, collapse = ", "))
  }

  file_matrix <- find_result_file(
    path = "results/empirical/matrix", model = "matrix_scalar", stage = "rt",
    Size = size, sel = selection
  )
  file_vecc <- find_result_file(
    path = "results/empirical/vector/EA", model = "vector", stage = "rt",
    Size = size, sel = selection
  )
  matrix_rt <- readRDS(file_matrix)
  vecc_rt <- readRDS(file_vecc)

  if (!identical(as.character(matrix_rt$proxy_mode), "scalar")) {
    stop("The Matrix EA result is not the baseline scalar/aggregate-target proxy design.")
  }
  if (!identical(as.character(matrix_rt$Size), size) ||
      !identical(as.character(vecc_rt$Size), size)) {
    stop("A loaded EA RDS does not match the expected information-set size.")
  }
  list(matrix = matrix_rt, vecc = vecc_rt,
       files = c(matrix = file_matrix, vecc = file_vecc))
}

# Build a strictly paired evaluation sample for the RMSFE table.
# Each row of 'matched' is one common target-quarter vintage; GDP is attached
# by reference quarter, not by the calendar month when the nowcast was issued.
prepare_ea_nowcast_data <- function(saved, cfg) {
  months <- c("M1", "M2", "M3")
  start_eval <- as.Date(cfg$start_eval)
  end_eval <- as.Date(cfg$end_eval)
  if (is.na(start_eval) || is.na(end_eval) || end_eval < start_eval) {
    stop("Invalid evaluation dates in replication_config.")
  }

  m <- extract_matrix_ea_nowcasts(saved$matrix, month_order = months)
  v <- extract_vector_ea_nowcasts(saved$vecc, month_order = months)
  m <- dplyr::filter(m, date >= start_eval, date <= end_eval,
                     type %in% months, is.finite(nowcast))
  v <- dplyr::filter(v, date >= start_eval, date <= end_eval,
                     type %in% months, is.finite(nowcast))

  # Normalise to first day of the calendar month before pairing.
  # The underlying nowcast vintage (M1/M2/M3) is not changed.
  m$date <- as.Date(format(as.Date(m$date), "%Y-%m-01"))
  v$date <- as.Date(format(as.Date(v$date), "%Y-%m-01"))
  for (nm in c("Matrix MF-TPRF", "VEC-C")) {
    obj <- if (nm == "Matrix MF-TPRF") m else v
    if (anyDuplicated(obj[c("date", "type")])) {
      stop("Duplicate EA nowcasts after date normalisation: ", nm)
    }
  }

  paired <- dplyr::inner_join(
    dplyr::select(m, date, type, matrix_nowcast = nowcast),
    dplyr::select(v, date, type, vecc_nowcast = nowcast),
    by = c("date", "type")
  )
  if (!nrow(paired)) stop("No common EA nowcasts for Matrix MF-TPRF and VEC-C.")
  paired$quarter_id <- make_ea_quarter_id(paired$date)
  paired <- dplyr::arrange(paired, date, factor(type, levels = months))

  gdp_matrix <- extract_matrix_ea_gdp(saved$matrix)
  gdp_vecc <- extract_vector_ea_gdp(saved$vecc)
  gdp_check <- dplyr::inner_join(
    dplyr::select(gdp_matrix, quarter_id, GDP_matrix = GDP),
    dplyr::select(gdp_vecc, quarter_id, GDP_vecc = GDP),
    by = "quarter_id"
  )
  disagreements <- dplyr::filter(gdp_check, is.finite(GDP_matrix),
                                 is.finite(GDP_vecc),
                                 abs(GDP_matrix - GDP_vecc) > 1e-10)
  if (nrow(disagreements)) {
    warning("Matrix and VEC-C EA GDP values differ in some quarters; ",
            "the Matrix GDP series is used for both benchmarks.")
  }

  gdp_valid <- dplyr::filter(gdp_matrix, is.finite(GDP))
  if (anyDuplicated(gdp_valid$quarter_id)) {
    stop("Repeated quarterly EA GDP targets in the Matrix RT result.")
  }
  gdp_valid <- dplyr::rename(gdp_valid, target_date = date)
  matched <- dplyr::inner_join(paired, gdp_valid, by = "quarter_id")
  if (!nrow(matched)) stop("No common EA nowcasts with evaluable quarterly GDP.")

  # Exactly the same target date, country, month and sample for both models.
  matched$period <- dplyr::case_when(
    matched$target_date < as.Date(cfg$covid_start) ~ "Pre-Covid",
    matched$target_date <= as.Date(cfg$covid_end) ~ "Covid",
    matched$target_date > as.Date(cfg$covid_end) ~ "Post-Covid",
    TRUE ~ NA_character_
  )
  if (anyNA(matched$period)) stop("Unclassified EA target dates.")
  matched$se_matrix <- (matched$matrix_nowcast - matched$GDP)^2
  matched$se_vecc <- (matched$vecc_nowcast - matched$GDP)^2
  matched <- dplyr::arrange(matched, date, factor(type, levels = months))

  # Plot the same common forecast dates as the source figure; when GDP for
  # the final target quarter is not yet available, show its nowcasts anyway.
  pred_long <- tidyr::pivot_longer(
    paired, cols = c("matrix_nowcast", "vecc_nowcast"),
    names_to = "source", values_to = "value"
  )
  pred_long$series <- ifelse(pred_long$source == "matrix_nowcast",
                             "Matrix MF-TPRF", "VEC-C")
  pred_long <- dplyr::select(pred_long, date, type, value, series)
  gdp_long <- dplyr::transmute(matched,
                               date, type, value = GDP,
                               series = "EA GDP")
  plot_data <- dplyr::bind_rows(pred_long, gdp_long)
  plot_data$type <- factor(plot_data$type, levels = months)
  plot_data$series <- factor(
    plot_data$series,
    levels = c("EA GDP", "Matrix MF-TPRF", "VEC-C")
  )
  plot_data <- dplyr::arrange(plot_data, type, series, date)

  # Diagnostics: rows in 'matched' are the common, evaluable vintage targets.
  coverage <- matched |>
    dplyr::count(period, type, name = "N") |>
    dplyr::arrange(factor(period, levels = c("Pre-Covid", "Covid", "Post-Covid")),
                   factor(type, levels = months))
  list(matched = matched, plot = plot_data, gdp = gdp_valid,
       coverage = coverage, paired = paired)
}

compute_ea_rmsfe_table <- function(data) {
  periods <- c("Pre-Covid", "Covid", "Post-Covid")
  months <- c("M1", "M2", "M3")
  x <- data$matched
  rmsfe <- x |>
    dplyr::group_by(period, type) |>
    dplyr::summarise(
      N = dplyr::n(),
      Matrix = sqrt(mean(se_matrix)),
      VecC = sqrt(mean(se_vecc)),
      .groups = "drop"
    ) |>
    dplyr::mutate(Relative = Matrix / VecC)
  if (any(!is.finite(rmsfe$Relative))) {
    stop("Non-finite EA relative RMSFE; check zero/missing benchmark errors.")
  }
  required <- expand.grid(period = periods, type = months,
                          stringsAsFactors = FALSE)
  if (nrow(dplyr::anti_join(required, rmsfe, by = c("period", "type")))) {
    stop("Not every evaluation period has all three EA vintages.")
  }
  rmsfe <- rmsfe |>
    dplyr::mutate(period = factor(period, levels = periods),
                  type = factor(type, levels = months)) |>
    dplyr::arrange(period, type)

  wide <- rmsfe |>
    dplyr::select(period, type, Matrix, VecC, Relative) |>
    tidyr::pivot_wider(names_from = type,
                       values_from = c(Matrix, VecC, Relative),
                       names_glue = "{.value}_{type}")
  list(long = rmsfe, wide = wide)
}

# Original three-block, booktabs-based EA RMSFE table, without obsolete
# LASSO/correlation categories. The nine numbers per period are from RDS files.
make_ea_rmsfe_latex <- function(table,
    caption = "EA GDP growth pseudo-real-time RMSFE: Matrix MF--TPRF versus VEC-C",
    label = "tab:ea_rmsfe_matrix_vector") {
  needed <- c(paste0("Matrix_M", 1:3), paste0("VecC_M", 1:3),
              paste0("Relative_M", 1:3))
  if (!all(needed %in% names(table$wide))) stop("Missing EA RMSFE columns.")
  x <- table$wide
  x <- x[order(x$period), , drop = FALSE]
  fmt <- function(value) sprintf("%.3f", value)
  row_end <- " \\\\"
  rows <- vapply(seq_len(nrow(x)), function(i) {
    values <- unlist(x[i, needed], use.names = FALSE)
    paste0(as.character(x$period[i]), " & ",
           paste(fmt(as.numeric(values)), collapse = " & "), row_end)
  }, character(1))

  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\scriptsize",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\setlength{\\tabcolsep}{2.5pt}",
    "\\definecolor{topgray}{gray}{0.92}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{0.88\\textwidth}{!}{%",
    "\\begin{tabular}{lccc@{\\hspace{0.40cm}}ccc@{\\hspace{0.40cm}}ccc}",
    "\\toprule",
    paste0(" & \\multicolumn{3}{c}{\\textbf{Matrix MF--TPRF}}",
           " & \\multicolumn{3}{c}{\\textbf{VEC-C}}",
           " & \\multicolumn{3}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-C}}$}} ", row_end),
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
    "\\rowcolor{topgray}",
    paste0("\\textbf{Period}",
           " & ", paste(rep("\\textbf{M1} & \\textbf{M2} & \\textbf{M3}", 3L),
                          collapse = " & "), row_end),
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.08cm}",
    paste0("\\parbox{0.90\\linewidth}{\\centering\\footnotesize ",
           "\\textit{Notes:} RMSFEs are computed on common evaluable EA ",
           "forecasts. Ratios below one favour Matrix MF--TPRF.}"),
    "\\endgroup",
    "\\end{table}"
  )
}

# Original EA figure design: three vintage panels, original red/green/blue
# palette and COVID band. The x-axis now displays EVERY year (2017--2026).
make_ea_nowcast_plot <- function(data, cfg, base_size = 13) {
  d <- data$plot
  dates <- seq(as.Date(paste0(format(as.Date(cfg$start_eval), "%Y"), "-01-01")),
               as.Date(paste0(format(as.Date(cfg$end_eval), "%Y"), "-01-01")),
               by = "1 year")
  limits_y <- get_ea_gdp_ylim(data$gdp, pad_frac = 0.15, min_pad = 0.05)

  ggplot2::ggplot(
    d, ggplot2::aes(x = date, y = value, colour = series,
                    linetype = series, group = interaction(type, series))
  ) +
    ggplot2::annotate("rect", xmin = as.Date(cfg$covid_start),
                      xmax = as.Date(cfg$covid_end),
                      ymin = -Inf, ymax = Inf, fill = "grey72", alpha = 0.11) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey55") +
    ggplot2::geom_line(linewidth = 1.00, alpha = 0.98,
                       lineend = "round", na.rm = TRUE) +
    ggplot2::facet_wrap(~ type, ncol = 3) +
    ggplot2::scale_colour_manual(
      values = c("EA GDP" = "#C00000",
                 "Matrix MF-TPRF" = "#228B22", "VEC-C" = "#1F77B4"),
      breaks = c("EA GDP", "Matrix MF-TPRF", "VEC-C"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c("EA GDP" = "solid",
                 "Matrix MF-TPRF" = "solid", "VEC-C" = "33"),
      breaks = c("EA GDP", "Matrix MF-TPRF", "VEC-C"),
      name = NULL
    ) +
    ggplot2::scale_x_date(
      limits = c(as.Date(cfg$start_eval), as.Date(cfg$end_eval)),
      breaks = dates, date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      limits = limits_y, oob = scales::squish,
      breaks = pretty(limits_y, n = 6),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(expand = FALSE, clip = "on") +
    ggplot2::labs(title = "Expanding pseudo-real-time EA nowcasts",
                  x = NULL, y = "EA GDP growth") +
    theme_country_compare(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )
  }

save_ea_nowcast_outputs <- function(data, table, plot,
    output_dir = "results/paper/appendix",
    figure_stem = "figure_ea_nowcasts_matrix_vecc",
    table_stem = "table_ea_rmsfe_matrix_vecc") {
  figure_dir <- file.path(output_dir, "figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  tex_path <- file.path(table_dir, paste0(table_stem, ".tex"))
  writeLines(make_ea_rmsfe_latex(table), tex_path, useBytes = TRUE)
  utils::write.csv(table$long,
                   file.path(table_dir, paste0(table_stem, "_detail.csv")),
                   row.names = FALSE)
  utils::write.csv(table$wide,
                   file.path(table_dir, paste0(table_stem, ".csv")),
                   row.names = FALSE)
  utils::write.csv(data$matched,
                   file.path(table_dir, "ea_nowcasts_common_evaluation.csv"),
                   row.names = FALSE)
  utils::write.csv(data$coverage,
                   file.path(table_dir, "ea_nowcasts_coverage.csv"),
                   row.names = FALSE)

  for (ext in c("pdf", "png")) {
    path <- file.path(figure_dir, paste0(figure_stem, ".", ext))
    ggplot2::ggsave(path, plot = plot, width = 12.5, height = 4.8,
                    units = "in", dpi = 500, bg = "white")
    message("Saved: ", path)
  }
  message("Saved: ", tex_path)
  invisible(list(figure_dir = figure_dir, table_dir = table_dir,
                 tex = tex_path))
}
