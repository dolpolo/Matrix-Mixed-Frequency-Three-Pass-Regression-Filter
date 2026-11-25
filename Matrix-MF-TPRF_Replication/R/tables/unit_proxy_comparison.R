# ==============================================================================
# Unit-target proxy robustness: RMSFE and relative performance
# Read only saved national nowcasts; do not re-estimate any model.
# Uses original Matrix MF-TPRF evaluation and Diebold-Mariano helpers.
# ==============================================================================

unit_comparison_countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
unit_comparison_vintages <- c("M1", "M2", "M3")
unit_comparison_periods <- c("Pre-COVID", "COVID period", "Post-COVID")

load_unit_comparison_results <- function(size = "fixed50", selection = "fixed") {
  required <- c("find_result_file", "find_dfm_summary_file")
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop("Source R/matrix/matrix.mf.tprf.utils.R first: ",
         paste(missing, collapse = ", "))
  }

  find_matrix <- function(model) {
    tryCatch(
      find_result_file(path = "results/empirical/matrix", model = model,
                       stage = "rt", Size = size, sel = selection),
      error = function(e) {
        stop("Missing ", model, " pseudo-real-time RDS (", size, "/", selection,
             "). The current 01_matrix.R defaults to proxy_mode = 'scalar'. ",
             "To obtain unit-target results, run 01_matrix.R with ",
             "proxy_mode = 'multivariate' on the current information set.\n",
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  files <- list(
    unit   = find_matrix("matrix_multivariate"),
    scalar = find_matrix("matrix_scalar"),
    dfm    = find_dfm_summary_file(path = "results/empirical/dfm/q0",
                                   Size = size, sel = selection)
  )
  out <- lapply(files, readRDS)
  for (nm in c("unit", "scalar")) {
    expected <- if (nm == "unit") "multivariate" else "scalar"
    if (!identical(as.character(out[[nm]]$proxy_mode), expected) ||
        !identical(as.character(out[[nm]]$stage), "pseudo_realtime") ||
        !identical(as.character(out[[nm]]$Size), size) ||
        !identical(as.character(out[[nm]]$sel), selection)) {
      stop("Invalid saved ", nm, " Matrix result: ", files[[nm]])
    }
  }
  if (!is.data.frame(out$dfm$df_rt_all) ||
      !is.data.frame(out$dfm$df_yq_eval_all)) {
    stop("DFM summary lacks df_rt_all or df_yq_eval_all: ", files$dfm)
  }
  if (!is.null(out$dfm$idio_spec) &&
      !identical(as.character(out$dfm$idio_spec), "q0")) {
    stop("Expected DFM with q0 idiosyncratic specification, found: ",
         as.character(out$dfm$idio_spec))
  }
  attr(out, "source_files") <- files
  out
}

unit_comparison_clean_eval <- function(x, model, countries,
                                       start_eval, end_eval) {
  required <- c("date", "country", "quarter_id", "type", "period", "GDP", "se")
  absent <- setdiff(required, names(x))
  if (length(absent)) stop("Missing evaluation columns for ", model, ": ",
                          paste(absent, collapse = ", "))
  x <- dplyr::transmute(
    x, date = as.Date(date), country = as.character(country),
    quarter_id = as.character(quarter_id),
    type = as.character(type), period = as.character(period),
    GDP = as.numeric(GDP), se = as.numeric(se)
  )
  x <- dplyr::filter(
    x, country %in% countries,
    type %in% unit_comparison_vintages,
    period %in% unit_comparison_periods,
    date >= as.Date(start_eval), date <= as.Date(end_eval),
    is.finite(GDP), is.finite(se)
  )
  key <- paste(x$country, x$quarter_id, x$type, sep = "|")
  if (anyDuplicated(key)) {
    stop("Duplicated country/quarter/vintage forecasts in ", model)
  }
  x
}

make_unit_comparison_data <- function(results,
                                      countries = unit_comparison_countries,
                                      start_eval = as.Date("2017-01-01"),
                                      end_eval = as.Date("2026-02-01")) {
  needed <- c("build_unit_proxy_evaluation", "build_eval_df", "dm_test_hac")
  if (!all(vapply(needed, exists, logical(1), mode = "function"))) {
    stop("The Matrix evaluation and DM helpers must be sourced first.")
  }

  eval_u <- build_unit_proxy_evaluation(results$unit, countries)
  eval_s <- build_unit_proxy_evaluation(results$scalar, countries)
  eval_d <- build_eval_df(results$dfm$df_rt_all,
                          results$dfm$df_yq_eval_all)
  eval_u <- unit_comparison_clean_eval(eval_u, "unit", countries,
                                        start_eval, end_eval)
  eval_s <- unit_comparison_clean_eval(eval_s, "scalar", countries,
                                        start_eval, end_eval)
  eval_d <- unit_comparison_clean_eval(eval_d, "DFM", countries,
                                        start_eval, end_eval)

  keys <- c("country", "quarter_id", "type")
  u <- dplyr::select(eval_u, dplyr::all_of(keys),
                     period_unit = period, GDP_unit = GDP, se_unit = se)
  s <- dplyr::select(eval_s, dplyr::all_of(keys),
                     period_scalar = period, GDP_scalar = GDP,
                     se_scalar = se)
  d <- dplyr::select(eval_d, dplyr::all_of(keys),
                     period_dfm = period, GDP_dfm = GDP, se_dfm = se)

  common <- dplyr::inner_join(dplyr::inner_join(u, s, by = keys), d, by = keys)
  if (!nrow(common)) stop("No common evaluable unit/scalar/DFM forecasts.")
  if (any(common$period_unit != common$period_scalar |
          common$period_unit != common$period_dfm)) {
    stop("Regime assignments differ across the three result files.")
  }
  if (any(abs(common$GDP_unit - common$GDP_scalar) > 1e-10 |
          abs(common$GDP_unit - common$GDP_dfm) > 1e-10)) {
    stop("GDP realizations differ across unit, scalar and DFM results. ",
         "Verify targets and vintages before comparing RMSFEs.")
  }
  common <- dplyr::transmute(
    common, country, quarter_id, type, period = period_unit,
    GDP = GDP_unit, se_unit, se_scalar, se_dfm
  )
  common <- dplyr::arrange(
    common, country, period,
    as.integer(substr(quarter_id, 1L, 4L)),
    as.integer(sub("^.*Q", "", quarter_id)), type
  )

  # All three models use exactly the SAME country/quarter/vintage observations.
  table <- common |>
    dplyr::group_by(country, period, type) |>
    dplyr::group_modify(~ {
      x <- .x
      rmse_u <- sqrt(mean(x$se_unit))
      rmse_s <- sqrt(mean(x$se_scalar))
      rmse_d <- sqrt(mean(x$se_dfm))
      dm_s <- dm_test_hac(x$se_unit - x$se_scalar, alternative = "less")
      dm_d <- dm_test_hac(x$se_unit - x$se_dfm, alternative = "less")
      data.frame(
        n = nrow(x), rmsfe_unit = rmse_u,
        rmsfe_scalar = rmse_s, rmsfe_dfm = rmse_d,
        ratio_scalar = if (rmse_s > 0) rmse_u / rmse_s else NA_real_,
        ratio_dfm = if (rmse_d > 0) rmse_u / rmse_d else NA_real_,
        dm_scalar_stat = dm_s$stat, dm_scalar_p = dm_s$p_value,
        dm_dfm_stat = dm_d$stat, dm_dfm_p = dm_d$p_value
      )
    }) |>
    dplyr::ungroup()

  list(table = as.data.frame(table), common = as.data.frame(common),
       n_unit = nrow(eval_u), n_scalar = nrow(eval_s), n_dfm = nrow(eval_d))
}

unit_comparison_latex <- function(data,
                                  countries = unit_comparison_countries,
                                  digits = 3,
                                  caption = "Unit-target proxy design: RMSFE and relative performance",
                                  label = "tab:unit_target_comparison") {
  tbl <- data$table
  rowbreak <- " \\\\"
  fmt <- function(x, ratio = FALSE, p = NA_real_) {
    if (length(x) != 1L || !is.finite(x)) return("--")
    out <- sprintf(paste0("%.", digits, "f"), x)
    if (ratio && x < 1 && is.finite(p)) {
      stars <- if (p < 0.01) "***" else if (p < 0.05) "**" else
               if (p < 0.10) "*" else ""
      # The superscript has zero width: stars cannot shift or misalign
      # the numerical entries, even when adjacent rows have no stars.
      if (nzchar(stars)) out <- paste0(out, "\\unitstars{", stars, "}")
    }
    out
  }
  output <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\begingroup",
    "\\fontsize{8.4}{9.1}\\selectfont",
    "\\renewcommand{\\arraystretch}{1.02}",
    "\\setlength{\\tabcolsep}{0pt}",
    # Fixed-width centered columns, as in the DFM comparison table.
    # Together with zero-width stars, numeric centering is unchanged.
    "\\newcommand{\\unitstars}[1]{\\rlap{\\textsuperscript{\\fontsize{4.2}{4.6}\\selectfont #1}}}",
    "\\definecolor{matrixgray}{gray}{0.965}",
    "\\definecolor{regimegray}{gray}{0.94}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{0.98\\textwidth}{!}{%",
    "\\begin{tabular}{@{}>{\\centering\\arraybackslash}p{1.35cm}@{\\hspace{0.18cm}}",
    "*{3}{>{\\columncolor{matrixgray}\\centering\\arraybackslash}p{1.40cm}}@{\\hspace{0.18cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{1.40cm}}@{\\hspace{0.18cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{1.40cm}}@{}}",
    "\\toprule",
    paste0("& \\multicolumn{3}{c}{\\textbf{Unit-target proxy design}}",
           " & \\multicolumn{3}{c}{\\textbf{$\\dfrac{\\text{Unit-target}}{\\text{Aggregate-target}}$}}",
           " & \\multicolumn{3}{c}{\\textbf{$\\dfrac{\\text{Unit-target}}{\\text{DFM}}$}}", rowbreak),
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
    paste0("\\textbf{Country} & M1 & M2 & M3 & M1 & M2 & M3",
           " & M1 & M2 & M3", rowbreak),
    "\\midrule"
  )

  regime_names <- c("Pre-Covid", "Covid", "Post-Covid")
  for (j in seq_along(unit_comparison_periods)) {
    period <- unit_comparison_periods[j]
    if (j != 1L) output <- c(output, "\\addlinespace[0.20em]")
    output <- c(output,
      "\\rowcolor{regimegray}",
      paste0("\\multicolumn{10}{@{}l}{\\textbf{", regime_names[j],
             "}}", rowbreak),
      "\\cmidrule{1-10}")
    for (cc in countries) {
      cells <- character(9)
      for (k in seq_along(unit_comparison_vintages)) {
        v <- unit_comparison_vintages[k]
        r <- tbl[tbl$country == cc & tbl$period == period & tbl$type == v,
                 , drop = FALSE]
        if (nrow(r) > 1L) stop("Duplicate summary for ", cc, "/", period, "/", v)
        if (nrow(r)) {
          cells[k]     <- fmt(r$rmsfe_unit)
          cells[k + 3] <- fmt(r$ratio_scalar, TRUE, r$dm_scalar_p)
          cells[k + 6] <- fmt(r$ratio_dfm, TRUE, r$dm_dfm_p)
        } else {
          cells[c(k, k + 3, k + 6)] <- "--"
        }
      }
      output <- c(output, paste0(cc, " & ", paste(cells, collapse = " & "),
                                  rowbreak))
    }
  }
  output <- c(output,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\par\\vspace{0.10cm}",
    paste0("\\parbox{0.96\\textwidth}{\\footnotesize ",
           "\\textit{Notes:} The shaded block reports RMSFEs for the ",
           "unit-target proxy design. The remaining blocks report RMSFE ratios ",
           "relative to the aggregate-target proxy design and DFM. ",
           "All statistics use common evaluable forecasts; ratios below one ",
           "favour unit-target. Asterisks indicate one-sided Diebold--Mariano ",
           "significance using a Newey--West HAC variance estimator ",
           "($^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$).}"),
    "\\endgroup",
    "\\end{table}"
  )
  paste(output, collapse = "\n")
}

save_unit_comparison <- function(result,
                                 output_dir = "results/paper/appendix/tables",
                                 stem = "table_unit_target_comparison") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path_tex <- file.path(output_dir, paste0(stem, ".tex"))
  path_csv <- file.path(output_dir, paste0(stem, ".csv"))
  path_dm  <- file.path(output_dir, paste0(stem, "_dm_tests.csv"))
  path_obs <- file.path(output_dir, paste0(stem, "_common_forecasts.csv"))
  writeLines(unit_comparison_latex(result), path_tex, useBytes = TRUE)
  utils::write.csv(result$table[, c("country", "period", "type", "n",
    "rmsfe_unit", "rmsfe_scalar", "rmsfe_dfm", "ratio_scalar", "ratio_dfm")],
    path_csv, row.names = FALSE, na = "")
  utils::write.csv(result$table[, c("country", "period", "type", "n",
    "dm_scalar_stat", "dm_scalar_p", "dm_dfm_stat", "dm_dfm_p")],
    path_dm, row.names = FALSE, na = "")
  utils::write.csv(result$common, path_obs, row.names = FALSE)
  message("Saved: ", path_tex)
  message("Saved: ", path_csv)
  message("Saved: ", path_dm)
  message("Saved: ", path_obs)
  invisible(c(path_tex, path_csv, path_dm, path_obs))
}
