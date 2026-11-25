# ==============================================================================
# DFM versus VEC-P and VEC-C | Appendix country-level RMSFE table
# Author: Davide Delfino
# Input: saved q0 DFM, VEC-P and VEC-C national summary RDS files
# Output: a LaTeX table and auditable CSV data; no model is re-estimated
# ==============================================================================
# Requires the helper functions in R/matrix/matrix.mf.tprf.utils.R:
#   get_rt_table_for_final, get_hyper_country_table, hp_cells_dfm,
#   hp_cells_tprf, build_eval_df, build_dm_wide, build_period_labels.
# Requires dplyr, tidyr, lubridate and sandwich (as used by those helpers).

load_dfm_benchmark_summaries <- function() {
  dfm_path <- find_dfm_summary_file(
    path = "results/empirical/dfm/q0", Size = "fixed50", sel = "fixed"
  )
  vecp_path <- find_result_file(
    path = "results/empirical/pooled", model = "vectensor",
    stage = "summary", Size = "fixed50", sel = "fixed"
  )
  vecc_path <- find_result_file(
    path = "results/empirical/vector", model = "vector",
    stage = "summary", Size = "fixed50", sel = "fixed"
  )
  summaries <- list(dfm = readRDS(dfm_path),
                    vecp = readRDS(vecp_path),
                    vecc = readRDS(vecc_path))
  attr(summaries, "paths") <- c(dfm = dfm_path,
                                  vecp = vecp_path,
                                  vecc = vecc_path)
  summaries
}

prepare_dfm_benchmark_comparison <- function(summaries,
    country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
    period_order = c("Pre-COVID", "COVID period", "Post-COVID")) {

  stopifnot(all(c("dfm", "vecp", "vecc") %in% names(summaries)))

  get_model_table <- function(summary_obj, model_tag) {
    # These are the original full-period, model-specific rolling RMSFEs,
    # consistent with the existing comparison-table convention.
    x <- get_rt_table_for_final(summary_obj, country_order, period_order)
    required <- c("country", "period", "M1", "M2", "M3")
    if (!all(required %in% names(x))) {
      stop(model_tag, ": incomplete tab_rt_all. Missing: ",
           paste(setdiff(required, names(x)), collapse = ", "))
    }
    x <- dplyr::transmute(x,
       country = as.character(country), period = as.character(period),
       M1 = as.numeric(M1), M2 = as.numeric(M2), M3 = as.numeric(M3))
    if (anyDuplicated(x[c("country", "period")])) {
      stop(model_tag, ": duplicated country-period in tab_rt_all.")
    }
    names(x)[3:5] <- paste0(model_tag, "_", c("M1", "M2", "M3"))
    x
  }

  dfm_rt <- get_model_table(summaries$dfm, "DFM")
  vp_rt  <- get_model_table(summaries$vecp, "VecP")
  vc_rt  <- get_model_table(summaries$vecc, "VecC")

  all_keys <- expand.grid(country = country_order, period = period_order,
                          stringsAsFactors = FALSE)
  table_data <- all_keys |>
    dplyr::left_join(dfm_rt, by = c("country", "period")) |>
    dplyr::left_join(vp_rt, by = c("country", "period")) |>
    dplyr::left_join(vc_rt, by = c("country", "period"))

  for (m in c("M1", "M2", "M3")) {
    den_p <- table_data[[paste0("VecP_", m)]]
    den_c <- table_data[[paste0("VecC_", m)]]
    num   <- table_data[[paste0("DFM_", m)]]
    table_data[[paste0("Rel_VecP_", m)]] <-
      ifelse(is.finite(den_p) & den_p > 0, num / den_p, NA_real_)
    table_data[[paste0("Rel_VecC_", m)]] <-
      ifelse(is.finite(den_c) & den_c > 0, num / den_c, NA_real_)
  }

  # The benchmark's published GDP observations are used as the common target.
  # build_eval_df aligns forecasts by country, quarter and within-quarter month.
  if (is.null(summaries$dfm$df_yq_eval_all)) {
    stop("DFM summary is missing df_yq_eval_all (observed GDP).")
  }
  y <- summaries$dfm$df_yq_eval_all
  e_dfm <- build_eval_df(summaries$dfm$df_rt_all, y)
  e_vp  <- build_eval_df(summaries$vecp$df_rt_all, y)
  e_vc  <- build_eval_df(summaries$vecc$df_rt_all, y)

  # d = SE_DFM - SE_benchmark; the one-sided alternative is d < 0.
  # The Newey-West bandwidth and finite-sample correction are inherited
  # from the author's existing dm_test_hac() implementation.
  dm_vecp <- build_dm_wide(e_dfm, e_vp, "VecP", alternative = "less")
  dm_vecc <- build_dm_wide(e_dfm, e_vc, "VecC", alternative = "less")

  table_data <- table_data |>
    dplyr::left_join(dm_vecp$wide |>
       dplyr::mutate(country = as.character(country),
                     period = as.character(period)),
       by = c("country", "period")) |>
    dplyr::left_join(dm_vecc$wide |>
       dplyr::mutate(country = as.character(country),
                     period = as.character(period)),
       by = c("country", "period"))

  for (tag in c("VecP", "VecC")) {
    for (m in c("M1", "M2", "M3")) {
      pcol <- paste0("p_", tag, "_", m)
      if (!pcol %in% names(table_data)) table_data[[pcol]] <- NA_real_
      rcol <- paste0("Rel_", tag, "_", m)
      scol <- paste0("star_", tag, "_", m)
      p <- table_data[[pcol]]
      ratio <- table_data[[rcol]]
      table_data[[scol]] <- ifelse(
        is.finite(ratio) & ratio < 1 & is.finite(p),
        ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                                    ifelse(p < 0.10, "*", ""))), ""
      )
    }
  }

  hp <- list(
    dfm_pre  = get_hyper_country_table(summaries$dfm,  "pre",  country_order),
    dfm_post = get_hyper_country_table(summaries$dfm,  "post", country_order),
    vecp_pre = get_hyper_country_table(summaries$vecp, "pre",  country_order),
    vecp_post= get_hyper_country_table(summaries$vecp, "post", country_order),
    vecc_pre = get_hyper_country_table(summaries$vecc, "pre",  country_order),
    vecc_post= get_hyper_country_table(summaries$vecc, "post", country_order)
  )

  list(data = table_data, hyper = hp,
       dm = list(vecp = dm_vecp, vecc = dm_vecc),
       evaluations = list(dfm = e_dfm, vecp = e_vp, vecc = e_vc))
}

build_dfm_relative_rmsfe_latex <- function(comparison, params,
    caption = "DFM RMSFE and relative RMSFE by country",
    label = "tab:dfm_relative_rmsfe",
    country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")) {

  # The DFM in this replication uses q=0, hence HP = (P_f, r).
  # VEC-P and VEC-C each use HP = (L, P_f, P_rho, r).
  d <- comparison$data
  hp <- comparison$hyper
  periods <- c("Pre-COVID", "COVID period", "Post-COVID")
  period_labels <- build_period_labels(params)

  num <- function(x) {
    if (length(x) == 0L || !is.finite(x[1L])) return("")
    sprintf("%.3f", x[1L])
  }
  ratio_cell <- function(x, star) {
    ans <- num(x)
    if (!nzchar(ans)) return(ans)
    if (length(star) && !is.na(star[1L]) && nzchar(star[1L])) {
      ans <- paste0(ans, "\\dfmsym{", star[1L], "}")
    }
    ans
  }
  hp_clean <- function(x) {
    x[is.na(x) | x == "NA" | x == "NaN"] <- ""
    unname(x)
  }
  hp_for_period <- function(cc, period_name, model) {
    regime <- if (identical(period_name, "Post-COVID")) "post" else "pre"
    hp_df <- hp[[paste0(model, "_", regime)]]
    if (identical(model, "dfm")) {
      hp_clean(hp_cells_dfm(hp_df, cc, include_pq = FALSE))
    } else {
      hp_clean(hp_cells_tprf(hp_df, cc))
    }
  }

  # Exactly 20 columns, following the original 3-period/8-country design:
  # Country | DFM: 3 + 2 HP | DFM/VEC-P: 3 + 4 HP | DFM/VEC-C: 3 + 4 HP.
  spec <- paste0(
    ">{\\centering\\arraybackslash}p{1.05cm}",
    "@{\\hspace{0.16cm}}",
    "*{3}{>{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.82cm}}",
    "*{2}{>{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.13cm}}",
    "@{\\hspace{0.16cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
    "*{4}{>{\\centering\\arraybackslash}p{0.13cm}}",
    "@{\\hspace{0.16cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
    "*{4}{>{\\centering\\arraybackslash}p{0.13cm}}"
  )

  row_end <- " \\\\" # Two LaTeX backslashes, not four.
  period_block <- function(period_name) {
    lab <- period_labels[[period_name]]
    lines <- c(
      "\\specialrule{0.08em}{0.15em}{0.08em}",
      paste0("\\multicolumn{20}{l}{{\\fontsize{5.8}{6.3}\\selectfont\\textbf{",
             lab[["title"]], "}}}\\\\[-0.60em]"),
      paste0("\\multicolumn{20}{l}{{\\fontsize{4.5}{4.8}\\selectfont ",
             lab[["date"]], "}}\\\\[-0.32em]"),
      "\\cmidrule{1-20}",
      "\\rowcolor{hypergray}",
      paste0(
        "\\textit{\\fontsize{5.0}{5.8}\\selectfont Hyper-params}",
        " & \\multicolumn{3}{c}{}",
        " & \\cellcolor{hypergray}\\dfmhdr $P_f$",
        " & \\cellcolor{hypergray}\\dfmhdr $r$",
        " & \\multicolumn{3}{c}{}",
        " & \\cellcolor{hypergray}\\dfmhdr $L$",
        " & \\cellcolor{hypergray}\\dfmhdr $P_f$",
        " & \\cellcolor{hypergray}\\dfmhdr $P_\\rho$",
        " & \\cellcolor{hypergray}\\dfmhdr $r$",
        " & \\multicolumn{3}{c}{}",
        " & \\cellcolor{hypergray}\\dfmhdr $L$",
        " & \\cellcolor{hypergray}\\dfmhdr $P_f$",
        " & \\cellcolor{hypergray}\\dfmhdr $P_\\rho$",
        " & \\cellcolor{hypergray}\\dfmhdr $r$",
        "\\\\[0.25em]"
      )
    )

    for (cc in country_order) {
      row <- d[d$country == cc & d$period == period_name, , drop = FALSE]
      if (nrow(row) != 1L) stop("Expected one row for ", cc, ", ", period_name)
      hp_dfm <- hp_for_period(cc, period_name, "dfm")
      hp_vp  <- hp_for_period(cc, period_name, "vecp")
      hp_vc  <- hp_for_period(cc, period_name, "vecc")

      cells <- c(
        cc,
        paste0("\\dfmcell{", vapply(c("M1", "M2", "M3"),
           function(m) num(row[[paste0("DFM_", m)]]), ""), "}"),
        paste0("\\dfmhp ", hp_dfm),
        vapply(c("M1", "M2", "M3"), function(m) {
          ratio_cell(row[[paste0("Rel_VecP_", m)]],
                     row[[paste0("star_VecP_", m)]])
        }, ""),
        paste0("\\dfmhp ", hp_vp),
        vapply(c("M1", "M2", "M3"), function(m) {
          ratio_cell(row[[paste0("Rel_VecC_", m)]],
                     row[[paste0("star_VecC_", m)]])
        }, ""),
        paste0("\\dfmhp ", hp_vc)
      )
      if (length(cells) != 20L) stop("Internal error: expected 20 table cells.")
      lines <- c(lines, paste0(paste(cells, collapse = " & "), row_end))
    }
    lines
  }

  note <- paste0(
    "The shaded block reports DFM RMSFEs and selected hyperparameters. ",
    "The other blocks report RMSFE ratios relative to VEC-P and VEC-C; ",
    "values below one favour DFM. Asterisks indicate one-sided ",
    "Diebold--Mariano significance in favour of DFM using squared errors ",
    "and a Newey--West HAC variance estimator ",
    "($^{*}\\,p\\!<\\!0.10$, $^{**}\\,p\\!<\\!0.05$, ",
    "$^{***}\\,p\\!<\\!0.01$)."
  )

  tex <- c(
    "%=========================================================",
    "% DFM benchmark; same layout as the earlier VEC-C comparison table",
    "% No information-set-size grouping: results read from saved RDS",
    "%=========================================================",
    "\\begin{table}[!ht]",
    "\\centering",
    "\\begingroup",
    "\\scriptsize",
    "\\renewcommand{\\arraystretch}{1.06}",
    "\\setlength{\\tabcolsep}{1.25pt}",
    "\\newcommand{\\dfmhp}{\\fontsize{4.2}{4.8}\\selectfont}",
    "\\newcommand{\\dfmhdr}{\\fontsize{4.3}{5.0}\\selectfont}",
    "\\newcommand{\\dfmcell}[1]{\\cellcolor{matrixgray}#1}",
    "\\newcommand{\\dfmsym}[1]{\\textsuperscript{\\fontsize{3.8}{3.8}\\selectfont#1}}",
    "\\definecolor{topgray}{gray}{0.86}",
    "\\definecolor{hypergray}{gray}{0.91}",
    "\\definecolor{matrixgray}{gray}{0.965}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", spec, "}"),
    "\\toprule",
    paste0(
      "\\multicolumn{1}{>{\\cellcolor{topgray}}c@{\\hspace{0.16cm}}}{}",
      " & \\multicolumn{5}{@{}>{\\columncolor{matrixgray}}c@{\\hspace{0.16cm}}}{\\textbf{DFM}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{DFM}}{\\text{VEC-P MF--TPRF}}$}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{DFM}}{\\text{VEC-C MF--TPRF}}$}}",
      row_end
    ),
    "\\cmidrule(lr){2-6}\\cmidrule(lr){7-13}\\cmidrule(lr){14-20}",
    paste0(
      "\\textbf{Country}",
      " & \\cellcolor{matrixgray}M1 & \\cellcolor{matrixgray}M2",
      " & \\cellcolor{matrixgray}M3",
      " & \\multicolumn{2}{>{\\columncolor{matrixgray}}c}{HP}",
      " & M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
      " & M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
      row_end
    ),
    period_block("Pre-COVID"),
    period_block("COVID period"),
    period_block("Post-COVID"),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.15cm}",
    paste0("\\parbox{0.96\\textwidth}{\\footnotesize ", note, "}"),
    "\\endgroup",
    "\\end{table}"
  )
  paste(tex, collapse = "\n")
}

save_dfm_benchmark_comparison <- function(comparison, params,
    out_dir = "results/paper/appendix/tables",
    stem = "table_dfm_relative_rmsfe") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_tex <- file.path(out_dir, paste0(stem, ".tex"))
  file_csv <- file.path(out_dir, paste0(stem, ".csv"))
  file_dm  <- file.path(out_dir, paste0(stem, "_dm_tests.csv"))

  writeLines(build_dfm_relative_rmsfe_latex(comparison, params), file_tex,
             useBytes = TRUE)
  utils::write.csv(comparison$data, file_csv, row.names = FALSE, na = "")
  dm1 <- comparison$dm$vecp$long
  dm2 <- comparison$dm$vecc$long
  dm1$benchmark <- "VEC-P"
  dm2$benchmark <- "VEC-C"
  utils::write.csv(dplyr::bind_rows(dm1, dm2), file_dm,
                   row.names = FALSE, na = "")
  invisible(c(latex = file_tex, data = file_csv, dm = file_dm))
}
