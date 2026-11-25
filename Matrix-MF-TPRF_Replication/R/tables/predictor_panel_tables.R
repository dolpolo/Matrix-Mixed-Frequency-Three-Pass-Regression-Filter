# ==============================================================================
# Appendix | Predictor panel tables
# Paper: Matrix Three-Pass Regression Filter
# Functions only: no execution, no model estimation, no predictor selection.
# Input: information-set CSV, national transformed workbooks, existing prep tools.
# ==============================================================================

# Internal helper: read and validate the 50-indicator dataset specification.
read_predictor_panel_spec <- function(
    csv_path = "data/information_set/selected_variables_50.csv") {
  if (!file.exists(csv_path)) stop("Missing indicator file: ", csv_path)
  selected <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  if (!all(c("base_name", "frequency") %in% names(selected)) ||
      nrow(selected) != 50L || anyNA(selected$base_name) ||
      anyDuplicated(selected$base_name) ||
      !all(selected$frequency %in% c("M", "Q")) ||
      sum(selected$frequency == "M") != 35L ||
      sum(selected$frequency == "Q") != 15L) {
    stop("Expected 50 distinct indicators: 35 monthly and 15 quarterly.")
  }
  selected
}

# Internal helper: use exactly the grouping and metadata of the old table.
get_predictor_panel_metadata <- function(selected) {
  if (!exists("build_metadata_table", mode = "function")) {
    stop("Please source R/matrix/matrix.mf.tprf.utils.R first.")
  }
  metadata <- as.data.frame(build_metadata_table(), stringsAsFactors = FALSE)
  idx <- match(selected$base_name, metadata$ID)
  if (anyNA(idx) || any(metadata$Fr[idx] != selected$frequency)) {
    stop("Some indicators have missing or inconsistent original metadata.")
  }
  x <- metadata[idx, c("ID", "Series", "Cl", "Cat", "Tr", "Fr", "Del", "Group")]
  groups <- c(
    "National Accounts / Real Economy", "Labor Market", "Credit Aggregates",
    "Labor Costs", "Financial Markets", "Interest Rates",
    "Industrial Production and Turnover", "Prices", "Confidence Indicators"
  )
  if (anyNA(match(x$Group, groups))) {
    stop("Unknown variable group in original metadata.")
  }
  x <- x[order(match(x$Group, groups), x$Fr, x$ID), , drop = FALSE]
  x$N <- seq_len(nrow(x))
  list(metadata = x, groups = groups)
}

# Internal helper: determine whether each national series exists in its workbook.
# Deliberately NOT based on the historical S/M/L selection cells.
get_country_indicator_coverage <- function(metadata, countries,
                                           transformed_dir = "data/transformed") {
  out <- matrix(FALSE, nrow = nrow(metadata), ncol = length(countries),
                dimnames = list(metadata$ID, countries))
  for (cc in countries) {
    files <- setNames(
      file.path(transformed_dir, paste0(cc, "data", c("M", "Q"), "_NA_TR2.xlsx")),
      c("M", "Q")
    )
    if (!all(file.exists(files))) {
      stop("Missing transformed workbook for ", cc, ": ",
           paste(files[!file.exists(files)], collapse = ", "))
    }
    headers <- lapply(files, function(f) names(readxl::read_excel(f, n_max = 0L)))
    out[, cc] <- vapply(seq_len(nrow(metadata)), function(j) {
      paste0(metadata$ID[j], "_", cc) %in% headers[[metadata$Fr[j]]]
    }, logical(1))
  }
  out
}

# ------------------------------------------------------------------------------
# TABLE 1: Composition and country coverage of the information set
# Same longtable layout, metadata, groups and dimensions as the paper's old table.
# ------------------------------------------------------------------------------
make_information_set_table <- function(
    countries = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
    csv_path = "data/information_set/selected_variables_50.csv",
    transformed_dir = "data/transformed",
    output_file = "results/paper/appendix/tables/information_set.tex") {

  selected <- read_predictor_panel_spec(csv_path)
  original <- get_predictor_panel_metadata(selected)
  tab <- original$metadata
  coverage <- get_country_indicator_coverage(tab, countries, transformed_dir)
  unavailable <- colSums(!coverage)

  # R string with four backslashes generates the required LaTeX two-backslash newline.
  row_end <- " \\\\"
  header_names <- c("N", "ID", "Series", "Cl.", "Cat.", "Tr.", "Fr.", "Del.", countries)
  header <- paste0(paste(sprintf("\\textbf{%s}", header_names), collapse = " & "), row_end)

  lines <- c(
    "%=========================================================",
    "% DATASET: INDICATOR COMPOSITION AND COUNTRY AVAILABILITY",
    "%=========================================================",
    "\\begingroup", "\\centering", "\\fontsize{7.7}{8.1}\\selectfont",
    "\\providecommand{\\setcell}[1]{#1}",
    "\\renewcommand{\\setcell}[1]{{\\fontsize{4.8}{6.1}\\selectfont #1}}",
    "\\renewcommand{\\arraystretch}{0.60}",
    "\\setlength{\\extrarowheight}{0pt}",
    "\\setlength{\\tabcolsep}{1.6pt}",
    "\\setlength{\\LTleft}{0pt plus 1fill}",
    "\\setlength{\\LTright}{0pt plus 1fill}",
    "\\setlength{\\LTpre}{0pt}",
    "\\setlength{\\LTpost}{0pt}",
    "\\begin{longtable}{@{}c p{1.55cm} p{7.90cm} *{13}{c}@{}}",
    paste0("\\caption{Composition and country coverage of the information set} ",
           "\\label{tab:fixed_information_set}", row_end),
    "\\toprule", header, "\\midrule \\midrule \\endfirsthead",
    paste0("\\multicolumn{16}{c}{\\scriptsize\\textit{Table \\thetable\\ continued}}", row_end),
    "\\toprule", header, "\\midrule \\midrule \\endhead",
    paste0("\\midrule \\multicolumn{16}{r}{\\scriptsize\\textit{Continued on next page}}", row_end),
    "\\endfoot", "\\bottomrule \\endlastfoot"
  )

  for (group in original$groups) {
    idx <- which(tab$Group == group)
    if (!length(idx)) next
    lines <- c(lines,
      paste0("\\multicolumn{16}{c}{\\textbf{", escape_latex(group), "}}", row_end),
      "\\midrule")
    for (j in idx) {
      cells <- c(
        as.character(tab$N[j]), escape_latex(tab$ID[j]),
        escape_latex(tab$Series[j]), as.character(tab$Cl[j]),
        as.character(tab$Cat[j]), as.character(tab$Tr[j]),
        as.character(tab$Fr[j]), as.character(tab$Del[j]),
        vapply(countries, function(cc) {
          if (coverage[j, cc]) "\\setcell{$\\checkmark$}" else "-"
        }, character(1))
      )
      lines <- c(lines, paste0(paste(cells, collapse = " & "), row_end))
    }
    lines <- c(lines, "\\midrule")
  }

  lines <- c(
    lines,
    paste0("\\multicolumn{8}{r}{\\textbf{Unavailable series}} & ",
           paste(sprintf("\\textbf{%d}", unavailable[countries]), collapse = " & "),
           row_end),
    "\\end{longtable}",
    "\\vspace{-0.40cm}",
    "\\begin{center}", "\\parbox{0.98\\textwidth}{%",
    "\\fontsize{7}{5}\\selectfont",
    paste0(
      "\\textit{Notes:} The table reports the 50 monthly and quarterly indicators ",
      "in the information set. A checkmark denotes an available country--indicator ",
      "series; a dash denotes an unavailable series. The final row counts unavailable ",
      "series by country. ",
      "\\textbf{Cl.}: $R$ = real activity, $N$ = nominal, $C$ = confidence, ",
      "$F$ = financial. \\textbf{Cat.}: $H$ = hard, $S$ = soft. ",
      "\\textbf{Tr.}, \\textbf{Fr.}, and \\textbf{Del.} denote transformation, ",
      "frequency, and approximate publication delay."
    ),
    "}", "\\end{center}", "\\endgroup"
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, output_file, useBytes = TRUE)
  invisible(list(path = output_file, missing_by_country = unavailable,
                 coverage = coverage, metadata = tab))
}

# ------------------------------------------------------------------------------
# TABLE 2: Predictor unavailability by period and frequency
# Compact: three periods x Monthly / Quarterly / Overall.
# Averaged over monthly vintages M1, M2, M3, before imputation.
# ------------------------------------------------------------------------------
compute_predictor_unavailability <- function(
    config,
    raw_dir = "data/raw",
    transformed_dir = "data/transformed") {

  params <- list(
    start_est = config$start_est,
    start_eval = config$start_eval,
    end_eval = config$end_eval,
    covid_start = config$covid_start,
    covid_end = config$covid_end,
    covid_mask_m = TRUE,
    covid_mask_q = TRUE,
    target = "GDP",
    target_cc = "EA",
    sel_method = "fixed",
    n_m = 20L,
    n_q = 5L
  )

  prep <- prepare_all_countries(
    countries = c(config$countries, params$target_cc),
    params = params, path_raw = raw_dir, path_adj = transformed_dir,
    covid_mask_m = params$covid_mask_m,
    covid_mask_q = params$covid_mask_q
  )
  tensor <- build_tensor(prep = prep, params = params, var_scope = "union")
  panel <- make_regime_data(tensor, label = "dataset",
                            selection_end = params$end_eval, params = params)

  d <- dim(panel$X_full)
  if (length(d) != 3L || d[2] != length(config$countries) ||
      d[3] != 50L || panel$N_m != 35L || panel$N_q != 15L) {
    stop("Unexpected predictor panel dimensions; expected 8 x 50 with 35M+15Q.")
  }

  dates_m <- as.Date(tensor$dates)
  eval_t <- which(dates_m >= params$start_eval & dates_m <= params$end_eval)
  pct_na <- function(x) 100 * sum(is.na(x)) / length(x)

  stats_by_month <- lapply(eval_t, function(t) {
    # A nowcast vintage is the first/second/third calendar month in its quarter.
    # Do not infer it from OBSERVED GDP dates: Covid masking removes GDP target
    # observations and would incorrectly suppress whole evaluation months.
    vintage <- (as.integer(format(dates_m[t], "%m")) - 1L) %% 3L + 1L
    raw_vintage <- unbalancedness_tensor(panel$X_full, panel$Unb, current_t = t)
    monthly <- raw_vintage[, , seq_len(panel$N_m), drop = FALSE]
    quarterly <- raw_vintage[, , panel$N_m + seq_len(panel$N_q), drop = FALSE]
    date <- dates_m[t]
    period <- if (date < params$covid_start) "Pre-Covid" else if (
      date <= params$covid_end) "Covid" else "Post-Covid"

    data.frame(
      date = date, period = period, vintage = paste0("M", vintage),
      monthly = pct_na(monthly), quarterly = pct_na(quarterly),
      overall = pct_na(raw_vintage)
    )
  })
  stats_by_month <- Filter(Negate(is.null), stats_by_month)
  if (!length(stats_by_month)) stop("No evaluation vintages available.")
  stats <- do.call(rbind, stats_by_month)

  periods <- c("Pre-Covid", "Covid", "Post-Covid")
  if (!all(periods %in% stats$period)) {
    stop("At least one evaluation period has no observed vintages.")
  }
  summary <- do.call(rbind, lapply(periods, function(p) {
    xx <- stats[stats$period == p, , drop = FALSE]
    data.frame(period = p,
               monthly = mean(xx$monthly),
               quarterly = mean(xx$quarterly),
               overall = mean(xx$overall),
               n_vintages = nrow(xx),
               stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL
  list(summary = summary, per_vintage = stats)
}

make_predictor_unavailability_table <- function(
    config,
    output_file = "results/paper/appendix/tables/realtime_missingness.tex",
    raw_dir = "data/raw",
    transformed_dir = "data/transformed") {

  calc <- compute_predictor_unavailability(config, raw_dir, transformed_dir)
  values <- calc$summary
  row_end <- " \\\\"
  fmt <- function(x) sprintf("%.1f", x)
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering", "\\small",
    "\\renewcommand{\\arraystretch}{1.12}",
    "\\setlength{\\tabcolsep}{12pt}",
    "\\definecolor{topgray}{gray}{0.92}",
    "\\caption{Predictor unavailability by period and frequency}",
    "\\label{tab:realtime_missingness}",
    "\\begin{tabular}{lccc}",
    "\\toprule", "\\rowcolor{topgray}",
    paste0("\\textbf{Period} & \\textbf{Monthly} & \\textbf{Quarterly} & \\textbf{Overall}", row_end),
    "\\midrule"
  )
  for (i in seq_len(nrow(values))) {
    lines <- c(lines, paste0(
      values$period[i], " & ",
      paste(fmt(c(values$monthly[i], values$quarterly[i], values$overall[i])),
            collapse = " & "), row_end
    ))
  }
  lines <- c(
    lines, "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.1cm}", "\\begin{minipage}{0.85\\textwidth}",
    "\\footnotesize",
    "\\textit{Notes:} Missing observations (\\%) before imputation, averaged ",
    "over M1--M3 and the eight countries. Quarterly indicators are ",
    "represented on the monthly grid.",
    "\\end{minipage}", "\\end{table}"
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, output_file, useBytes = TRUE)
  invisible(list(path = output_file, summary = values, per_vintage = calc$per_vintage))
}
