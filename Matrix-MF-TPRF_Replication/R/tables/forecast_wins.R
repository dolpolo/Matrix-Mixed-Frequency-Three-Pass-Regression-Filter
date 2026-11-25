# ==============================================================================
# Forecast performance summaries | Matrix MF-TPRF vs VEC-P, VEC-C and DFM
# Functions only: called from scripts/03_outputs/12_forecast_wins_main.R
# Inputs: saved baseline (scalar-proxy) national summary RDS objects.
# No re-estimation and no model-specific data transformations.
# ==============================================================================

fw_country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
fw_vintage_order <- c("M1", "M2", "M3")
fw_period_order <- c("Pre-COVID", "COVID period", "Post-COVID")
fw_benchmarks <- c("VEC-P" = "vecp", "VEC-C" = "vecc", "DFM" = "dfm")

fw_read_summaries <- function(size = "fixed50", selection = "fixed") {
  if (!exists("find_result_file", mode = "function") ||
      !exists("find_dfm_summary_file", mode = "function")) {
    stop("Source R/matrix/matrix.mf.tprf.utils.R before calling fw_read_summaries().")
  }
  read_one <- function(path, model = NULL, dfm = FALSE) {
    file <- if (dfm) {
      find_dfm_summary_file(path = path, Size = size, sel = selection)
    } else {
      find_result_file(path = path, model = model, stage = "summary",
                       Size = size, sel = selection)
    }
    if (length(file) != 1L || !file.exists(file)) {
      stop("Missing or ambiguous saved summary: ", path)
    }
    readRDS(file)
  }
  ans <- list(
    matrix = read_one("results/empirical/matrix", "matrix_scalar"),
    vecp   = read_one("results/empirical/pooled", "vectensor"),
    vecc   = read_one("results/empirical/vector", "vector"),
    dfm    = read_one("results/empirical/dfm/q0", dfm = TRUE)
  )
  if (!identical(as.character(ans$matrix$proxy_mode), "scalar")) {
    stop("The baseline Matrix summary must use the aggregate (scalar) proxy.")
  }
  if (!is.data.frame(ans$matrix$df_yq_eval_all)) {
    stop("The Matrix summary has no df_yq_eval_all GDP evaluation data.")
  }
  invisible(lapply(names(ans), function(k) {
    if (!is.data.frame(ans[[k]]$df_rt_all))
      stop("Missing df_rt_all in ", k, " summary.")
  }))
  ans
}

# The original project helper build_eval_df() attaches the SAME quarterly GDP
# observations/period labels to each model and calculates squared errors.
# Pair strictly on (country, reference quarter, vintage), excluding missing
# squared errors. This makes both country-level and forecast-level comparisons
# internally consistent on the common four-model evaluation sample.
fw_prepare_evaluation <- function(summaries,
                                  countries = fw_country_order,
                                  start_eval = as.Date("2017-01-01"),
                                  end_eval = as.Date("2026-02-01")) {
  if (!exists("build_eval_df", mode = "function")) {
    stop("The original build_eval_df() helper has not been sourced.")
  }
  start_eval <- as.Date(start_eval)
  end_eval <- as.Date(end_eval)
  if (is.na(start_eval) || is.na(end_eval) || start_eval > end_eval)
    stop("Invalid evaluation dates.")
  if (anyDuplicated(countries) || !all(countries %in% fw_country_order))
    stop("Unknown or duplicate country codes.")
  keys <- c("country", "quarter_id", "type")
  assembled <- list()
  counts <- list()

  for (model in c("matrix", "vecp", "vecc", "dfm")) {
    raw <- build_eval_df(summaries[[model]]$df_rt_all,
                         summaries$matrix$df_yq_eval_all)
    need <- c(keys, "date", "period", "se")
    if (!all(need %in% names(raw)))
      stop("Evaluation columns missing for model ", model, ".")
    raw$country <- as.character(raw$country)
    raw$quarter_id <- as.character(raw$quarter_id)
    raw$type <- as.character(raw$type)
    raw$date <- as.Date(raw$date)
    raw$se <- as.numeric(raw$se)
    raw <- raw[
      !is.na(raw$date) & raw$date >= start_eval & raw$date <= end_eval &
        raw$country %in% countries & raw$type %in% fw_vintage_order &
        is.finite(raw$se), , drop = FALSE
    ]
    key <- do.call(paste, c(raw[keys], sep = "|"))
    if (anyDuplicated(key)) {
      stop("Duplicate country/quarter/vintage in ", model, " evaluation; ",
           "inspect its df_rt_all and GDP join before computing win shares.")
    }
    counts[[model]] <- nrow(raw)
    keep <- raw[c(keys, "se")]
    names(keep)[ncol(keep)] <- paste0("se_", model)
    if (model == "matrix") {
      p <- trimws(as.character(raw$period))
      # Preserve the original evaluation-period classification from the
      # quarterly target (do not classify individual M1/M2/M3 months).
      p[tolower(p) %in% c("pre-covid", "pre-covid period")] <- "Pre-COVID"
      p[tolower(p) %in% c("covid", "covid period")] <- "COVID period"
      p[tolower(p) %in% c("post-covid", "post-covid period")] <- "Post-COVID"
      keep$period <- p
      assembled[[model]] <- keep
    } else {
      assembled[[model]] <- keep
    }
  }
  joined <- Reduce(function(a, b) dplyr::inner_join(a, b, by = keys), assembled)
  joined <- joined[joined$period %in% fw_period_order, , drop = FALSE]
  if (!nrow(joined))
    stop("No forecasts evaluable across all four models; inspect summary dates and keys.")
  joined$country <- factor(joined$country, levels = countries)
  joined$type <- factor(joined$type, levels = fw_vintage_order)
  joined$period <- factor(joined$period, levels = fw_period_order)
  joined <- dplyr::arrange(joined, period, country, type, quarter_id)
  for (k in c("vecp", "vecc", "dfm")) {
    joined[[paste0("win_", k)]] <- as.integer(joined$se_matrix < joined[[paste0("se_", k)]])
  }
  if (anyNA(joined[, c("win_vecp", "win_vecc", "win_dfm")]))
    stop("Unexpected missing win indicators after evaluation matching.")
  attr(joined, "valid_forecasts_by_model") <- counts
  joined
}

fw_country_rmsfe <- function(evaluation) {
  # With the same observations for every model, the country-level RMSFE
  # comparison is equivalent to the sign of the mean squared-error difference.
  evaluation |>
    dplyr::group_by(period, country, type) |>
    dplyr::summarise(
      n = dplyr::n(),
      rmsfe_matrix = sqrt(mean(se_matrix)),
      rmsfe_vecp = sqrt(mean(se_vecp)),
      rmsfe_vecc = sqrt(mean(se_vecc)),
      rmsfe_dfm = sqrt(mean(se_dfm)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      win_vecp = as.integer(rmsfe_matrix < rmsfe_vecp),
      win_vecc = as.integer(rmsfe_matrix < rmsfe_vecc),
      win_dfm = as.integer(rmsfe_matrix < rmsfe_dfm)
    ) |>
    dplyr::arrange(period, country, type)
}

fw_country_rmsfe_share <- function(country_rmsfe) {
  country_rmsfe |>
    dplyr::group_by(period, type) |>
    dplyr::summarise(
      countries = dplyr::n(),
      vecp = 100 * mean(win_vecp),
      vecc = 100 * mean(win_vecc),
      dfm = 100 * mean(win_dfm),
      .groups = "drop"
    ) |>
    dplyr::arrange(period, type)
}

fw_forecast_summary <- function(evaluation) {
  one_group <- function(df, label) {
    data.frame(
      Group = label,
      Total = nrow(df),
      vecp_wins = sum(df$win_vecp),
      vecp_share = 100 * mean(df$win_vecp),
      vecc_wins = sum(df$win_vecc),
      vecc_share = 100 * mean(df$win_vecc),
      dfm_wins = sum(df$win_dfm),
      dfm_share = 100 * mean(df$win_dfm),
      stringsAsFactors = FALSE
    )
  }
  out <- list(one_group(evaluation, "Overall"))
  for (p in fw_period_order) {
    d <- evaluation[as.character(evaluation$period) == p, , drop = FALSE]
    out[[length(out) + 1L]] <- one_group(d, p)
  }
  for (v in fw_vintage_order) {
    d <- evaluation[as.character(evaluation$type) == v, , drop = FALSE]
    out[[length(out) + 1L]] <- one_group(d, v)
  }
  answer <- dplyr::bind_rows(out)
  if (sum(answer$Total[2:4]) != answer$Total[1] ||
      sum(answer$Total[5:7]) != answer$Total[1]) {
    stop("Internal inconsistency in forecast-level denominators.")
  }
  answer
}

fw_forecast_country <- function(evaluation) {
  evaluation |>
    dplyr::group_by(period, country, type) |>
    dplyr::summarise(
      n = dplyr::n(),
      vecp = 100 * mean(win_vecp),
      vecc = 100 * mean(win_vecc),
      dfm = 100 * mean(win_dfm),
      .groups = "drop"
    ) |>
    dplyr::arrange(period, country, type)
}

# ------------------------------------------------------------------------------
# LaTeX formatting: no old Small/Medium/Large and no screening terminology.
# Bold values exceed 50% (strict majority). Ties count as no win.
# ------------------------------------------------------------------------------
fw_pct <- function(x, bold = TRUE) {
  if (length(x) != 1L || !is.finite(x)) return("--")
  s <- sprintf("%.1f", x)
  if (isTRUE(bold) && x > 50) paste0("\\textbf{", s, "}") else s
}
fw_integer <- function(x) {
  if (length(x) != 1L || !is.finite(x)) "--" else as.character(as.integer(x))
}
fw_tex_row <- function(items) paste0(paste(items, collapse = " & "), " \\\\")
fw_period_title <- function(p) {
  switch(as.character(p), "Pre-COVID" = "Pre-Covid",
         "COVID period" = "Covid", "Post-COVID" = "Post-Covid")
}
fw_model_header <- function() {
  c(
    " & \\multicolumn{3}{c}{\\textbf{VEC-P}} & \\multicolumn{3}{c}{\\textbf{VEC-C}} & \\multicolumn{3}{c}{\\textbf{DFM}} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
    "\\rowcolor{topgray}",
    fw_tex_row(c("\\textbf{Period}", rep(c("\\textbf{M1}", "\\textbf{M2}", "\\textbf{M3}"), 3)))
  )
}
fw_table_start <- function(caption, label, size = "\\footnotesize", stretch = "1.05", tabsep = "4pt") {
  c(
    "\\begin{table}[!htbp]", "\\centering", "\\begingroup",
    size, paste0("\\renewcommand{\\arraystretch}{", stretch, "}"),
    paste0("\\setlength{\\tabcolsep}{", tabsep, "}"),
    "\\definecolor{topgray}{gray}{0.92}",
    "\\definecolor{regimegray}{gray}{0.965}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}")
  )
}
fw_table_end <- function(note) {
  c(
    "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.10cm}",
    paste0("\\parbox{0.95\\linewidth}{\\footnotesize\\textit{Notes:} ", note, "}"),
    "\\endgroup", "\\end{table}"
  )
}

fw_latex_country_rmsfe <- function(shares) {
  hdr <- fw_model_header()
  body <- character(0)
  for (p in fw_period_order) {
    x <- shares[as.character(shares$period) == p, , drop = FALSE]
    pick <- function(model, v) {
      y <- x[as.character(x$type) == v, , drop = FALSE]
      if (!nrow(y)) return("--")
      fw_pct(y[[model]][1])
    }
    body <- c(body, fw_tex_row(c(
      fw_period_title(p),
      vapply(c("vecp", "vecc", "dfm"), function(model)
        vapply(fw_vintage_order, function(v) pick(model, v), character(1)),
        character(3)) |> as.vector()
    )))
  }
  c(fw_table_start(
      "Country-level RMSFE win shares of Matrix MF--TPRF",
      "tab:country_rmsfe_win_shares", size = "\\small", tabsep = "5.4pt"
    ),
    "\\begin{tabular}{l@{\\hspace{0.40cm}}ccc@{\\hspace{0.30cm}}ccc@{\\hspace{0.30cm}}ccc}",
    "\\toprule", hdr, "\\midrule", body,
    fw_table_end(paste0(
      "Entries are percentages of countries with lower Matrix MF--TPRF RMSFE ",
      "than the benchmark in each period and vintage (common evaluable forecasts). ",
      "Bold indicates a share above 50\\%."
    ))
  )
}

fw_latex_forecast_summary <- function(x) {
  row <- function(r) fw_tex_row(c(
    r$Group, fw_integer(r$Total),
    fw_integer(r$vecp_wins), fw_pct(r$vecp_share),
    fw_integer(r$vecc_wins), fw_pct(r$vecc_share),
    fw_integer(r$dfm_wins), fw_pct(r$dfm_share)
  ))
  c(
    fw_table_start("Forecast-level wins of Matrix MF--TPRF against the benchmarks",
                   "tab:forecast_wins_summary", size = "\\small", tabsep = "5pt"),
    "\\begin{tabular}{l c@{\\hspace{0.30cm}}cc@{\\hspace{0.25cm}}cc@{\\hspace{0.25cm}}cc}",
    "\\toprule",
    " & & \\multicolumn{2}{c}{\\textbf{VEC-P}} & \\multicolumn{2}{c}{\\textbf{VEC-C}} & \\multicolumn{2}{c}{\\textbf{DFM}} \\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
    "\\rowcolor{topgray}",
    fw_tex_row(c("\\textbf{Group}", "\\textbf{Total}",
                 rep(c("\\textbf{Wins}", "\\textbf{Share}"), 3))),
    "\\midrule", row(x[1, ]),
    "\\addlinespace[0.25em]",
    "\\rowcolor{regimegray}",
    "\\multicolumn{8}{@{}l}{\\textbf{By evaluation period}} \\\\",
    vapply(seq.int(2L, 4L), function(i) row(x[i, ]), character(1)),
    "\\addlinespace[0.25em]",
    "\\rowcolor{regimegray}",
    "\\multicolumn{8}{@{}l}{\\textbf{By nowcast vintage}} \\\\",
    vapply(seq.int(5L, 7L), function(i) row(x[i, ]), character(1)),
    fw_table_end(paste0(
      "Wins count matched forecasts with strictly lower squared error for Matrix MF--TPRF. ",
      "Shares are percentages of the corresponding total; bold indicates a share above 50\\%."
    ))
  )
}

fw_latex_forecast_country <- function(x, countries = fw_country_order) {
  body <- character(0)
  for (p in fw_period_order) {
    body <- c(body, "\\addlinespace[0.18em]",
              "\\rowcolor{regimegray}",
              paste0("\\multicolumn{10}{@{}l}{\\textbf{", fw_period_title(p), "}} \\\\"),
              "\\cmidrule(lr){1-10}")
    for (cc in countries) {
      y <- x[as.character(x$period) == p & as.character(x$country) == cc,
             , drop = FALSE]
      pick <- function(model, v) {
        z <- y[as.character(y$type) == v, , drop = FALSE]
        if (!nrow(z)) return("--")
        fw_pct(z[[model]][1])
      }
      values <- unlist(lapply(c("vecp", "vecc", "dfm"), function(model) {
        vapply(fw_vintage_order, function(v) pick(model, v), character(1))
      }), use.names = FALSE)
      body <- c(body, fw_tex_row(c(cc, values)))
    }
  }
  header <- fw_model_header()
  header[length(header)] <- sub("Period", "Country", header[length(header)], fixed = TRUE)
  c(
    fw_table_start(
      "Forecast-level win shares of Matrix MF--TPRF by country, period and vintage",
      "tab:forecast_wins_country_period_month", size = "\\footnotesize",
      stretch = "0.98", tabsep = "3.7pt"
    ),
    "\\begin{tabular}{l@{\\hspace{0.30cm}}ccc@{\\hspace{0.32cm}}ccc@{\\hspace{0.32cm}}ccc}",
    "\\toprule", header, "\\midrule", body,
    fw_table_end(paste0(
      "Entries are percentages of matched forecasts with lower Matrix MF--TPRF squared error. ",
      "Bold indicates a share above 50\\%."
    ))
  )
}

fw_save_all <- function(evaluation,
                        output_dir = "results/paper/appendix/tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  country_rmsfe <- fw_country_rmsfe(evaluation)
  period_shares <- fw_country_rmsfe_share(country_rmsfe)
  forecast_summary <- fw_forecast_summary(evaluation)
  forecast_country <- fw_forecast_country(evaluation)

  files <- c(
    country_rmsfe = "table_country_rmsfe_wins.tex",
    forecast_summary = "table_forecast_win_summary.tex",
    forecast_country = "table_forecast_win_country.tex"
  )
  tex <- list(
    fw_latex_country_rmsfe(period_shares),
    fw_latex_forecast_summary(forecast_summary),
    fw_latex_forecast_country(forecast_country)
  )
  for (i in seq_along(files)) {
    writeLines(tex[[i]], file.path(output_dir, files[[i]]), useBytes = TRUE)
  }
  utils::write.csv(country_rmsfe, file.path(output_dir, "country_rmsfe_wins.csv"),
                   row.names = FALSE)
  utils::write.csv(period_shares, file.path(output_dir, "country_rmsfe_win_shares.csv"),
                   row.names = FALSE)
  utils::write.csv(forecast_summary, file.path(output_dir, "forecast_win_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(forecast_country, file.path(output_dir, "forecast_win_country.csv"),
                   row.names = FALSE)
  utils::write.csv(evaluation, file.path(output_dir, "forecast_win_common_errors.csv"),
                   row.names = FALSE)
  message("Saved three LaTeX tables and five audit CSVs to ", output_dir)
  invisible(list(country_rmsfe = country_rmsfe,
                 country_rmsfe_share = period_shares,
                 forecast_summary = forecast_summary,
                 forecast_country = forecast_country,
                 filenames = file.path(output_dir, unname(files))))
}
