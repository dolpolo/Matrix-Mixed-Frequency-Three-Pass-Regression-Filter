# ==============================================================================
# Country nowcasts | national panels for the appendix
# Input: saved Matrix MF-TPRF, VEC-P, VEC-C and DFM results.
# Output: 8 separate country figures (each with M1, M2, M3), plus the
#         UNCHANGED post-Covid 8-country within-quarter figure.
# No re-estimation. Called by scripts/03_outputs/09_country_nowcasts_main.R.
# ==============================================================================

# Shared model styles: also used by the unchanged post-Covid figure.
nowcast_model_names <- c(
  matrix = "Matrix MF-TPRF", vecp = "VEC-P", vecc = "VEC-C", dfm = "DFM"
)
nowcast_model_colors <- c(
  "Matrix MF-TPRF" = "#1B6CA8", "VEC-P" = "#55A868",
  "VEC-C" = "#C44E52", "DFM" = "grey35", "GDP" = "black"
)
nowcast_model_linetypes <- c(
  "Matrix MF-TPRF" = "solid", "VEC-P" = "11",
  "VEC-C" = "33", "DFM" = "44", "GDP" = "solid"
)

# ----------------------------------------------------------------------------
# Full evaluation: one country per figure, with M1 / M2 / M3 side by side.
# ----------------------------------------------------------------------------
load_country_nowcasts_summaries <- function(size = "fixed50", selection = "fixed") {
  # The two file-locator helpers are defined in R/matrix/matrix.mf.tprf.utils.R
  # and are sourced by the main script.
  if (!exists("find_result_file", mode = "function") ||
      !exists("find_dfm_summary_file", mode = "function")) {
    stop("Source R/matrix/matrix.mf.tprf.utils.R before loading summaries.")
  }

  load_one <- function(path, model = NULL, is_dfm = FALSE) {
    file <- if (is_dfm) {
      find_dfm_summary_file(path = path, Size = size, sel = selection)
    } else {
      find_result_file(path = path, model = model, stage = "summary",
                       Size = size, sel = selection)
    }
    if (length(file) != 1L || !file.exists(file)) {
      stop("Summary RDS not found in: ", path)
    }
    result <- readRDS(file)
    if (!is.data.frame(result$df_rt_all)) {
      stop("Saved summary has no df_rt_all data frame: ", file)
    }
    result
  }

  out <- list(
    matrix = load_one("results/empirical/matrix", "matrix_scalar"),
    vecp   = load_one("results/empirical/pooled", "vectensor"),
    vecc   = load_one("results/empirical/vector", "vector"),
    dfm    = load_one("results/empirical/dfm/q0", is_dfm = TRUE)
  )
  if (!identical(as.character(out$matrix$proxy_mode), "scalar")) {
    stop("Expected the baseline scalar-proxy Matrix MF-TPRF summary.")
  }
  if (!is.data.frame(out$matrix$df_yq_eval_all)) {
    stop("Matrix summary does not contain df_yq_eval_all (observed GDP).")
  }
  out
}

prepare_country_nowcast_data <- function(summaries, countries,
                                               start_date, end_date) {
  stopifnot(length(countries) == 1L, !anyDuplicated(countries))
  models <- nowcast_model_names
  vintages <- c("M1", "M2", "M3")
  needed_rt <- c("date", "country", "type", "nowcast")
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  model_data <- lapply(names(models), function(k) {
    x <- summaries[[k]]$df_rt_all
    if (!all(needed_rt %in% names(x))) {
      stop("Missing nowcast columns in summary: ", k, ".")
    }
    x <- dplyr::filter(x, as.character(country) %in% countries)
    x <- dplyr::transmute(
      x,
      date = as.Date(date),
      country = as.character(country),
      type = as.character(type),
      value = as.numeric(nowcast),
      model = unname(models[[k]])
    )
    x <- dplyr::filter(x, date >= start_date, date <= end_date,
                       type %in% vintages, is.finite(value))
    x
  })
  model_df <- dplyr::bind_rows(model_data)

  expected <- expand.grid(country = countries, type = vintages,
                          model = unname(models), stringsAsFactors = FALSE)
  seen <- unique(model_df[c("country", "type", "model")])
  absent <- dplyr::anti_join(expected, seen,
                             by = c("country", "type", "model"))
  if (nrow(absent)) {
    stop("No saved nowcasts for these country/vintage/model combinations: ",
         paste(apply(absent, 1, paste, collapse = "/"), collapse = ", "))
  }

  gdp <- summaries$matrix$df_yq_eval_all
  if (!all(c("date", "country", "GDP") %in% names(gdp))) {
    stop("The matrix GDP data must have date, country and GDP columns.")
  }
  gdp <- dplyr::transmute(gdp, date = as.Date(date),
                           country = as.character(country),
                           value = as.numeric(GDP))
  gdp <- dplyr::filter(gdp, country %in% countries,
                       date >= start_date, date <= end_date,
                       is.finite(value))
  gdp <- dplyr::distinct(gdp, country, date, .keep_all = TRUE)
  missing_gdp <- setdiff(countries, unique(gdp$country))
  if (length(missing_gdp)) {
    stop("No observed GDP in matrix summary for: ",
         paste(missing_gdp, collapse = ", "))
  }

  # Identical observed quarterly GDP is displayed in M1, M2, M3,
  # matching the paper's original per-country graphics.
  gdp <- tidyr::crossing(gdp, type = vintages)
  gdp$model <- "GDP"

  d <- dplyr::bind_rows(model_df, gdp)
  d$country <- factor(d$country, levels = countries)
  d$type    <- factor(d$type, levels = vintages)
  d$model   <- factor(d$model,
                      levels = c(unname(models), "GDP"))
  dplyr::arrange(d, country, type, model, date)
}

# Each single-country figure matches the original three-panel layout.
# Annual labels are vertical so 2026 cannot overlap 2017 in the next facet.
# The three vintages use a COMMON y scale, independently fitted by country.
make_country_nowcast <- function(
    summaries, country,
    start_date = as.Date("2017-01-01"),
    end_date = as.Date("2026-02-01"),
    covid_start = as.Date("2020-03-01"),
    covid_end = as.Date("2021-07-01"),
    base_size = 11,
    model_linewidth = 0.78,
    gdp_linewidth = 1.10
) {
  country_names <- c(
    DE = "Germany", FR = "France", IT = "Italy", ES = "Spain",
    NL = "Netherlands", BE = "Belgium", AT = "Austria", PT = "Portugal"
  )
  if (length(country) != 1L || !country %in% names(country_names)) {
    stop("country must be one of: ", paste(names(country_names), collapse = ", "))
  }
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  covid_start <- as.Date(covid_start)
  covid_end <- as.Date(covid_end)

  d <- prepare_country_nowcast_data(summaries, country, start_date, end_date)
  display_models <- c(unname(nowcast_model_names), "GDP")
  covid_band <- data.frame(
    xmin = covid_start,
    xmax = seq(covid_end, by = "1 month", length.out = 2L)[2L]
  )
  year_breaks <- seq(
    as.Date(sprintf("%s-01-01", format(start_date, "%Y"))),
    as.Date(sprintf("%s-01-01", format(end_date, "%Y"))),
    by = "1 year"
  )

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = covid_band,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "grey83", alpha = 0.38, colour = NA
    ) +
    ggplot2::geom_hline(
      yintercept = 0, colour = "grey78", linewidth = 0.28
    ) +
    ggplot2::geom_line(
      data = d[as.character(d$model) != "GDP", , drop = FALSE],
      ggplot2::aes(x = date, y = value, colour = model,
                   linetype = model, group = model),
      linewidth = model_linewidth, lineend = "round", na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = d[as.character(d$model) == "GDP", , drop = FALSE],
      ggplot2::aes(x = date, y = value, colour = model,
                   linetype = model, group = model),
      linewidth = gdp_linewidth, lineend = "round", na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~type, nrow = 1L, drop = FALSE) +
    ggplot2::scale_x_date(
      limits = c(start_date, end_date),
      breaks = year_breaks, date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.012, 0.012))
    ) +
    ggplot2::scale_y_continuous(
      labels = function(y) formatC(y, format = "f", digits = 2),
      expand = ggplot2::expansion(mult = c(0.065, 0.065))
    ) +
    ggplot2::scale_colour_manual(
      values = nowcast_model_colors, breaks = display_models,
      limits = display_models, name = NULL, drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = nowcast_model_linetypes, breaks = display_models,
      limits = display_models, name = NULL, drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1, nrow = 1, override.aes = list(linewidth = 1.2)
      ),
      linetype = ggplot2::guide_legend(order = 1, nrow = 1)
    ) +
    ggplot2::labs(
      title = paste("Expanding pseudo-real-time nowcasts for",
                    unname(country_names[[country]])),
      x = NULL, y = "GDP growth"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = 14,
        margin = ggplot2::margin(b = 6)
      ),
      plot.margin = ggplot2::margin(t = 6, r = 7, b = 4, l = 3),
      panel.grid.major.x = ggplot2::element_line(
        colour = "grey90", linewidth = 0.25
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey88", linewidth = 0.30
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey76", fill = NA, linewidth = 0.35
      ),
      panel.spacing.x = grid::unit(0.14, "cm"),
      strip.background = ggplot2::element_rect(
        fill = "grey89", colour = "grey69", linewidth = 0.4
      ),
      strip.text = ggplot2::element_text(
        face = "bold", size = 11, colour = "black",
        margin = ggplot2::margin(t = 4, b = 4)
      ),
      axis.text.x = ggplot2::element_text(
        size = 8.2, colour = "grey20",
        angle = 90, hjust = 1, vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 8.3, colour = "grey25"
      ),
      axis.title.y = ggplot2::element_text(
        size = 10, margin = ggplot2::margin(r = 5)
      ),
      axis.ticks.x = ggplot2::element_line(
        colour = "grey45", linewidth = 0.28
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.text = ggplot2::element_text(size = 9.5),
      legend.key.width = grid::unit(0.78, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.06, "cm"),
      legend.box.margin = ggplot2::margin(t = 5)
    )
}

# Export with the exact country filenames referenced in the user's LaTeX.
save_country_nowcast <- function(
    plot, country,
    output_dir = "results/paper/appendix/figures",
    width = 12.8, height = 4.85, dpi = 450
) {
  if (length(country) != 1L) stop("Supply one country code.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("figure_07_country_", country)
  paths <- file.path(output_dir, paste0(stem, c(".pdf", ".png")))
  for (path in paths) {
    ggplot2::ggsave(
      filename = path, plot = plot,
      width = width, height = height, units = "in",
      dpi = dpi, bg = "white", limitsize = FALSE
    )
    message("Saved: ", path)
  }
  invisible(paths)
}

# ----------------------------------------------------------------------------
# Post-Covid: monthly M1-M2-M3 nowcast paths and quarterly GDP at M3 only.
# ----------------------------------------------------------------------------
postcovid_quarter_m3_date <- function(dates) {
  d <- as.Date(as.character(dates))
  yy <- as.integer(format(d, "%Y"))
  mm <- as.integer(format(d, "%m"))
  last_month <- 3L * ((mm - 1L) %/% 3L + 1L)
  as.Date(sprintf("%04d-%02d-01", yy, last_month))
}

prepare_postcovid_country_data <- function(
  summaries,
  countries = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  start_date = as.Date("2021-08-01"),
  end_date = as.Date("2026-02-01")
) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  if (is.na(start_date) || is.na(end_date) || start_date > end_date) {
    stop("Invalid post-Covid evaluation dates.")
  }
  if (length(countries) != 8L || anyDuplicated(countries)) {
    stop("Supply exactly eight different countries.")
  }

  needed_rt <- c("date", "country", "type", "nowcast")
  model_data <- lapply(names(nowcast_model_names), function(key) {
    raw <- summaries[[key]]$df_rt_all
    missing_cols <- setdiff(needed_rt, names(raw))
    if (length(missing_cols)) {
      stop("Missing ", paste(missing_cols, collapse = ", "),
           " in ", key, " nowcasts.")
    }
    # Align each within-quarter forecast to its calendar month, regardless
    # of whether the saved date is month-start or month-end.
    saved_date <- as.Date(as.character(raw$date))
    calendar_month <- as.Date(format(saved_date, "%Y-%m-01"), format = "%Y-%m-%d")
    x <- data.frame(
      date = calendar_month,
      country = as.character(raw$country),
      type = as.character(raw$type),
      value = as.numeric(raw$nowcast),
      series = unname(nowcast_model_names[[key]]),
      stringsAsFactors = FALSE
    )
    x <- x[x$country %in% countries &
             !is.na(x$date) & x$date >= start_date & x$date <= end_date &
             x$type %in% c("M1", "M2", "M3"), , drop = FALSE]
    x
  })
  nowcasts <- do.call(rbind, model_data)
  rownames(nowcasts) <- NULL

  # Every available nowcast belongs to one model, country and monthly date.
  # A repeated date would spuriously connect several vintages at the same x.
  key <- paste(nowcasts$series, nowcasts$country, nowcasts$date, sep = "|")
  if (anyDuplicated(key)) {
    stop("Duplicated (model, country, date) in saved nowcasts; ",
         "check the underlying df_rt_all before plotting.")
  }
  present <- unique(nowcasts[, c("country", "series"), drop = FALSE])
  expected <- expand.grid(country = countries,
                          series = unname(nowcast_model_names),
                          stringsAsFactors = FALSE)
  missing <- expected[!paste(expected$country, expected$series) %in%
                        paste(present$country, present$series), , drop = FALSE]
  if (nrow(missing)) {
    stop("No saved post-Covid nowcasts for: ",
         paste(paste(missing$country, missing$series), collapse = "; "))
  }

  # Retain missing monthly nowcasts as NA gaps rather than interpolate across
  # unpublished/unevaluable months. Within-quarter M1--M2--M3 are linked.
  months <- seq(as.Date(format(start_date, "%Y-%m-01")),
                as.Date(format(end_date, "%Y-%m-01")), by = "1 month")
  grid <- expand.grid(
    date = months,
    country = countries,
    series = unname(nowcast_model_names),
    stringsAsFactors = FALSE
  )
  nowcasts <- merge(grid, nowcasts,
                    by = c("date", "country", "series"),
                    all.x = TRUE, sort = FALSE)

  # Preserve gaps when a monthly vintage is unavailable: do not draw a
  # line across intervening NA months. Each nonmissing run is a segment.
  nowcasts <- nowcasts[order(nowcasts$country, nowcasts$series, nowcasts$date), ]
  run_id <- interaction(nowcasts$country, nowcasts$series, drop = TRUE)
  nowcasts$segment <- ave(
    as.integer(is.na(nowcasts$value)), run_id,
    FUN = function(x) cumsum(x)
  )

  # The GDP summary is quarterly. Plot GDP at M3 only, with one point per
  # quarter, then connect quarterly GDP points (never fill M1/M2).
  raw_gdp <- summaries$matrix$df_yq_eval_all
  if (!all(c("date", "country", "GDP") %in% names(raw_gdp))) {
    stop("Matrix df_yq_eval_all must contain date, country and GDP.")
  }
  gdp <- data.frame(
    date = postcovid_quarter_m3_date(raw_gdp$date),
    country = as.character(raw_gdp$country),
    value = as.numeric(raw_gdp$GDP),
    type = "M3",
    series = "GDP",
    stringsAsFactors = FALSE
  )
  gdp <- gdp[!is.na(gdp$date) & gdp$date >= start_date &
               gdp$date <= end_date & gdp$country %in% countries &
               is.finite(gdp$value), , drop = FALSE]
  if (anyDuplicated(paste(gdp$country, gdp$date))) {
    stop("Repeated country/quarter observations in saved GDP series.")
  }
  missing_gdp <- setdiff(countries, unique(gdp$country))
  if (length(missing_gdp)) {
    stop("Missing quarterly GDP series for ", paste(missing_gdp, collapse = ", "))
  }

  country_order <- countries
  nowcasts$country <- factor(nowcasts$country, levels = country_order)
  gdp$country <- factor(gdp$country, levels = country_order)
  nowcasts$series <- factor(nowcasts$series,
                            levels = c(unname(nowcast_model_names), "GDP"))
  gdp$series <- factor(gdp$series,
                       levels = c(unname(nowcast_model_names), "GDP"))
  nowcasts <- nowcasts[order(nowcasts$country, nowcasts$series, nowcasts$date), ]
  gdp <- gdp[order(gdp$country, gdp$date), ]

  list(nowcasts = nowcasts, gdp = gdp,
       combined = rbind(nowcasts[, c("date", "country", "type", "value", "series")],
                        gdp[, c("date", "country", "type", "value", "series")]))
}

make_postcovid_country_plot <- function(
  data,
  start_date = as.Date("2021-08-01"),
  end_date = as.Date("2026-02-01"),
  free_y = TRUE,
  base_size = 11,
  show_monthly_points = TRUE,
  x_labels = c("annual", "quarterly")
) {
  x_labels <- match.arg(x_labels)
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  model_levels <- c(unname(nowcast_model_names), "GDP")

  # Annual labels are vertical. Quarterly boundaries appear as light gridlines.
  # Set x_labels = "quarterly" to label each quarter instead (optional).
  first_next_year <- as.Date(sprintf("%d-01-01", as.integer(format(start_date, "%Y")) + 1L))
  annual_ticks <- c(
    start_date,
    if (first_next_year <= end_date) seq(
      first_next_year,
      as.Date(sprintf("%d-01-01", as.integer(format(end_date, "%Y")))),
      by = "1 year"
    ) else as.Date(character())
  )
  quarterly_ticks <- seq(
    as.Date(sprintf("%s-01-01", format(start_date, "%Y"))),
    as.Date(sprintf("%s-01-01", as.integer(format(end_date, "%Y")) + 1L)),
    by = "3 months"
  )
  quarterly_ticks <- quarterly_ticks[
    quarterly_ticks >= start_date & quarterly_ticks <= end_date
  ]
  if (identical(x_labels, "quarterly")) {
    axis_ticks <- sort(unique(c(start_date, quarterly_ticks)))
    axis_labels <- function(x) {
      q <- (as.integer(format(x, "%m")) - 1L) %/% 3L + 1L
      show_year <- as.integer(format(x, "%m")) == 1L | x == start_date
      ifelse(show_year, paste0(format(x, "%Y"), " Q", q), paste0("Q", q))
    }
  } else {
    axis_ticks <- annual_ticks
    axis_labels <- function(x) format(x, "%Y")
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey77", linewidth = 0.40) +
    ggplot2::geom_line(
      data = data$nowcasts,
      ggplot2::aes(x = date, y = value, colour = series,
                   linetype = series, group = interaction(series, segment)),
      linewidth = 0.8, lineend = "round", na.rm = FALSE
    )

  if (isTRUE(show_monthly_points)) {
    p <- p + ggplot2::geom_point(
      data = subset(data$nowcasts, is.finite(value)),
      ggplot2::aes(x = date, y = value, colour = series),
      size = 1.1, alpha = 0.85, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_line(
      data = data$gdp,
      ggplot2::aes(x = date, y = value, colour = series,
                   linetype = series, group = series),
      linewidth = 1.72, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = data$gdp, ggplot2::aes(x = date, y = value),
      colour = "black", size = 1.20, show.legend = FALSE
    ) +
    ggplot2::facet_wrap(
      ~country, ncol = 2, scales = if (free_y) "free_y" else "fixed",
      drop = FALSE
    ) +
    ggplot2::scale_x_date(
      limits = c(start_date, end_date),
      breaks = axis_ticks, minor_breaks = quarterly_ticks,
      labels = axis_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) formatC(x, format = "f", digits = 2),
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::scale_colour_manual(
      values = nowcast_model_colors, breaks = model_levels,
      limits = model_levels, name = NULL, drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = nowcast_model_linetypes, breaks = model_levels,
      limits = model_levels, name = NULL, drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, nrow = 1,
                                     override.aes = list(linewidth = 1.35)),
      linetype = ggplot2::guide_legend(order = 1, nrow = 1)
    ) +
    ggplot2::labs(x = NULL, y = "GDP growth") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 4, r = 7, b = 4, l = 4),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.31),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_line(
        colour = "grey92", linewidth = 0.32
      ),
      panel.border = ggplot2::element_rect(
        colour = "grey72", fill = NA, linewidth = 0.42
      ),
      panel.spacing.x = grid::unit(0.45, "cm"),
      panel.spacing.y = grid::unit(0.40, "cm"),
      strip.background = ggplot2::element_rect(
        fill = "grey88", colour = "grey68", linewidth = 0.42
      ),
      strip.text = ggplot2::element_text(
        face = "bold", colour = "black", size = 11,
        margin = ggplot2::margin(t = 4, b = 4)
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.50,
        size = 9, colour = "grey20"
      ),
      axis.text.y = ggplot2::element_text(size = 8.5, colour = "grey22"),
      axis.title.y = ggplot2::element_text(size = 10,
                                          margin = ggplot2::margin(r = 5)),
      axis.ticks.x = ggplot2::element_line(colour = "grey42", linewidth = 0.32),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.text = ggplot2::element_text(size = 9.6),
      legend.key.width = grid::unit(0.80, "cm"),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.spacing.x = grid::unit(0.09, "cm"),
      legend.box.margin = ggplot2::margin(t = 8)
    )
  p
}

save_postcovid_country_plot <- function(
  plot, data,
  output_dir = "results/paper/appendix/figures",
  stem = "figure_08_postcovid_country_nowcasts_4x2",
  width = 8.1, height = 10.4, dpi = 450
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    path <- file.path(output_dir, paste0(stem, ".", ext))
    ggplot2::ggsave(filename = path, plot = plot,
                    width = width, height = height, units = "in",
                    dpi = dpi, bg = "white", limitsize = FALSE)
    message("Saved: ", path)
  }
  csv <- file.path(output_dir, paste0(stem, "_values.csv"))
  utils::write.csv(data$combined, file = csv, row.names = FALSE,
                   na = "")
  message("Saved: ", csv)
  invisible(c(file.path(output_dir, paste0(stem, c(".pdf", ".png"))), csv))
}
