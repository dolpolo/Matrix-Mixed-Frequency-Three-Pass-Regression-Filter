# Author: Davide Delfino
# Input: data/information_set/selected_variables_50.csv
# Output: fixed national predictor lists (requested names, with coverage report).

load_fixed_selection <- function(countries, path_raw = NULL, path_adj = NULL,
                                 file = "data/information_set/selected_variables_50.csv") {
  if (!file.exists(file)) {
    stop("Missing fixed information set: ", file,
         ". Supply the fixed CSV before running the models.")
  }
  selected <- utils::read.csv(file, stringsAsFactors = FALSE)
  required <- c("base_name", "frequency")
  if (!all(required %in% names(selected))) stop("Invalid information-set columns.")
  if (nrow(selected) != 50L || anyDuplicated(selected$base_name) ||
      !all(selected$frequency %in% c("M", "Q")) || anyNA(selected$base_name) ||
      sum(selected$frequency == "M") != 35L ||
      sum(selected$frequency == "Q") != 15L) {
    stop("Information set must contain 35 monthly and 15 quarterly distinct predictors.")
  }
  if ("GDP" %in% toupper(selected$base_name)) stop("GDP target cannot be a predictor.")

  m <- selected$base_name[selected$frequency == "M"]
  q <- selected$base_name[selected$frequency == "Q"]
  selections <- list(
    m = setNames(lapply(countries, function(cc) paste0(m, "_", cc)), countries),
    q = setNames(lapply(countries, function(cc) paste0(c("GDP", q), "_", cc)), countries)
  )

  if (!is.null(path_adj)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl first.")
    missing_by_country <- list()
    for (cc in countries) {
      dm <- names(readxl::read_excel(file.path(path_adj, paste0(cc, "dataM_NA_TR2.xlsx")),
                                       n_max = 0))
      dq <- names(readxl::read_excel(file.path(path_adj, paste0(cc, "dataQ_NA_TR2.xlsx")),
                                       n_max = 0))
      missing <- c(setdiff(selections$m[[cc]], dm), setdiff(selections$q[[cc]], dq))
      if (length(missing)) missing_by_country[[cc]] <- missing
    }
    if (length(missing_by_country)) {
      msg <- paste(vapply(names(missing_by_country), function(cc) {
        paste0(cc, ": ", paste(missing_by_country[[cc]], collapse = ", "))
      }, character(1)), collapse = " | ")
      warning("Some union predictors are not present in national/EA data: ", msg,
              ". They must remain missing in the matrix; no substitute predictors are selected.",
              call. = FALSE)
    }
  }
  selections
}
