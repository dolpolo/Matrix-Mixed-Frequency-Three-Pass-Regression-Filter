# ==============================================================================
# Input validation | Fixed information set
# Author: Davide Delfino
# Institution: Aarhus University
# Paper: Matrix Three-Pass Regression Filter: A Nowcasting Framework for Mixed-Frequency Multidimensional Data
# Last updated: 2026-11-24
# Input: fixed predictor CSV; raw and transformed EA-MD-QD workbooks
# Output: input checks in the console
# ==============================================================================

file <- "data/information_set/selected_variables_50.csv"
if (!file.exists(file)) stop("Missing fixed predictor CSV: ", file)
a <- utils::read.csv(file, stringsAsFactors = FALSE)
if (!all(c("base_name", "frequency") %in% names(a)) || nrow(a) != 50L ||
    anyDuplicated(a$base_name) || sum(a$frequency == "M") != 35L ||
    sum(a$frequency == "Q") != 15L) {
  stop("Expected 50 distinct predictors (35 monthly and 15 quarterly).")
}
cc <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
required <- c(
  file.path("data/raw", paste0(cc, "data.xlsx")),
  file.path("data/transformed", paste0(cc, "dataM_NA_TR2.xlsx")),
  file.path("data/transformed", paste0(cc, "dataQ_NA_TR2.xlsx"))
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing data inputs:\n", paste(missing, collapse = "\n"))
cat("Input check passed: 50 fixed predictors; 27 required workbooks.\n")
