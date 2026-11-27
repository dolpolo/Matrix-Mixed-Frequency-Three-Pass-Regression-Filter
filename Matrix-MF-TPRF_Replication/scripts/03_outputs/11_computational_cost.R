# ==============================================================================
# MATRIX THREE-PASS REGRESSION FILTER:
# A NOWCASTING FRAMEWORK FOR MIXED-FREQUENCY MULTIDIMENSIONAL DATA
# ==============================================================================
# Authors:      Matteo Barigozzi, Davide Delfino, Massimiliano Marcellino
# Institutions: University of Bologna, Aarhus University, Bocconi University
# Contact:      davide.delfino@econ.au.dk
# File:         scripts/03_outputs/11_computational_cost.R
# Input:        results/audit/computational_cost/timing_records.rds
# Output:       Appendix LaTeX table and audit CSV/system information.
# ==============================================================================

source("R/common/computational_cost.R")

if (!file.exists(computational_records_file)) {
  stop("No recorded timings. Run: Rscript run_all.R benchmark")
}

records <- readRDS(computational_records_file)
if (!identical(names(records), computational_components)) {
  stop("Invalid computational timing records.")
}

# 1. Construct the comparison --------------------------------------------------

elapsed <- vapply(records, function(x) x$elapsed_seconds, numeric(1))
workers <- vapply(records, function(x) x$workers, integer(1))
budgets <- vapply(records, function(x) x$worker_budget, integer(1))
if (any(!is.finite(elapsed) | elapsed <= 0) ||
    length(unique(budgets)) != 1L || workers[["matrix"]] != 1L ||
    any(workers[c("vecp", "dfm")] != budgets[c("vecp", "dfm")])) {
  stop("Invalid or inconsistent computational timings.")
}

comparison <- data.frame(
  model = c("Matrix MF-TPRF", "VEC-P", "DFM"),
  workers = unname(workers),
  elapsed_seconds = unname(elapsed),
  elapsed_minutes = unname(elapsed / 60),
  relative_time = unname(elapsed / elapsed[["matrix"]]),
  row.names = NULL
)

# 2. Save numerical results and system information -----------------------------

audit_dir <- dirname(computational_records_file)
utils::write.csv(comparison, file.path(audit_dir, "computational_times.csv"),
                 row.names = FALSE)
writeLines(capture.output(str(records, max.level = 3)),
           file.path(audit_dir, "system_information.txt"))

# 3. Generate the appendix table ------------------------------------------------

out_dir <- "results/paper/appendix/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file_table <- file.path(out_dir, "table_computational_cost.tex")

rows <- sprintf("%s & %d & %.2f & %.2f %s",
                comparison$model, comparison$workers,
                comparison$elapsed_minutes, comparison$relative_time,
                "\\\\")

notes <- paste0(
  "Notes: Wall-clock minutes for one run of the empirical estimation pipeline, ",
  "including full-sample estimation, hyperparameter calibration, rolling nowcasts ",
  "and estimation-result serialization. Initial panel construction and common ",
  "data preparation, publication tables and figures, and the matrix fixed-rank ",
  "appendix exercise are excluded. Vector worker setup and country-specific ",
  "input preparation within the parallel pipeline are included. VEC-P and DFM ",
  "use ", budgets[[1L]], " parallel R workers; the matrix implementation uses ",
  "one R process. Numerical-library thread counts are requested to be one per ",
  "process. The matrix model also produces the EA aggregate nowcast, whereas ",
  "VEC-P and DFM produce national nowcasts; DFM uses country-specific rather ",
  "than pooled predictors. The comparison concerns the implemented workflows, ",
  "not an isolated measure of the effect of parameterization. Relative times ",
  "are normalized by the Matrix MF-TPRF time."
)

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Computational costs of the nowcasting methods}",
  "\\label{tab:computational_cost}",
  "\\begin{tabular}{lrrr}",
  "\\toprule",
  paste("Model & R workers & Time (min.) & Relative time", "\\\\"),
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\medskip",
  "\\begin{minipage}{0.95\\linewidth}\\footnotesize",
  notes,
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(lines, file_table)
print(comparison, row.names = FALSE)
message("Saved computational-cost table: ", file_table)
