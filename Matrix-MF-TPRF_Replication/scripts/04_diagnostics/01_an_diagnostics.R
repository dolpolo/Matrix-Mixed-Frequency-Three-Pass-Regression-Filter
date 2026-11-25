# Mandatory AN diagnostic output (not necessarily a paper figure).
# Four distinct plots: row/column histogram and normal Q-Q, each PNG+PDF.
scenario <- "results/monte_carlo/an_pass1/T200_p1100_p2100"
plots <- c("histogram_standardized_error", "qqplot_standardized_error")
axes <- c("row_loadings", "column_loadings")
expected <- unlist(lapply(axes, function(a) unlist(lapply(plots, function(p)
  file.path(scenario, a, paste0(p, c(".png", ".pdf")))))))
if (!all(file.exists(expected))) {
  message("AN figures missing; running the original Pass 1 Monte Carlo experiment.")
  source("scripts/02_monte_carlo/03_an_pass1.R", local=new.env(parent=globalenv()))
}
if (!all(file.exists(expected)))
  stop("AN diagnostics missing after simulation: ", paste(expected[!file.exists(expected)],collapse=", "))
for (ax in axes) {
  dst_dir <- file.path("results/diagnostics/asymptotic_normality", ax)
  dir.create(dst_dir, recursive=TRUE, showWarnings=FALSE)
  for (p in plots) for (ext in c("png", "pdf")) {
    src <- file.path(scenario, ax, paste0(p,".",ext))
    dst <- file.path(dst_dir, basename(src))
    if (!file.copy(src,dst,overwrite=TRUE)) stop("Cannot publish AN diagnostics: ", src)
  }
}
message("Mandatory AN histograms and QQ plots validated (row and column).")
