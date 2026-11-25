# Publish the loading-structure figure and the baseline row-Wald table.
# Scientific simulations/results retain their own separate scripts.
if (!file.exists("config/settings.R")) stop("Run from project root")
loading_dir <- "results/monte_carlo/loading_r2x2_T201_p20x50_rG0.95_rE0.9"
source_stem <- file.path(loading_dir, "01_matrix_vs_vector_2x2")
if (!file.exists(paste0(source_stem, ".pdf")) ||
    !file.exists(paste0(source_stem, ".png"))) {
  message("Loading-structure figure not found. Executing its existing simulation script.")
  source("scripts/02_monte_carlo/02_loading_maps.R", local = new.env(parent=globalenv()))
}
fig_dir <- "results/paper/main/figures"
tab_dir <- "results/paper/main/tables"
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(tab_dir, recursive=TRUE, showWarnings=FALSE)
for (ext in c("pdf", "png")) {
  src <- paste0(source_stem, ".", ext)
  dst <- file.path(fig_dir, paste0("figure_02_matrix_vs_vector_loading_maps.", ext))
  if (!file.exists(src)) stop("Missing required loading figure: ", src)
  if (!file.copy(src, dst, overwrite=TRUE)) stop("Unable to copy ", src)
}
src <- "results/empirical/matrix/tables/row_wald_scalar_fixed50_fixed.tex"
dst <- file.path(tab_dir, "table_row_wald_single_restrictions.tex")
if (file.exists(src)) {
  if (!file.copy(src, dst, overwrite=TRUE)) stop("Unable to copy ", src)
  message("Published loading figure and baseline Wald table.")
} else {
  message("Published loading figure; Wald table is unavailable until Matrix estimation.")
}
