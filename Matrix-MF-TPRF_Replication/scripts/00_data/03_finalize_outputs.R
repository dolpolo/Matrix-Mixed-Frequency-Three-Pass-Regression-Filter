# Move technical CSV audit files out of results/paper without losing information.
# Preserve exactly the approved publication files (PDF, PNG, .tex).
manifest <- utils::read.csv("config/output_manifest.csv", stringsAsFactors=FALSE)
allowed <- manifest$path[startsWith(manifest$path, "results/paper/")]
files <- list.files("results/paper", recursive=TRUE, full.names=TRUE)
files <- files[file.info(files)$isdir == FALSE]
extra <- setdiff(files, allowed)
if (length(extra)) {
  for (src in extra) {
    rel <- substring(src, nchar("results/paper/") + 1L)
    dst <- file.path("results/audit/paper_intermediates", rel)
    dir.create(dirname(dst), recursive=TRUE, showWarnings=FALSE)
    if (!file.copy(src, dst, overwrite=TRUE)) stop("Cannot archive intermediate: ",src)
    if (!file.remove(src)) stop("Cannot remove extra paper file: ",src)
    message("Archived non-publication output: ",src)
  }
}
