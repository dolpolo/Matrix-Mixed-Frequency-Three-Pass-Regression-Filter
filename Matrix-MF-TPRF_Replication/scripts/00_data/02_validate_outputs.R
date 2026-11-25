# Validate paper+mandatory AN outputs against the explicit source manifest.
manifest <- utils::read.csv("config/output_manifest.csv", stringsAsFactors=FALSE)
required <- manifest$path[manifest$required]
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing publication/diagnostic outputs: ",
                          paste(missing, collapse="; "))
# Intermediate CSVs are allowed in results/audit, but not publication folders.
paper_dir <- "results/paper"
present <- list.files(paper_dir, full.names=TRUE, recursive=TRUE,
                      all.files=FALSE)
present <- present[file.info(present)$isdir == FALSE]
expected_paper <- manifest$path[grepl("^results/paper/",manifest$path)]
unexpected <- setdiff(present, expected_paper)
if (length(unexpected)) stop("Unexpected file(s) in results/paper: ",
                             paste(unexpected, collapse="; "))
message("Verified ", length(required), " required publication and AN files.")
