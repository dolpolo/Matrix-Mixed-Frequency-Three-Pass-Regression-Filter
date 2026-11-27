# Validate paper and mandatory AN outputs against the explicit source manifest.
manifest <- utils::read.csv("config/output_manifest.csv", stringsAsFactors = FALSE)

# Convert the CSV required flag to a logical vector before indexing.
required_flag <- tolower(trimws(as.character(manifest$required)))
if (length(required_flag) != nrow(manifest) ||
    anyNA(required_flag) ||
    any(!required_flag %in% c("true", "false"))) {
  stop("Invalid required flag in config/output_manifest.csv")
}

required <- manifest$path[required_flag == "true"]
missing <- required[!file.exists(required)]

if (length(missing)) {
  stop("Missing publication/diagnostic outputs: ",
       paste(missing, collapse = "; "))
}

# Intermediate CSVs are allowed in results/audit, but not publication folders.
paper_dir <- "results/paper"
present <- list.files(paper_dir, full.names = TRUE, recursive = TRUE,
                      all.files = FALSE)
present <- present[file.info(present)$isdir == FALSE]
expected_paper <- manifest$path[grepl("^results/paper/", manifest$path)]
unexpected <- setdiff(present, expected_paper)

if (length(unexpected)) {
  stop("Unexpected file(s) in results/paper: ",
       paste(unexpected, collapse = "; "))
}

message("Verified ", length(required), " required publication and AN files.")
