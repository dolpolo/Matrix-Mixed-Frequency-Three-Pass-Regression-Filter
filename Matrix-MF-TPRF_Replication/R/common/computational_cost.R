# ==============================================================================
# MATRIX THREE-PASS REGRESSION FILTER:
# A NOWCASTING FRAMEWORK FOR MIXED-FREQUENCY MULTIDIMENSIONAL DATA
# ==============================================================================
# Authors:      Matteo Barigozzi, Davide Delfino, Massimiliano Marcellino
# Institutions: University of Bologna, Aarhus University, Bocconi University
# Contact:      davide.delfino@econ.au.dk
# File:         R/common/computational_cost.R
# Input:        Runner environment and elapsed times from empirical scripts.
# Output:       Timing records for Matrix MF-TPRF, VEC-P and DFM.
# ==============================================================================

# 1. Benchmark configuration ---------------------------------------------------

computational_components <- c("matrix", "vecp", "dfm")
computational_records_file <- file.path(
  "results", "audit", "computational_cost", "timing_records.rds"
)

computational_active <- function() {
  nzchar(Sys.getenv("MF_TPRF_TIMING_DIR", unset = ""))
}

computational_worker_budget <- function() {
  value <- Sys.getenv("MF_TPRF_WORKERS", unset = "8")
  if (!grepl("^[2-8]$", value)) {
    stop("MF_TPRF_WORKERS must be an integer between 2 and 8.")
  }
  workers <- as.integer(value)
  available <- parallelly::availableCores()
  if (workers > available) {
    stop("Requested ", workers, " workers but only ", available,
         " CPU slots are available. Set MF_TPRF_WORKERS to a feasible value.")
  }
  workers
}

computational_workers <- function() {
  if (computational_active()) {
    return(as.integer(Sys.getenv("MF_TPRF_WORKERS")))
  }
  # Use up to eight workers in standalone model runs.
  detected <- parallel::detectCores()
  if (is.na(detected)) detected <- 2L
  max(1L, min(8L, detected - 1L))
}

# 2. Timing environment --------------------------------------------------------

# These settings are inherited by fresh R sessions and their future workers.
# They request one numerical-library thread per process; actual support depends
# on the BLAS/OpenMP implementation. They are not a CPU-affinity restriction.
computational_thread_settings <- c(
  OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1", MKL_DYNAMIC = "FALSE"
)

with_computational_environment <- function(directory, workers, code) {
  settings <- c(
    MF_TPRF_TIMING_DIR = directory,
    MF_TPRF_WORKERS = as.character(workers),
    computational_thread_settings
  )
  old <- Sys.getenv(names(settings), unset = NA_character_)
  names(old) <- names(settings)
  on.exit({
    Sys.unsetenv(names(settings))
    present <- !is.na(old)
    if (any(present)) do.call(Sys.setenv, as.list(old[present]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(settings))
  force(code)
}

# 3. Record elapsed wall-clock time -------------------------------------------

computational_clock <- function() {
  unname(proc.time()[["elapsed"]])
}

computational_record <- function(model, elapsed, workers, phases,
                                 start_eval, end_eval, targets,
                                 predictor_counts, end_by_target = NULL) {
  if (!computational_active()) return(invisible(NULL))
  if (!is.finite(elapsed) || elapsed <= 0) stop("Invalid timing for ", model)
  directory <- Sys.getenv("MF_TPRF_TIMING_DIR")
  if (!dir.exists(directory)) stop("Missing benchmark timing directory.")
  budget <- as.integer(Sys.getenv("MF_TPRF_WORKERS"))
  if (model != "matrix" && workers != budget) {
    stop("Unexpected number of future workers for ", model)
  }
  info <- Sys.info()
  record <- list(
    model = model,
    elapsed_seconds = elapsed,
    workers = as.integer(workers),
    worker_budget = budget,
    phases = phases,
    start_eval = as.character(start_eval),
    end_eval = as.character(end_eval),
    targets = targets,
    predictor_counts = predictor_counts,
    end_by_target = end_by_target,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    platform = unname(info[c("sysname", "release", "machine")]),
    cpu = Sys.getenv("PROCESSOR_IDENTIFIER", unset = ""),
    logical_cpus = parallel::detectCores(logical = TRUE),
    physical_cpus = parallel::detectCores(logical = FALSE),
    r_version = R.version.string,
    blas = sessionInfo()$BLAS,
    lapack = sessionInfo()$LAPACK,
    thread_settings = Sys.getenv(names(computational_thread_settings))
  )
  saveRDS(record, file.path(directory, paste0(model, ".rds")))
  message("Computational cost (", model, "): ", round(elapsed, 2), " seconds")
  invisible(record)
}

# 4. Collect one completed run -------------------------------------------------

collect_computational_records <- function(directory) {
  paths <- file.path(directory, paste0(computational_components, ".rds"))
  if (!all(file.exists(paths))) {
    stop("Missing timing records: ", paste(basename(paths[!file.exists(paths)]),
                                           collapse = ", "))
  }
  records <- setNames(lapply(paths, readRDS), computational_components)
  stopifnot(identical(vapply(records, function(x) x$model, character(1)),
                      setNames(computational_components, computational_components)))
  budgets <- vapply(records, function(x) x$worker_budget, integer(1))
  workers <- vapply(records, function(x) x$workers, integer(1))
  elapsed <- vapply(records, function(x) x$elapsed_seconds, numeric(1))
  starts <- vapply(records, function(x) x$start_eval, character(1))
  ends <- vapply(records, function(x) x$end_eval, character(1))
  if (length(unique(budgets)) != 1L || workers[["matrix"]] != 1L ||
      any(workers[c("vecp", "dfm")] != budgets[c("vecp", "dfm")]) ||
      any(!is.finite(elapsed) | elapsed <= 0) ||
      length(unique(starts)) != 1L || length(unique(ends)) != 1L) {
    stop("Timing records do not describe a comparable completed run.")
  }
  dir.create(dirname(computational_records_file), recursive = TRUE,
             showWarnings = FALSE)
  saveRDS(records, computational_records_file)
  invisible(records)
}
