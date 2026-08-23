# ---------------------------------------------------------------------------
# R/utils.R -- console logging, session capture, and the automatic run-log
# append.
#
# Sourced by every stage script. Nothing here is analysis-specific.
# ---------------------------------------------------------------------------

# ---- console logging ------------------------------------------------------

.log <- function(level, fmt, ...) {
  msg <- if (length(list(...))) sprintf(fmt, ...) else fmt
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%H:%M:%S"), level, msg))
  invisible(msg)
}

log_info <- function(fmt, ...) .log("INFO", fmt, ...)
log_warn <- function(fmt, ...) .log("WARN", fmt, ...)
log_step <- function(fmt, ...) {
  cat("\n", strrep("-", 74), "\n", sep = "")
  .log("STEP", fmt, ...)
  cat(strrep("-", 74), "\n", sep = "")
}

# ---- reproducibility ------------------------------------------------------

# Call at the end of every script. Writes sessionInfo() plus the seed to
# 08_logs/ so the exact package versions behind each result are recoverable
# years later, which is what makes a published figure auditable.
log_session_info <- function(script_name) {
  path <- file.path(PATHS$logs, paste0("sessionInfo_", script_name, ".txt"))
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(c(
    paste("script      :", script_name),
    paste("run at      :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("seed        :", if (exists("SEED")) SEED else NA),
    paste("R home      :", R.home()),
    paste("libPaths    :", paste(.libPaths(), collapse = " ; ")),
    paste("conda prefix:", Sys.getenv("CONDA_PREFIX", "<none>")),
    "", strrep("=", 70), ""
  ), con)
  capture.output(sessionInfo(), file = con)
  log_info("sessionInfo -> %s", path)
  invisible(path)
}

# ---- automatic run-log append ---------------------------------------------

# Every stage that produces or changes a file in 04_processed/, 05_results/,
# 06_figures/ or 07_final_figures/ appends a row to 08_logs/RUN_LOG.md. Doing
# it in code rather than by hand means it also happens on an unattended re-run,
# so the log is always a true record of what produced the current outputs.
#
# The notes field is expected to carry the numbers a reader would otherwise
# have to recompute: matrix dimensions, group sizes, whether the values were
# raw counts, and anything surprising about the dataset.
RUN_LOG_HEADER <- c(
  "# Run log",
  "",
  "Appended automatically by the pipeline. One row per completed step.",
  "",
  "| Date | Dataset | Step | Status | Notes | Output file path |",
  "|---|---|---|---|---|---|"
)

append_progress <- function(dataset, step, status, notes, output,
                            date = format(Sys.Date(), "%Y-%m-%d")) {
  stopifnot(status %in% c("DONE", "IN PROGRESS", "BLOCKED",
                          "NEEDS REVIEW", "SKIPPED"))
  path <- file.path(PATHS$logs, "RUN_LOG.md")
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(RUN_LOG_HEADER, path)
  }
  clean <- function(x) gsub("\\|", "/", gsub("[\r\n]+", " ", as.character(x)))
  row <- sprintf("| %s | %s | %s | %s | %s | %s |",
                 date, clean(dataset), clean(step), status,
                 clean(notes), clean(output))

  # Idempotent: re-running a dataset should UPDATE its row, not append a near
  # duplicate with slightly different numbers. Match on dataset + step so the
  # log stays a current-state summary rather than a growing pile.
  lines <- readLines(path, warn = FALSE)
  key <- sprintf("| %s | %s | ", clean(dataset), clean(step))
  is_old <- startsWith(lines, "| ") &
    grepl(key, lines, fixed = TRUE) &
    !startsWith(lines, "| Date")
  n_replaced <- sum(is_old)
  if (n_replaced > 0) lines <- lines[!is_old]

  writeLines(c(lines, row), path)
  log_info("RUN_LOG %s %s / %s [%s]",
           if (n_replaced > 0) sprintf("~ (replaced %d)", n_replaced) else "+=",
           dataset, step, status)
  invisible(row)
}

# ---- small helpers --------------------------------------------------------

# Format an integer with thousands separators, for log lines and captions.
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# Write a data.frame to 05_results/... with a consistent naming convention.
save_result <- function(df, subdir, name) {
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(subdir, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  log_info("saved %s (%d rows)", path, nrow(df))
  invisible(path)
}

# Guard used by pipeline stages that require an optional package. Gives an
# actionable message instead of a bare "there is no package called ..." error.
require_pkg <- function(pkg, needed_for) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required for %s but is not installed.\n  Run: Rscript 01_scripts/00_install_packages.R",
      pkg, needed_for), call. = FALSE)
  }
  invisible(TRUE)
}
