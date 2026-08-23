#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00b_install_extras.R -- second-pass installer for the two packages that
# routinely fail in 00_install_packages.R. Run it after that script.
#
# 1. org.Hs.eg.db  -- a ~100 MB annotation tarball that regularly exceeds the
#                     default download timeout. Retried here with a 2-hour one.
#
# 2. xCell         -- install_github(dependencies = TRUE) also pulls Suggests
#                     (fs, magick, rmarkdown, testthat), and fs/magick need
#                     system libraries that are often absent and need root to
#                     add. Those are only used for vignettes, so GSVA is
#                     installed from Bioconductor directly and xCell with hard
#                     dependencies only.
#
# If GSVA still fails with "dependency not available", install the system
# library it needs first:  sudo apt-get install -y libmagick++-dev
# ---------------------------------------------------------------------------

options(timeout = 7200)
options(Ncpus = max(1, parallel::detectCores() - 2))
options(repos = c(CRAN = "https://cloud.r-project.org"))
LIB <- .libPaths()[1]

have <- function(p) requireNamespace(p, quietly = TRUE)

# ---- 1. org.Hs.eg.db ------------------------------------------------------

if (!have("org.Hs.eg.db")) {
  cat("\n=== installing org.Hs.eg.db (large annotation package) ===\n")
  tryCatch(BiocManager::install("org.Hs.eg.db", lib = LIB,
                                update = FALSE, ask = FALSE),
           error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
} else {
  cat("org.Hs.eg.db already present.\n")
}

# ---- 2. GSVA then xCell ---------------------------------------------------

if (!have("GSVA")) {
  cat("\n=== installing GSVA (xCell's engine) ===\n")
  tryCatch(BiocManager::install("GSVA", lib = LIB, update = FALSE, ask = FALSE),
           error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

if (!have("xCell")) {
  cat("\n=== installing xCell from GitHub (hard dependencies only) ===\n")
  tryCatch(
    remotes::install_github("dviraran/xCell", lib = LIB, upgrade = "never",
                            dependencies = c("Depends", "Imports"),
                            build_vignettes = FALSE),
    error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

# ---- verify ---------------------------------------------------------------

cat("\n", strrep("=", 60), "\nVERIFICATION\n", strrep("=", 60), "\n", sep = "")
for (p in c("org.Hs.eg.db", "GSVA", "xCell")) {
  cat(sprintf("  %-16s %s\n", p,
              if (have(p)) as.character(packageVersion(p)) else "** MISSING **"))
}

dir.create("08_logs", showWarnings = FALSE)
capture.output(sessionInfo(), file = "08_logs/sessionInfo_00b_install_extras.txt")
