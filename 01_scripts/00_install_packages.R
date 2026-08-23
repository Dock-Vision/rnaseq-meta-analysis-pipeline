#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 00_install_packages.R
#
# Installs every package required by the analysis standard (docs/METHODS.md)
# into your system R user library.
#
# Run it OUTSIDE any conda environment. Installing the Bioconductor stack via
# conda shadows the system library and breaks DESeq2/sva/GSVA in ways that are
# hard to diagnose.
#
# Idempotent: already-installed packages are skipped, so this is safe to re-run
# on a fresh machine to reproduce the environment.
# ---------------------------------------------------------------------------

options(Ncpus = max(1, parallel::detectCores() - 2))
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)

LIB <- .libPaths()[1]
cat("Installing into:", LIB, "\n")
cat("R version:", R.version.string, "\n\n")

# ---- package manifest -----------------------------------------------------
# Each entry: package -> which requirement in docs/METHODS.md it serves.

cran_pkgs <- c(
  "remotes",            # infrastructure for the GitHub install below
  "metafor",            # METHODS S4 random-effects meta-analysis
  "ashr",               # METHODS S1 alternative LFC shrinkage
  "patchwork",          # METHODS S6 figure composition
  "cowplot",            # METHODS S6 figure composition
  "ggrepel",            # METHODS S6 non-overlapping gene labels
  "RColorBrewer",       # METHODS S6 consistent palettes
  "pheatmap"            # METHODS S6 (already present, kept for completeness)
)

bioc_pkgs <- c(
  "org.Hs.eg.db",       # symbol <-> Ensembl mapping (offline, no biomaRt calls)
  "apeglm",             # METHODS S1 primary LFC shrinkage estimator
  "variancePartition",  # METHODS S3 % variance batch vs biological group
  "EnhancedVolcano",    # METHODS S6 volcano plots
  "ComplexHeatmap",     # METHODS S6 heatmaps
  "PCAtools"            # METHODS S3/METHODS S6 PCA + variance explained
)

github_pkgs <- c(
  xCell = "dviraran/xCell"   # METHODS S2 deconvolution; not on CRAN/Bioconductor
)

# ---- helpers ---------------------------------------------------------------

have <- function(p) requireNamespace(p, quietly = TRUE)

report <- function(pkgs) {
  for (p in pkgs) {
    cat(sprintf("  %-20s %s\n", p,
                if (have(p)) as.character(packageVersion(p)) else "** FAILED **"))
  }
}

install_if_missing <- function(pkgs, installer, label) {
  todo <- pkgs[!vapply(pkgs, have, logical(1))]
  cat("\n===", label, "===\n")
  if (length(todo) == 0) {
    cat("All present, nothing to do.\n")
    return(invisible())
  }
  cat("To install:", paste(todo, collapse = ", "), "\n\n")
  for (p in todo) {
    cat("---- installing", p, "----\n")
    tryCatch(installer(p),
             error = function(e) cat("ERROR installing ", p, ": ",
                                     conditionMessage(e), "\n", sep = ""))
  }
}

# ---- install ---------------------------------------------------------------

install_if_missing(cran_pkgs,
                   function(p) install.packages(p, lib = LIB),
                   "CRAN")

install_if_missing(bioc_pkgs,
                   function(p) BiocManager::install(p, lib = LIB,
                                                    update = FALSE, ask = FALSE),
                   "Bioconductor")

cat("\n=== GitHub ===\n")
for (nm in names(github_pkgs)) {
  if (have(nm)) {
    cat(nm, "already present.\n")
    next
  }
  cat("---- installing", nm, "from", github_pkgs[[nm]], "----\n")
  tryCatch(
    remotes::install_github(github_pkgs[[nm]], lib = LIB,
                            upgrade = "never", dependencies = TRUE),
    error = function(e) cat("ERROR installing ", nm, ": ",
                            conditionMessage(e), "\n", sep = "")
  )
}

# ---- verify ----------------------------------------------------------------

cat("\n\n", strrep("=", 70), "\nFINAL VERIFICATION\n", strrep("=", 70), "\n",
    sep = "")

cat("\nCore stack (pre-existing):\n")
report(c("GEOquery", "DESeq2", "sva", "limma", "biomaRt", "data.table",
         "readxl", "ggplot2"))

cat("\nAnalysis-standard additions:\n")
report(c(cran_pkgs, bioc_pkgs, names(github_pkgs)))

all_req <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))
failed <- all_req[!vapply(all_req, have, logical(1))]
cat("\n")
if (length(failed) == 0) {
  cat("SUCCESS: all required packages installed.\n")
} else {
  cat("STILL MISSING:", paste(failed, collapse = ", "), "\n")
}

# ---- reproducibility: log sessionInfo -------------------------------------
si_path <- file.path("08_logs", "sessionInfo_00_install_packages.txt")
dir.create("08_logs", showWarnings = FALSE)
capture.output(sessionInfo(), file = si_path)
cat("\nsessionInfo written to", si_path, "\n")
