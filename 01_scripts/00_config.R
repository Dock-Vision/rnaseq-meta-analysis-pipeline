# ---------------------------------------------------------------------------
# 00_config.R
#
# THE ONLY FILE YOU NORMALLY NEED TO EDIT.
#
# Single source of truth for paths, the dataset registry, the gene panels and
# the analysis constants. Every other script starts with:
#
#     source("01_scripts/00_config.R")
#
# Nothing here runs an analysis; it only defines constants and creates the
# output folder tree.
#
# To point the pipeline at a new question you change three things and nothing
# else:
#   1. GENES_OF_INTEREST  -- the gene(s) the whole analysis is about
#   2. MARKERS            -- the cell-identity panel used to test whether a
#                            change is regulation or a composition shift
#   3. config/datasets.csv -- the studies to analyse (see datasets.example.csv)
#
# DockVision bulk RNA-seq meta-analysis pipeline.
# ---------------------------------------------------------------------------

# ---- reproducibility ------------------------------------------------------
# One seed for the whole project. Every stochastic step (sva/ComBat, ggrepel
# label placement, any bootstrap) inherits it, and it is written into every
# sessionInfo log so a result can always be traced back to the run that made it.

SEED <- 20260816L
set.seed(SEED)

# ---- paths ----------------------------------------------------------------
# Numbered so the folder listing reads in pipeline order. Raw downloads are
# never written into 04_processed/ -- keeping the as-downloaded file in
# 03_raw_data/ is what makes every downstream step reproducible from source.

PROJ <- normalizePath(getwd())
if (!dir.exists(file.path(PROJ, "01_scripts"))) {
  stop("Run scripts from the project root: ", PROJ, " is not it.")
}

PATHS <- list(
  scripts    = file.path(PROJ, "01_scripts"),
  config     = file.path(PROJ, "config"),
  metadata   = file.path(PROJ, "02_metadata"),
  raw        = file.path(PROJ, "03_raw_data"),
  raw_case   = file.path(PROJ, "03_raw_data", "case_control"),
  raw_atlas  = file.path(PROJ, "03_raw_data", "atlas"),
  processed  = file.path(PROJ, "04_processed"),
  results    = file.path(PROJ, "05_results"),
  de         = file.path(PROJ, "05_results", "deseq2"),
  coexpr     = file.path(PROJ, "05_results", "coexpression"),
  cellcomp   = file.path(PROJ, "05_results", "cell_composition"),
  fig_png    = file.path(PROJ, "06_figures", "png"),
  fig_pdf    = file.path(PROJ, "06_figures", "pdf"),
  fig_pub    = file.path(PROJ, "07_final_figures"),
  logs       = file.path(PROJ, "08_logs")
)

invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- genes of interest ----------------------------------------------------
# Stored as HGNC SYMBOLS only. Ensembl IDs are resolved at runtime via
# org.Hs.eg.db (see R/annotate.R) and cached, so no identifier is ever
# hardcoded from memory and silently wrong when a series ships a different ID
# space.
#
# The pipeline is written for ONE or TWO genes. With two, the first is treated
# as the primary gene and the second as its comparator/paralog, which is what
# the co-expression atlas (30_) and the coupling test (40_) correlate.
#
# EDIT ME -- the example below is the FLI1/ERG ETS-factor pair the pipeline was
# first built for.

GENES_OF_INTEREST <- c("FLI1", "ERG")

# ---- cell-identity marker panels ------------------------------------------
#
# THE CENTRAL SCIENTIFIC CONTROL OF THIS PIPELINE.
#
# In bulk RNA-seq a gene can appear to fall simply because the cells that
# express it are less abundant in the sample -- not because it is
# downregulated per cell. Every differential-expression result is therefore
# reported next to a panel of markers for the cell type that carries the gene,
# and the pipeline states which of the two explanations the data supports:
#
#   gene down AND markers all down together -> likely COMPOSITION SHIFT
#   gene down WHILE markers stay flat       -> likely GENUINE REGULATION
#
# `tissue_class` in the dataset registry decides which panel a study gets.
#
# EDIT ME -- the defaults are endothelial content markers for solid tissue and
# immune composition markers for whole blood.

MARKERS <- list(
  tissue = c("PECAM1", "CDH5", "VWF", "CLDN5"),   # endothelial content
  blood  = c("PTPRC", "FCGR3B", "CD3E", "PF4")    # immune composition
)

# Broader panels, used only for the composition sanity-check heatmaps.
MARKERS_EXTENDED <- list(
  endothelial = c("PECAM1", "CDH5", "VWF", "CLDN5", "TEK", "KDR", "ERG", "SOX17"),
  fibroblast  = c("COL1A1", "COL1A2", "ACTA2", "FN1", "TAGLN"),
  epithelial  = c("EPCAM", "KRT8", "KRT18", "SFTPC"),
  immune      = c("PTPRC", "CD3E", "CD19", "FCGR3B", "LYZ", "CD68"),
  platelet    = c("PF4", "PPBP", "ITGA2B")
)

# ---- dataset registry -----------------------------------------------------
# Read from config/datasets.csv (falling back to the shipped example). One row
# per GEO series; add a study by adding a row, never by writing a one-off
# per-dataset script.
#
# Columns:
#   gse           GEO series accession, e.g. "GSE112087"
#   condition     grouping label used for facets, palettes and subgroups
#                 (a disease name, a treatment, a cohort -- your choice)
#   use           PRIMARY / SECONDARY / VALIDATION / OPTIONAL. Free text; kept
#                 in the outputs so a reader can see which studies carried the
#                 result and which only replicated it.
#   tissue_class  "tissue" or "blood" -- selects the MARKERS panel above and
#                 forms the meta-analysis subgroups. Getting this wrong puts
#                 the wrong marker panel on a study, so check the GEO
#                 source_name field rather than assuming from the disease.
#   arm           "case_control" (differential expression) or
#                 "atlas" (co-expression across baseline samples only)
#   de_capable    TRUE/FALSE. FALSE for series with no control group -- they
#                 cannot support a case/control contrast and are skipped by
#                 10_run_dataset.R rather than silently producing a contrast
#                 against an arbitrary reference.

DATASET_REGISTRY <- local({
  user_csv    <- file.path(PATHS$config, "datasets.csv")
  example_csv <- file.path(PATHS$config, "datasets.example.csv")
  if (file.exists(user_csv)) {
    user_csv
  } else if (file.exists(example_csv)) {
    message("[config] config/datasets.csv not found -- using the shipped ",
            "example registry. Copy it to config/datasets.csv and edit it.")
    example_csv
  } else {
    stop("No dataset registry found. Expected ", user_csv)
  }
})

DATASETS <- utils::read.csv(DATASET_REGISTRY, stringsAsFactors = FALSE,
                            strip.white = TRUE, comment.char = "#")

local({
  required <- c("gse", "condition", "use", "tissue_class", "arm", "de_capable")
  missing  <- setdiff(required, names(DATASETS))
  if (length(missing)) {
    stop("Dataset registry ", DATASET_REGISTRY, " is missing column(s): ",
         paste(missing, collapse = ", "))
  }
  bad_tc <- setdiff(unique(DATASETS$tissue_class), names(MARKERS))
  if (length(bad_tc)) {
    stop("tissue_class value(s) with no matching MARKERS panel: ",
         paste(bad_tc, collapse = ", "))
  }
  bad_arm <- setdiff(unique(DATASETS$arm), c("case_control", "atlas"))
  if (length(bad_arm)) {
    stop("arm must be 'case_control' or 'atlas'; got: ",
         paste(bad_arm, collapse = ", "))
  }
  if (anyDuplicated(DATASETS$gse)) {
    stop("Duplicate GSE in the registry: ",
         paste(unique(DATASETS$gse[duplicated(DATASETS$gse)]), collapse = ", "))
  }
})

DATASETS$de_capable <- as.logical(DATASETS$de_capable)

# Backwards-compatible alias: several helpers still speak of "disease".
DATASETS$disease <- DATASETS$condition

ds_info <- function(gse) {
  row <- DATASETS[DATASETS$gse == gse, ]
  if (nrow(row) != 1) stop("Unknown GSE in registry: ", gse)
  as.list(row)
}

# Where a given GSE's raw files belong.
raw_dir_for <- function(gse) {
  info <- ds_info(gse)
  base <- if (info$arm == "case_control") PATHS$raw_case else PATHS$raw_atlas
  file.path(base, gse)
}

# Display labels for `condition`, used for figure facets and legends. Override
# any entry to get a prettier label than the registry code.
CONDITION_LABELS <- local({
  lv <- unique(DATASETS$condition[DATASETS$arm == "case_control"])
  setNames(lv, lv)
})
# e.g. CONDITION_LABELS["ILD"] <- "ILD / IPF"

# ---- manual count-column -> GSM maps --------------------------------------
# Last resort for series whose supplementary matrix uses sample names that
# cannot be matched to the GEO metadata by any string rule. Each entry must be
# justified by something other than column order alone -- record the reasoning,
# because a wrong mapping here silently mislabels samples and corrupts both the
# batch correction and the per-group statistics.
#
# Worked example of the standard expected (from a real series whose matrix used
# lab shorthand while the metadata used cell-type names -- the mapping was
# corroborated by cell-type semantics, not just position):
#
#   MANUAL_COLUMN_MAP <- list(
#     GSE128179 = c(
#       # dLyNeo = dermal LYMPHATIC neonatal -> the HLEC samples
#       hmvec_dLyNeo_1 = "GSM3666322", hmvec_dLyNeo_2 = "GSM3666323",
#       hmvec_dLyNeo_3 = "GSM3666324",
#       # haec / heac (a typo in the source file) -> the HAEC samples
#       haec_1         = "GSM3666325", heac_2         = "GSM3666326",
#       heac_3         = "GSM3666327")
#   )

MANUAL_COLUMN_MAP <- list()

# ---- explicit contrast definitions ----------------------------------------
# Many series carry MORE THAN TWO phenotype levels. A binary
# "control vs everything else" rule silently pools unrelated conditions -- a
# sepsis series that also contains a COVID-19 arm would put COVID patients into
# the sepsis case group; an IPF series with an acute-lung-injury arm would mix
# two pathologies.
#
# Where the phenotype column has more than two levels, state the contrast here.
# `case`/`control` are matched against the values of `col` (exact, after
# trimming). Samples matching NEITHER side are DROPPED, never absorbed into the
# case group.
#
# Worked examples:
#
#   CONTRAST_SPEC <- list(
#     # sepsis cohort that also contains a COVID-19 arm -> COVID excluded
#     GSE243217 = list(col  = "char_disease.state",
#                      case = "Sepsis", control = "healthy donor"),
#
#     # PAH validation cohort spanning WHO groups 1-5; only Group 1 is PAH
#     GSE243193 = list(col  = "char_ph.group.control",
#                      case = "Group 1 PAH", control = "Control")
#   )

CONTRAST_SPEC <- list()

# ---- analysis constants ---------------------------------------------------

ALPHA        <- 0.05      # FDR threshold, Benjamini-Hochberg
LFC_SHRINK   <- "apeglm"  # effect-size shrinkage; falls back to "ashr" when
                          # the contrast cannot be expressed as a model
                          # coefficient, and for the limma-trend path
MIN_COUNT    <- 10        # gene filter: >= MIN_COUNT reads ...
MIN_SAMPLES  <- NULL      # ... in >= smallest-group-size samples (computed)
NCPUS        <- max(1, parallel::detectCores() - 2)

# ---- misc -----------------------------------------------------------------

options(stringsAsFactors = FALSE)
options(timeout = max(1800, getOption("timeout")))   # GEO supplementary files
                                                     # are large and slow

message(sprintf("[config] project=%s  seed=%d  genes=%s  datasets=%d  alpha=%.2f",
                basename(PROJ), SEED, paste(GENES_OF_INTEREST, collapse = "/"),
                nrow(DATASETS), ALPHA))
