#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 60_celltype_atlas.R -- which cell populations does the gene track with?
#
# WHY THIS EXISTS
# ---------------
# 10_run_dataset.R asks whether a gene changes between groups, and checks that
# change against a hand-picked marker panel. This script asks the complementary
# question across the whole cell-type space at once: as samples get richer in
# each of ~64 cell populations, does the gene go up or down?
#
# Bulk RNA-seq cannot measure expression per cell type -- the library averages
# over every cell in the sample. The deconvolution route can get close, and the
# inputs already exist: 10_run_dataset.R runs xCell on every cohort (METHODS
# S2.4), giving a per-sample enrichment score for each cell type. Correlating
# expression against each of those scores, within a cohort, answers "the gene
# rises and falls with which cell populations?" -- a weaker route to the same
# biological question, and one whose weakness is stated rather than hidden.
#
# WHAT THE NUMBER MEANS -- read this before interpreting the output
#   r > 0 : samples richer in that cell type express more of the gene.
#   It is an ASSOCIATION ACROSS SAMPLES, not expression measured inside that
#   cell type. A high r is consistent with the gene being expressed in that
#   population, but composition, disease status and any correlated cell type
#   can all produce it. This is the composition confound of METHODS S2 turned
#   into the measurement, instead of being treated as noise.
#
# METHOD
#   per cohort  : Pearson r between VST expression and each xCell score across
#                 all samples. Cell types with more than 25% exact-zero scores
#                 or sd < 1e-3 are flagged unreliable and excluded from testing
#                 (the existing rule in R/deconv.R). BH-adjusted within cohort.
#   across      : one random-effects model per gene x cell type (metafor,
#                 Fisher z, REML), tau^2 / I^2 / Q reported (METHODS S4). Never
#                 a simple average. BH-adjusted across the cell types tested
#                 for that gene.
#
# Requires 10_run_dataset.R to have been run WITH xCell (i.e. without
# --skip-xcell) on at least MIN_K cohorts.
#
# Run from the project root:
#   Rscript 01_scripts/60_celltype_atlas.R 2>&1 \
#     | tee 08_logs/60_celltype_atlas.log
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source("01_scripts/00_config.R")
for (f in c("utils", "annotate", "deconv", "meta")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

SCRIPT <- "60_celltype_atlas"
set.seed(SEED)

OUT_STUDY  <- file.path(PATHS$cellcomp, "celltype_atlas_per_study.csv")
OUT_POOLED <- file.path(PATHS$results,  "meta_celltype_atlas_pooled.csv")

GENES <- GENES_OF_INTEREST

# A cell type needs at least this many reliable cohort estimates to be pooled.
MIN_K <- 3

# xCell emits three COMPOSITE scores alongside the cell types. They are sums of
# other columns, so including them would double-count and produce a guaranteed
# strong "correlation" with whatever dominates. Cell types only.
COMPOSITE <- c("ImmuneScore", "StromaScore", "MicroenvironmentScore")

# ---- lineage families ------------------------------------------------------
# xCell's cell types grouped into five lineages, used for the colour blocks in
# 61_celltype_atlas_figure.R. Five, not nine, because a categorical palette only
# stays colourblind-safe to about five hues; the finer split survives as the
# axis labels, which name every individual cell type.
#
# Override by defining CELL_FAMILIES in 00_config.R before this runs -- e.g. if
# you swap xCell for another deconvolution method with a different vocabulary.
# Any cell type not listed here is grouped as "Other" rather than aborting the
# run, so a new xCell release cannot break the script.
if (!exists("CELL_FAMILIES")) {
  CELL_FAMILIES <- list(
    `Structural & parenchymal` = c(
      "Astrocytes", "Neurons", "Epithelial cells", "Keratinocytes", "Sebocytes",
      "Hepatocytes", "Melanocytes", "Fibroblasts", "MSC", "Pericytes",
      "Chondrocytes", "Osteoblast", "Adipocytes", "Preadipocytes",
      "Mesangial cells", "Myocytes", "Skeletal muscle", "Smooth muscle"),
    Endothelial = c(
      "Endothelial cells", "ly Endothelial cells", "mv Endothelial cells"),
    Myeloid = c(
      "aDC", "cDC", "iDC", "pDC", "DC", "Macrophages", "Macrophages M1",
      "Macrophages M2", "Monocytes", "Neutrophils", "Eosinophils", "Basophils",
      "Mast cells"),
    Lymphoid = c(
      "B-cells", "naive B-cells", "Memory B-cells", "Class-switched memory B-cells",
      "pro B-cells", "Plasma cells", "CD4+ T-cells", "CD4+ naive T-cells",
      "CD4+ memory T-cells", "CD4+ Tcm", "CD4+ Tem", "CD8+ T-cells",
      "CD8+ naive T-cells", "CD8+ Tcm", "CD8+ Tem", "Th1 cells", "Th2 cells",
      "Tregs", "Tgd cells", "NK cells", "NKT"),
    `HSC, progenitor & Mk/erythroid` = c(
      "HSC", "MPP", "CLP", "CMP", "GMP", "MEP",
      "Megakaryocytes", "Erythrocytes", "Platelets"))
}
OTHER_FAMILY <- "Other"

FAMILY_LEVELS <- names(CELL_FAMILIES)
CELL_FAMILY <- setNames(rep(FAMILY_LEVELS, lengths(CELL_FAMILIES)),
                        unlist(CELL_FAMILIES, use.names = FALSE))
if (anyDuplicated(names(CELL_FAMILY))) {
  stop("CELL_FAMILIES lists the same cell type in more than one family: ",
       paste(unique(names(CELL_FAMILY)[duplicated(names(CELL_FAMILY))]),
             collapse = ", "))
}

# ---- per-cohort correlations ----------------------------------------------

rds <- sort(Sys.glob(file.path(PATHS$processed, "*_pipeline_result.rds")))
if (!length(rds)) {
  stop("No *_pipeline_result.rds found in ", PATHS$processed,
       ".\nRun 01_scripts/10_run_dataset.R on your datasets first.")
}
log_step("scanning %d pipeline result(s)", length(rds))

per_study <- list()
seen_types <- character(0)
for (f in rds) {
  r <- readRDS(f)
  if (is.null(r$xcell) || is.null(r$vst)) {
    log_warn("%s: no xCell scores or no VST matrix -- skipped", r$gse); next
  }

  ct_here <- setdiff(rownames(r$xcell), COMPOSITE)
  seen_types <- union(seen_types, ct_here)

  tab <- correlate_gene_with_xcell(r$vst, r$xcell, genes = GENES,
                                   cell_types = ct_here)
  if (is.null(tab)) { log_warn("%s: too few samples to correlate", r$gse); next }

  tab$study <- r$gse
  tab$tissue_class <- r$info$tissue_class
  tab$condition <- r$info$condition
  per_study[[r$gse]] <- tab

  log_info("%s: n=%d samples, %d cell types tested, %d reliable",
           r$gse, ncol(r$vst), length(ct_here),
           sum(tab$reliable[tab$gene == GENES[1]]))
}

if (!length(per_study)) {
  stop("No cohort produced xCell correlations. Re-run 10_run_dataset.R without ",
       "--skip-xcell.")
}

# Cell types the family map does not know about are grouped rather than fatal:
# a new xCell release adding a population should not stop the analysis.
unknown <- setdiff(seen_types, names(CELL_FAMILY))
if (length(unknown)) {
  log_warn("%d cell type(s) not in CELL_FAMILIES -- grouped as '%s': %s",
           length(unknown), OTHER_FAMILY, paste(unknown, collapse = ", "))
  CELL_FAMILY <- c(CELL_FAMILY, setNames(rep(OTHER_FAMILY, length(unknown)),
                                         unknown))
  FAMILY_LEVELS <- c(FAMILY_LEVELS, OTHER_FAMILY)
}

# Plot/report order: family blocks contiguous, alphabetical inside each block,
# and restricted to what was actually scored in this run.
CELL_ORDER <- unlist(lapply(FAMILY_LEVELS, function(fam) {
  sort(intersect(names(CELL_FAMILY)[CELL_FAMILY == fam], seen_types))
}), use.names = FALSE)
log_info("%d cell type(s) scored across %d family blocks",
         length(CELL_ORDER), length(unique(CELL_FAMILY[CELL_ORDER])))

study_dt <- rbindlist(per_study, fill = TRUE)
study_dt[, family := CELL_FAMILY[cell_type]]
study_dt[, family := factor(family, levels = FAMILY_LEVELS)]
stopifnot(!anyNA(study_dt$family))

fwrite(study_dt, OUT_STUDY)
log_info("wrote %s (%d rows, %d studies)", OUT_STUDY, nrow(study_dt),
         uniqueN(study_dt$study))

# ---- pool across cohorts ---------------------------------------------------
# One random-effects model per gene x cell type. Only reliable cohort-level
# estimates enter; a cell type needs at least MIN_K of them to be pooled at all.

pooled <- list()
for (g in GENES) {
  for (ct in CELL_ORDER) {
    d <- study_dt[gene == g & cell_type == ct & reliable == TRUE & !is.na(r)]
    if (nrow(d) < MIN_K) next
    fit <- pool_correlations(
      data.frame(subgroup = "all", r = d$r, n = d$n, study = d$study),
      subgroup = "all")
    if (is.null(fit)) next
    s <- fit$summary
    s$gene <- g; s$cell_type <- ct; s$family <- CELL_FAMILY[[ct]]
    pooled[[paste(g, ct)]] <- s
  }
}

if (!length(pooled)) {
  stop("No cell type reached k >= ", MIN_K, " cohorts. Run more datasets, or ",
       "lower MIN_K if you accept the weaker pooling.")
}

pooled_dt <- rbindlist(pooled, fill = TRUE)
# Multiple testing is across the cell types tested FOR A GENE, so adjust within
# gene -- not across the combined table for every gene at once.
pooled_dt[, p_adj := p.adjust(p_value, method = "BH"), by = gene]
setorder(pooled_dt, gene, -pooled_r)

fwrite(pooled_dt, OUT_POOLED)
log_info("wrote %s (%d gene x cell type models)", OUT_POOLED, nrow(pooled_dt))

for (g in GENES) {
  d <- pooled_dt[gene == g]
  if (!nrow(d)) { log_warn("%s: nothing pooled", g); next }
  cat(sprintf("\n==== %s: top pooled associations (k = cohorts) ====\n", g))
  print(head(d[, .(cell_type, family, k,
                   r = round(pooled_r, 3),
                   ci = sprintf("[%.2f, %.2f]", ci_low_r, ci_high_r),
                   padj = signif(p_adj, 2),
                   I2 = round(I2))], 12))
  cat(sprintf("%s: %d/%d cell types significant at BH-FDR < %.2f\n",
              g, sum(d$p_adj < ALPHA, na.rm = TRUE), nrow(d), ALPHA))
}

append_progress(
  "ALL (cell-type atlas)", "Gene vs cell-type association atlas", "DONE",
  sprintf("%d cohort(s), %d cell type(s) scored, %d gene x cell type models pooled at k >= %d. %s",
          uniqueN(study_dt$study), length(CELL_ORDER), nrow(pooled_dt), MIN_K,
          paste(vapply(GENES, function(g) {
            d <- pooled_dt[gene == g]
            sprintf("%s: %d/%d significant at BH-FDR<%.2f",
                    g, sum(d$p_adj < ALPHA, na.rm = TRUE), nrow(d), ALPHA)
          }, character(1)), collapse = "; ")),
  OUT_POOLED)

log_session_info(SCRIPT)
log_step("COMPLETE: cell-type association atlas")
