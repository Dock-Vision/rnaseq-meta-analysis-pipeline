#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 40_sensitivity_tests.R
#
# Is a pooled effect genuine biology, or an artifact of measuring a gene that
# sits near the detection floor in that tissue?
#
# This is the script to run when the meta-analysis hands you a significant
# result you do not believe -- typically a gene that is barely expressed in the
# tissue being tested (an endothelial transcription factor in whole blood, say)
# yet comes out strongly and significantly changed. Three tests, all reported
# with statistics rather than asserted:
#
#   TEST 1  Detection-restricted re-test. Keep only samples where the gene is
#           reliably detected, re-run the contrast, and see whether the effect
#           survives. A real signal persists, or strengthens, once unreliable
#           near-zero measurements are dropped. If too few samples remain to
#           re-test at all, that is itself the answer.
#
#   TEST 2  Which compartment does the gene track? Correlate it against each
#           cell-identity panel in MARKERS in the same samples. A signal that
#           is genuinely coming from a cell type should track that cell type's
#           markers and not a competing compartment's.
#
#   TEST 3  Paralog / comparator coupling. Where two genes are tightly
#           co-expressed in the reference tissue (30_coexpression_atlas.R
#           measures exactly this), a change driven by cell-type abundance must
#           move BOTH. One gene jumping while its partner stays flat is not
#           compatible with a composition explanation.
#
# Usage:
#   Rscript 01_scripts/40_sensitivity_tests.R
#   Rscript 01_scripts/40_sensitivity_tests.R --gene ERG --partner FLI1 \
#                                             --tissue-class blood
#
# Requires 10_run_dataset.R to have been run on the datasets being tested.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(DESeq2); library(SummarizedExperiment)
})

source("01_scripts/00_config.R")
source("01_scripts/00_theme.R")
for (f in c("utils", "annotate", "prepare", "de", "meta", "figures")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

SCRIPT <- "40_sensitivity_tests"
suppressPackageStartupMessages(library(data.table))

# ---- arguments -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

GENE    <- arg_val("--gene",    GENES_OF_INTEREST[1])
PARTNER <- arg_val("--partner", if (length(GENES_OF_INTEREST) > 1)
                                  GENES_OF_INTEREST[2] else NA_character_)
TCLASS  <- arg_val("--tissue-class", "blood")

if (is.na(PARTNER)) {
  stop("TEST 3 needs a comparator gene. Pass --partner <SYMBOL>, or define a ",
       "second gene in GENES_OF_INTEREST.")
}
if (!TCLASS %in% DATASETS$tissue_class) {
  stop("No datasets in the registry with tissue_class = '", TCLASS, "'.")
}

# Every case/control dataset of that tissue class that has actually been run.
TEST_GSES <- DATASETS$gse[DATASETS$tissue_class == TCLASS &
                          DATASETS$arm == "case_control" &
                          DATASETS$de_capable]

log_info("testing %s (partner %s) across %d '%s' dataset(s): %s",
         GENE, PARTNER, length(TEST_GSES), TCLASS,
         paste(TEST_GSES, collapse = ", "))

# The coupling benchmark for TEST 3 is MEASURED by 30_coexpression_atlas.R, not
# assumed. If the atlas has not been run, the test still reports the observed
# correlation, just without a reference value to compare it against.
ATLAS_R <- local({
  f <- file.path(PATHS$coexpr, "atlas_pooled_correlation.csv")
  if (!file.exists(f)) return(NA_real_)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  col <- intersect(c("pooled_r", "pooled_estimate", "estimate"), names(d))
  if (!length(col) || !nrow(d)) NA_real_ else as.numeric(d[[col[1]]][1])
})
if (!is.na(ATLAS_R)) {
  log_info("atlas coupling benchmark: pooled r = %.3f", ATLAS_R)
} else {
  log_warn("no atlas correlation found -- TEST 3 will report r without a benchmark")
}

results_t1 <- list(); results_t2 <- list(); results_t3 <- list()

for (gse in TEST_GSES) {
  rds <- file.path(PATHS$processed, paste0(gse, "_pipeline_result.rds"))
  if (!file.exists(rds)) { log_warn("missing %s", rds); next }
  r <- readRDS(rds)
  vsd <- r$vst
  meta <- r$metadata
  if (is.null(meta)) meta <- r$effect   # fallback
  # rebuild the group vector from the stored effect table's design
  dds <- attr(r$de, "dds")
  if (!is.null(dds)) {
    grp <- SummarizedExperiment::colData(dds)$group
    names(grp) <- colnames(dds)
  } else {
    # limma-path datasets store no DESeqDataSet; rebuild the grouping from the
    # saved metadata using the same column the pipeline chose
    md <- utils::read.csv(file.path(PATHS$metadata, paste0(gse, "_metadata.csv")),
                          stringsAsFactors = FALSE)
    if (is.null(md$sample_id)) md$sample_id <- md$gsm_id
    gcol <- r$group_col
    if (is.null(gcol) || is.null(md[[gcol]])) {
      log_warn("[%s] no group column recoverable -- skipping", gse); next
    }
    grp <- assign_groups(md, gcol, verbose = FALSE)
    names(grp) <- md$sample_id
    log_info("[%s] groups rebuilt from metadata column '%s'", gse, gcol)
  }
  common <- intersect(colnames(vsd), names(grp))
  vsd <- vsd[, common, drop = FALSE]
  grp <- droplevels(factor(grp[common]))

  panel <- resolve_panel(c(GENE, PARTNER), rownames(vsd))
  if (!all(panel$found)) {
    log_warn("[%s] %s/%s not in the matrix -- skipping", gse, GENE, PARTNER)
    next
  }
  expr_gene    <- as.numeric(vsd[panel$id[panel$symbol == GENE], ])
  expr_partner <- as.numeric(vsd[panel$id[panel$symbol == PARTNER], ])

  # ---- TEST 1: detection-restricted re-test ------------------------------
  log_step("[%s] TEST 1 -- detection-restricted re-test", gse)

  # Use the SAME detection basis as the per-run gene detection report: the
  # linear-scale expression matrix (CPM for count data, RPKM/TPM as shipped),
  # with a threshold of 1 unit. Anchoring to the published detection numbers
  # keeps this test interpretable rather than inventing a new floor.
  cpath <- file.path(PATHS$processed, paste0(gse, "_counts_collapsed.csv"))
  if (!file.exists(cpath)) cpath <- file.path(PATHS$processed,
                                              paste0(gse, "_counts.csv"))
  cm <- data.table::fread(cpath, data.table = FALSE)
  rownames(cm) <- cm[[1]]; cm <- as.matrix(cm[, -1, drop = FALSE])
  epanel <- resolve_panel(GENE, rownames(cm))
  if (!epanel$found[1]) { log_warn("[%s] %s absent from counts", gse, GENE); next }
  gene_lin <- cm[epanel$id[1], ]
  if (r$de_engine == "deseq2") {
    gene_lin <- gene_lin / (colSums(cm) / 1e6)   # CPM
  }
  gene_lin <- gene_lin[intersect(names(gene_lin), common)]
  detected <- rep(FALSE, length(common))
  names(detected) <- common
  detected[names(gene_lin)] <- gene_lin > 1
  log_info("%s detected (>1 %s) in %d/%d samples (%.0f%%)", GENE,
           if (r$de_engine == "deseq2") "CPM" else "RPKM/TPM",
           sum(detected), length(detected), 100 * mean(detected))

  t1 <- data.frame(study = gse, gene = GENE, n_total = length(expr_gene),
                   n_detected = sum(detected),
                   frac_detected = mean(detected),
                   stringsAsFactors = FALSE)

  if (sum(detected) >= 8 && nlevels(droplevels(grp[detected])) == 2 &&
      min(table(grp[detected])) >= 3) {
    tt <- stats::t.test(expr_gene[detected] ~ grp[detected])
    tt_all <- stats::t.test(expr_gene ~ grp)
    t1$diff_all <- unname(tt_all$estimate[2] - tt_all$estimate[1])  # Case - Control
    t1$p_all <- tt_all$p.value
    t1$diff_detected <- unname(tt$estimate[2] - tt$estimate[1])     # Case - Control
    t1$p_detected <- tt$p.value
    t1$n_case_detected <- sum(grp[detected] == "Case")
    t1$n_control_detected <- sum(grp[detected] == "Control")
    log_info("all samples : diff = %+.3f, p = %.3g", t1$diff_all, t1$p_all)
    log_info("detected only: diff = %+.3f, p = %.3g  (Case n=%d, Control n=%d)",
             t1$diff_detected, t1$p_detected,
             t1$n_case_detected, t1$n_control_detected)
  } else {
    t1$diff_all <- NA; t1$p_all <- NA
    t1$diff_detected <- NA; t1$p_detected <- NA
    t1$n_case_detected <- sum(grp[detected] == "Case")
    t1$n_control_detected <- sum(grp[detected] == "Control")
    log_warn("too few detected samples (or one group empty) for a re-test -- this itself is evidence the measurement is unreliable")
  }
  results_t1[[gse]] <- t1

  # ---- TEST 2: which cell compartment does the gene track? ---------------
  log_step("[%s] TEST 2 -- which marker panel does %s track?", gse, GENE)

  score_of <- function(symbols) {
    s <- marker_score(vsd, symbols)
    if (is.null(s)) return(NULL)
    as.numeric(s[colnames(vsd)])
  }
  # Every panel in MARKERS, so the gene's own compartment is always compared
  # against the competing one rather than tested in isolation.
  panel_scores <- lapply(MARKERS, score_of)

  for (nm in names(panel_scores)) {
    sc <- panel_scores[[nm]]
    if (is.null(sc)) next
    for (g in c(GENE, PARTNER)) {
      v <- if (g == GENE) expr_gene else expr_partner
      ct <- stats::cor.test(v, sc)
      results_t2[[length(results_t2) + 1]] <- data.frame(
        study = gse, gene = g, panel = nm, n = length(v),
        r = unname(ct$estimate),
        ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
        p_value = ct$p.value, stringsAsFactors = FALSE)
      log_info("%s vs %-12s markers: r = %+.3f [%.3f, %.3f], p = %.3g",
               g, nm, ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value)
    }
  }

  # ---- TEST 3: gene-partner coupling ------------------------------------
  log_step("[%s] TEST 3 -- %s-%s coupling%s", gse, GENE, PARTNER,
           if (is.na(ATLAS_R)) "" else sprintf(" (atlas benchmark r = %.3f)", ATLAS_R))
  ct <- stats::cor.test(expr_partner, expr_gene)
  eff <- read.csv(file.path(PATHS$de, paste0(gse, "_effect_for_meta.csv")))
  lfc <- setNames(eff$log2FC_shrunk, eff$symbol)
  padj <- setNames(eff$padj, eff$symbol)

  results_t3[[gse]] <- data.frame(
    study = gse, gene = GENE, partner = PARTNER, n = length(expr_gene),
    r_gene_partner = unname(ct$estimate),
    ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
    p_value = ct$p.value,
    atlas_benchmark_r = ATLAS_R,
    gene_log2FC = unname(lfc[GENE]), gene_padj = unname(padj[GENE]),
    partner_log2FC = unname(lfc[PARTNER]), partner_padj = unname(padj[PARTNER]),
    stringsAsFactors = FALSE)
  log_info("%s-%s r in %s = %+.3f [%.3f, %.3f], p = %.3g",
           GENE, PARTNER, TCLASS,
           ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value)
  log_info("%s log2FC = %+.3f (padj %.3g) vs %s log2FC = %+.3f (padj %.3g)",
           GENE, lfc[GENE], padj[GENE], PARTNER, lfc[PARTNER], padj[PARTNER])
}

# ---- save ------------------------------------------------------------------

out_dir <- PATHS$cellcomp
stem <- sprintf("%s_%s", GENE, TCLASS)
if (length(results_t1)) save_result(do.call(rbind, results_t1), out_dir,
                                    paste0(stem, "_TEST1_detection_restricted"))
if (length(results_t2)) save_result(do.call(rbind, results_t2), out_dir,
                                    paste0(stem, "_TEST2_marker_tracking"))
if (length(results_t3)) save_result(do.call(rbind, results_t3), out_dir,
                                    paste0(stem, "_TEST3_partner_coupling"))

# ---- verdict ---------------------------------------------------------------

log_step("VERDICT")
t3 <- do.call(rbind, results_t3)
t2 <- do.call(rbind, results_t2)

if (!is.null(t3)) {
  # If the two genes are tightly coupled at baseline, a change driven by cell
  # abundance has to move both. Count the studies where it moved only one.
  # 0.75 and 0.25 log2 units are a deliberately blunt "clearly moved" vs
  # "clearly did not" split -- the CSV carries the exact estimates.
  decoupled <- sum(abs(t3$partner_log2FC) < 0.25 & abs(t3$gene_log2FC) > 0.75,
                   na.rm = TRUE)
  log_info("studies where %s moves >0.75 log2 while %s stays flat (<0.25): %d/%d",
           GENE, PARTNER, decoupled, nrow(t3))
}
if (!is.null(t2)) {
  for (nm in unique(t2$panel)) {
    sub <- t2[t2$gene == GENE & t2$panel == nm, ]
    log_info("mean %s-%s marker r = %+.3f (k = %d)",
             GENE, nm, mean(sub$r, na.rm = TRUE), nrow(sub))
  }
}

log_session_info(SCRIPT)
log_step("COMPLETE: sensitivity tests")
