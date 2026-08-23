#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_run_dataset.R -- the reusable per-dataset pipeline. This is the workhorse.
#
# One dataset, end to end, to the standard in docs/METHODS.md:
#
#   fetch -> collapse technical replicates -> assign groups -> filter
#         -> DESeq2 + shrinkage (METHODS S1)
#         -> cell-composition quantification (METHODS S2: adjusted model,
#            correlation, xCell)
#         -> PCA + variance decomposition (METHODS S3)
#         -> figures (METHODS S6)
#         -> sessionInfo + run-log row (METHODS S5)
#
# Usage:
#   Rscript 01_scripts/10_run_dataset.R GSE112087
#   Rscript 01_scripts/10_run_dataset.R GSE112087 --group-col char_phenotype
#
# Saves everything under 04_processed/, 05_results/ and 06_figures/, and
# appends a row to 08_logs/RUN_LOG.md recording what it found.
#
# Re-runnable: cached GEO downloads in 03_raw_data/ are reused, so a re-run
# costs minutes rather than a fresh download of every supplementary file.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery); library(DESeq2); library(data.table)
  library(SummarizedExperiment); library(Biobase)
})

source("01_scripts/00_config.R")
source("01_scripts/00_theme.R")
for (f in c("utils", "annotate", "download", "prepare", "de", "deconv",
            "meta", "figures")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

# ---- arguments ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript 01_scripts/10_run_dataset.R <GSE> [--group-col <col>] [--skip-xcell]")
GSE <- args[1]
opt_group_col <- if ("--group-col" %in% args) args[which(args == "--group-col") + 1] else NULL
opt_skip_xcell <- "--skip-xcell" %in% args

SCRIPT <- paste0("10_run_dataset_", GSE)

# ---------------------------------------------------------------------------
run_dataset <- function(gse, group_col = NULL, skip_xcell = FALSE) {

  info <- ds_info(gse)
  results <- list(gse = gse, info = info)

  # ---- 1. fetch ----------------------------------------------------------
  fetched <- fetch_gse(gse)
  if (is.null(fetched$counts)) {
    append_progress(gse, "Fetch", "BLOCKED",
                    "No count matrix obtainable from GEO", "-")
    stop("No counts for ", gse)
  }
  counts <- fetched$counts
  meta   <- fetched$metadata

  # DESeq2 needs raw counts; normalized matrices (RPKM/FPKM/TPM) take the
  # limma-trend + ashr route instead. Both produce the same result schema, so
  # everything downstream (composition, figures, meta-analysis) is unchanged.
  de_engine <- if (fetched$scale$is_raw_counts) "deseq2" else "limma"
  if (de_engine == "limma") {
    log_warn("values are %s -- using limma-trend instead of DESeq2.",
             fetched$scale$verdict)
    append_progress(gse, "Value scale check", "DONE",
                    sprintf("Matrix is %s, not raw counts -> limma-trend + ashr path selected (DESeq2 not valid on normalized data).",
                            fetched$scale$verdict),
                    file.path("04_processed", paste0(gse, "_counts.csv")))
  }
  results$de_engine <- de_engine

  # ---- 2. collapse technical replicates -----------------------------------
  log_step("technical replicate check")
  rep_key <- detect_replicate_key(meta)
  if (!is.null(rep_key)) {
    col <- collapse_technical_replicates(counts, meta, rep_key)
    counts <- col$counts
    meta   <- col$metadata
    path <- file.path(PATHS$processed, paste0(gse, "_counts_collapsed.csv"))
    data.table::fwrite(cbind(gene_id = rownames(counts), as.data.frame(counts)), path)
    rep_type <- attr(rep_key, "type") %||% "technical"
    varying <- attr(rep_key, "varying_fields")
    append_progress(gse,
                    sprintf("Collapse replicates (%s)", rep_type), "DONE",
                    sprintf("Collapsed on '%s': %d GSMs -> %d donor-level samples (counts summed). Type: %s.%s",
                            as.character(rep_key), nrow(fetched$metadata), ncol(counts),
                            rep_type,
                            if (length(varying))
                              sprintf(" Within-donor samples differed on: %s -- distinct specimens aggregated to donor level; within-donor contrast discarded.",
                                      paste(varying, collapse = ", ")) else ""),
                    file.path("04_processed", basename(path)))
    results$replicate_type <- rep_type
    results$replicate_varying <- varying
  } else {
    meta$sample_id <- meta$gsm_id
    log_info("no technical replicates detected; %d samples", ncol(counts))
  }
  if (is.null(meta$sample_id)) meta$sample_id <- colnames(counts)

  # Align counts and metadata by shared IDs. Direct indexing fails whenever the
  # count columns could not all be mapped to GSM IDs, or a series matrix lists
  # samples that the supplementary matrix does not contain.
  common <- intersect(meta$sample_id, colnames(counts))
  if (length(common) < 2) {
    log_warn("counts and metadata share only %d sample ID(s)", length(common))
    log_warn("  count columns: %s", paste(head(colnames(counts), 4), collapse = ", "))
    log_warn("  metadata IDs : %s", paste(head(meta$sample_id, 4), collapse = ", "))
    append_progress(gse, "Sample alignment", "BLOCKED",
                    sprintf("Count columns could not be matched to GEO samples (%d shared IDs). Count columns look like: %s. Needs a manual mapping rule.",
                            length(common),
                            paste(head(colnames(counts), 3), collapse = ", ")),
                    file.path("04_processed", paste0(gse, "_counts.csv")))
    results$counts <- counts; results$metadata <- meta
    results$blocked <- "sample alignment"
    return(results)
  }
  if (length(common) < nrow(meta) || length(common) < ncol(counts)) {
    log_warn("aligning on %d shared samples (metadata had %d, counts had %d)",
             length(common), nrow(meta), ncol(counts))
  }
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]
  counts <- counts[, common, drop = FALSE]

  # ---- 3. groups ----------------------------------------------------------
  if (!info$de_capable) {
    log_warn("%s is registered as not DE-capable (%s) -- stopping after QC.",
             gse, info$use)
    append_progress(gse, "DE contrast", "SKIPPED",
                    "Registry marks this dataset as lacking a case/control contrast.", "-")
    results$counts <- counts; results$metadata <- meta
    return(results)
  }

  # An explicit contrast overrides automatic detection wherever the phenotype
  # column has more than two levels (see CONTRAST_SPEC in 00_config.R).
  spec <- CONTRAST_SPEC[[gse]]
  if (!is.null(spec)) {
    if (is.null(meta[[spec$col]])) {
      log_warn("CONTRAST_SPEC column '%s' absent from %s metadata", spec$col, gse)
    } else {
      v <- trimws(as.character(meta[[spec$col]]))
      keep_lv <- v %in% c(spec$case, spec$control)
      n_drop <- sum(!keep_lv)
      log_info("explicit contrast for %s: '%s' vs '%s' on column '%s'",
               gse, spec$case, spec$control, spec$col)
      if (n_drop > 0) {
        log_warn("dropping %d sample(s) in neither arm (levels: %s)", n_drop,
                 paste(sort(unique(v[!keep_lv])), collapse = ", "))
      }
      meta <- meta[keep_lv, , drop = FALSE]
      counts <- counts[, meta$sample_id, drop = FALSE]
      meta$group <- factor(ifelse(trimws(as.character(meta[[spec$col]])) == spec$control,
                                  "Control", "Case"),
                           levels = c("Control", "Case"))
      group_col <- spec$col
      results$contrast_spec <- sprintf("%s vs %s (%s); %d sample(s) excluded",
                                       spec$case, spec$control, spec$col, n_drop)
      tb <- table(meta$group)
      log_info("group assignment (explicit): %s",
               paste(sprintf("%s=%d", names(tb), tb), collapse = ", "))
    }
  }

  if (is.null(group_col)) {
    # exclude any field that varied within a donor -- it cannot be phenotype
    group_col <- pick_group_column(meta, exclude = results$replicate_varying %||% character(0))
  }
  if (is.null(group_col)) {
    log_warn("could not identify a phenotype column; pass --group-col")
    append_progress(gse, "Group assignment", "BLOCKED",
                    paste("No phenotype column found. Available:",
                          paste(grep("^char_", colnames(meta), value = TRUE),
                                collapse = ", ")), "-")
    results$metadata <- meta
    return(results)
  }
  if (is.null(meta$group)) meta$group <- assign_groups(meta, group_col)
  results$group_col <- group_col

  # ---- 4. differential expression (METHODS S1) ----------------------------------
  # include a batch term when one exists and is balanced across the contrast
  batch_col <- select_batch_covariate(meta, "group")
  if (!is.null(batch_col)) {
    meta$batch <- factor(as.character(meta[[batch_col]]))
    de_design <- ~ batch + group
  } else {
    de_design <- ~ group
  }
  results$batch_col <- batch_col

  if (de_engine == "deseq2") {
    dds <- build_dds(counts, meta, design = de_design)
    de  <- run_deseq(dds)
  } else {
    # normalized input: work on log2(x+1) throughout
    logmat <- log2(as.matrix(counts) + 1)
    de <- run_limma(logmat, meta, design_formula = de_design)
    dds <- NULL
  }
  save_result(de, PATHS$de, paste0(gse, "_DE_full"))

  panel_type <- info$tissue_class
  markers <- MARKERS[[panel_type]]
  goi <- extract_genes(de, GENES_OF_INTEREST)
  mrk <- extract_genes(de, markers)
  save_result(rbind(goi, mrk), PATHS$de, paste0(gse, "_DE_genes_of_interest"))

  results$de <- de; results$goi <- goi; results$markers <- mrk
  # Store the per-sample grouping so downstream scripts never have to re-derive
  # it. Re-deriving from metadata breaks for replicate-collapsed datasets, whose
  # sample IDs are donor-level and no longer match gsm_id (e.g. GSE134692).
  results$group <- setNames(as.character(meta$group), meta$sample_id)

  # record expression of every gene of interest, tested or not (METHODS S7 limitations)
  expr_for_detect <- if (de_engine == "deseq2") {
    t(t(as.matrix(counts)) / (colSums(as.matrix(counts)) / 1e6))   # CPM
  } else as.matrix(counts)
  det <- detection_report(expr_for_detect, de, GENES_OF_INTEREST)
  det$study <- gse; det$unit <- if (de_engine == "deseq2") "CPM" else "RPKM/FPKM/TPM"
  save_result(det, PATHS$de, paste0(gse, "_gene_detection_report"))
  results$detection <- det

  append_progress(gse, "DESeq2 + LFC shrinkage", "DONE",
                  sprintf("%s genes tested, %s significant at BH-FDR<%.2f (shrinkage: %s). %s",
                          fmt_n(nrow(de)), fmt_n(attr(de, "n_sig")), ALPHA,
                          attr(de, "shrink_type"),
                          paste(sprintf("%s log2FC=%.3f padj=%.3g", goi$symbol,
                                        goi$log2FC_shrunk, goi$padj), collapse = "; ")),
                  file.path("05_results/deseq2", paste0(gse, "_DE_full.csv")))

  # ---- 5. expression matrix + composition (METHODS S2) --------------------------
  # VST for counts; the already-normalized log matrix otherwise. Both are on a
  # log2-like scale, so the composition and figure code is identical.
  vsd <- if (de_engine == "deseq2") vst_matrix(dds, blind = TRUE) else
    log2(as.matrix(counts) + 1)
  results$vst <- vsd

  adj <- if (de_engine == "deseq2") {
    adjust_for_composition(dds, meta, vsd, markers)
  } else {
    adjust_for_composition_limma(vsd, meta, markers)
  }
  if (!is.null(adj)) save_result(adj, PATHS$cellcomp,
                                 paste0(gse, "_composition_adjusted_model"))
  corr <- correlate_with_composition(vsd, meta, markers, de_res = de)
  if (!is.null(corr)) save_result(corr, PATHS$cellcomp,
                                  paste0(gse, "_gene_vs_markerscore_correlation"))
  results$adjusted <- adj; results$marker_corr <- corr

  verdict <- composition_verdict(goi, mrk, adj)
  log_info("composition verdict -> %s", verdict)
  results$verdict <- verdict

  append_progress(gse, "Cell-composition quantification", "DONE",
                  sprintf("Marker panel: %s. %s", paste(markers, collapse = "/"), verdict),
                  file.path("05_results/cell_composition",
                            paste0(gse, "_composition_adjusted_model.csv")))

  # ---- 6. xCell (METHODS S2.4) --------------------------------------------------
  if (!skip_xcell) {
    xc <- tryCatch(run_xcell(counts, gse, linear_input = (de_engine == "limma")), error = function(e) {
      log_warn("xCell failed: %s", conditionMessage(e)); NULL })
    if (!is.null(xc)) {
      cmp <- compare_xcell_groups(xc, meta)
      if (!is.null(cmp)) save_result(cmp, PATHS$cellcomp,
                                     paste0(gse, "_xcell_group_comparison"))
      xcorr <- correlate_gene_with_xcell(vsd, xc, tissue_class = info$tissue_class)
      if (!is.null(xcorr)) save_result(xcorr, PATHS$cellcomp,
                                       paste0(gse, "_gene_vs_xcell_correlation"))
      results$xcell <- xc; results$xcell_cmp <- cmp
      append_progress(gse, "xCell deconvolution", "DONE",
                      sprintf("%d cell types scored; %d differ between groups at BH-FDR<%.2f.",
                              nrow(xc), if (is.null(cmp)) 0 else sum(cmp$p_adj < ALPHA, na.rm = TRUE), ALPHA),
                      file.path("05_results/cell_composition",
                                paste0(gse, "_xcell_scores.csv")))
    }
  }

  # ---- 7. variance structure (METHODS S3) --------------------------------------
  # always include group; add demographics plus any batch-like field the series
  # records (GSE202625 has char_sequencing.batch across 5 runs), since Quality
  # Bar 4.3 asks specifically for batch vs biological-group variance
  batch_cols <- grep("batch|run|lane|flowcell|site|centre|center|platform|processingdate",
                     grep("^char_", colnames(meta), value = TRUE),
                     value = TRUE, ignore.case = TRUE)
  covars <- unique(c("group",
                     intersect(c("char_gender", "char_Sex", "char_sex",
                                 "char_ancestry", "char_age",
                                 "n_technical_reps"), colnames(meta)),
                     batch_cols))
  covars <- intersect(covars, colnames(meta))
  if (length(batch_cols)) log_info("batch covariates detected: %s",
                                   paste(batch_cols, collapse = ", "))
  pcv <- pc_variance_explained(vsd, meta, covars)
  if (!is.null(pcv)) {
    save_result(pcv, PATHS$results, paste0(gse, "_pc_variance_by_covariate"))
    vsum <- summarise_variance(pcv)
    save_result(vsum, PATHS$results, paste0(gse, "_variance_summary"))
    results$variance <- vsum
    log_info("variance explained: %s",
             paste(sprintf("%s=%.1f%%", vsum$covariate,
                           vsum$total_var_explained_pct), collapse = ", "))
  }

  # ---- 8. figures (METHODS S6) --------------------------------------------------
  log_step("figures")
  fig_prefix <- paste0(gse, "_")
  try({
    p <- fig_volcano(de, gse, GENES_OF_INTEREST, markers)
    save_fig(p, paste0(fig_prefix, "volcano"), width = 7.5, height = 5.5)
  }, silent = FALSE)
  try({
    p <- fig_boxplot(vsd, meta, gse, GENES_OF_INTEREST, "group", goi)
    if (!is.null(p)) save_fig(p, paste0(fig_prefix, "boxplot_GOI"), width = 6.5, height = 4)
  }, silent = FALSE)
  try({
    p <- fig_heatmap(vsd, meta, gse)
    if (!is.null(p)) save_fig(p, paste0(fig_prefix, "heatmap_markers"), width = 8, height = 4.5)
  }, silent = FALSE)
  try({
    p <- fig_pca(vsd, meta, gse, "group")
    save_fig(p, paste0(fig_prefix, "pca"), width = 6.5, height = 5)
  }, silent = FALSE)
  try({
    if (!is.null(results$variance)) {
      p <- fig_variance_explained(results$variance, gse)
      save_fig(p, paste0(fig_prefix, "variance_explained"), width = 6.5, height = 4)
    }
  }, silent = FALSE)
  try({
    p <- fig_correlation(vsd, meta, gse, GENES_OF_INTEREST[1],
                         GENES_OF_INTEREST[2], "group")
    if (!is.null(p)) save_fig(p, paste0(fig_prefix, "gene_pair_correlation"), width = 6, height = 5)
  }, silent = FALSE)

  append_progress(gse, "Figures", "DONE",
                  "Volcano, boxplot, marker heatmap, PCA, variance, gene-pair correlation (PNG + PDF).",
                  "06_figures/png/, 06_figures/pdf/")

  # ---- 9. per-study row for the meta-analysis ----------------------------
  eff <- goi
  eff$study <- gse
  eff$disease <- info$disease
  eff$n_samples <- ncol(counts)
  eff$n_case <- sum(meta$group == "Case")
  eff$n_control <- sum(meta$group == "Control")
  save_result(eff, PATHS$de, paste0(gse, "_effect_for_meta"))
  results$effect <- eff

  results
}

# ---------------------------------------------------------------------------

res <- run_dataset(GSE, opt_group_col, opt_skip_xcell)
saveRDS(res, file.path(PATHS$processed, paste0(GSE, "_pipeline_result.rds")))
log_info("pipeline object -> %s",
         file.path(PATHS$processed, paste0(GSE, "_pipeline_result.rds")))

log_session_info(SCRIPT)
log_step("COMPLETE: %s", GSE)
