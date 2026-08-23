# ---------------------------------------------------------------------------
# R/deconv.R -- cell-type deconvolution (METHODS S2.4).
#
# xCell enrichment scores for every bulk dataset, saved as supplementary
# results. This is the quantitative counterpart to the marker panel: the panel
# says "endothelial genes moved", xCell says "the endothelial fraction moved by
# this much, and the difference between groups is significant at p = ...".
#
# xCell expects a SYMBOL-indexed expression matrix on a roughly linear scale
# (TPM/RPKM/normalized counts), not raw counts and not log-transformed.
# ---------------------------------------------------------------------------

# Cell types most relevant to this project. xCell returns 64 + 3 scores; we
# save all of them but summarise on these.
XCELL_FOCUS <- c("Endothelial cells", "mv Endothelial cells",
                 "ly Endothelial cells", "Fibroblasts",
                 "Neutrophils", "Monocytes", "Macrophages",
                 "CD4+ T-cells", "CD8+ T-cells", "B-cells", "NK cells",
                 "Platelets", "Epithelial cells",
                 "ImmuneScore", "StromaScore", "MicroenvironmentScore")

# ---- input preparation ----------------------------------------------------

# Counts -> CPM on a linear scale, indexed by symbol. xCell scores are rank
# based, so CPM is an appropriate input; log transformation would distort ranks
# less but xCell's signatures were built on linear-scale data.
prepare_xcell_input <- function(counts, linear_input = FALSE) {
  mat <- as.matrix(counts)
  storage.mode(mat) <- "double"
  if (linear_input) {
    # already RPKM/FPKM/TPM -- on a linear scale, so use as-is
    log_info("xCell input is already normalized (linear scale); no CPM step")
    cpm <- mat
  } else {
    lib <- colSums(mat, na.rm = TRUE)
    cpm <- sweep(mat, 2, lib, "/") * 1e6
  }
  cpm <- matrix_to_symbols(cpm)
  # xCell needs reasonable gene coverage to be meaningful
  log_info("xCell input: %s symbols x %d samples", fmt_n(nrow(cpm)), ncol(cpm))
  if (nrow(cpm) < 5000) {
    log_warn("only %s genes after symbol mapping -- xCell scores may be unstable",
             fmt_n(nrow(cpm)))
  }
  cpm
}

# ---- scoring --------------------------------------------------------------

run_xcell <- function(counts, gse, save = TRUE, linear_input = FALSE) {
  require_pkg("xCell", "cell-type deconvolution (METHODS S2.4)")
  log_step("xCell deconvolution: %s", gse)

  # xCell ships its signature matrices as a lazy-loaded dataset (xCell.data).
  # Calling xCell::xCellAnalysis() without attaching the package leaves that
  # object unreachable ("object 'xCell.data' not found"), so attach it here.
  if (!"package:xCell" %in% search()) {
    suppressPackageStartupMessages(library(xCell))
  }
  # If attaching alone did not expose it, load it into the global environment.
  # Package functions resolve through namespace -> imports -> base -> globalenv,
  # so this is reachable from inside xCellAnalysis().
  if (!exists("xCell.data")) {
    utils::data("xCell.data", package = "xCell", envir = globalenv())
  }

  input <- prepare_xcell_input(counts, linear_input = linear_input)
  set.seed(SEED)   # xCell's significance routine samples randomly

  scores <- xCell::xCellAnalysis(input, rnaseq = TRUE, parallel.sz = NCPUS)
  scores <- as.data.frame(scores)
  log_info("xCell returned %d cell-type scores x %d samples",
           nrow(scores), ncol(scores))

  if (save) {
    path <- file.path(PATHS$cellcomp, paste0(gse, "_xcell_scores.csv"))
    utils::write.csv(cbind(cell_type = rownames(scores), scores),
                     path, row.names = FALSE)
    log_info("xCell scores -> %s", path)
  }
  scores
}

# ---- group comparison -----------------------------------------------------

# Tests whether each cell-type score differs between Case and Control. This is
# what makes the deconvolution a *result* rather than a decoration: it states
# whether composition actually shifted, with a statistic and BH-adjusted p.
compare_xcell_groups <- function(scores, meta, group_col = "group",
                                 sample_col = "sample_id") {
  ids <- meta[[sample_col]]
  common <- intersect(colnames(scores), ids)
  if (length(common) < 4) {
    log_warn("too few matched samples (%d) for xCell group comparison",
             length(common))
    return(NULL)
  }
  grp <- droplevels(factor(meta[[group_col]][match(common, ids)]))
  if (nlevels(grp) < 2) {
    log_warn("only one group -- skipping xCell group comparison")
    return(NULL)
  }
  s <- scores[, common, drop = FALSE]

  rows <- lapply(rownames(s), function(ct) {
    v <- as.numeric(s[ct, ])
    if (stats::sd(v, na.rm = TRUE) == 0) return(NULL)
    wt <- suppressWarnings(stats::wilcox.test(v ~ grp))
    means <- tapply(v, grp, mean, na.rm = TRUE)
    data.frame(
      cell_type = ct,
      n_control = sum(grp == levels(grp)[1]),
      n_case    = sum(grp == levels(grp)[2]),
      mean_control = unname(means[1]),
      mean_case    = unname(means[2]),
      difference   = unname(means[2] - means[1]),
      W = unname(wt$statistic),
      p_value = wt$p.value,
      test = "Wilcoxon rank-sum",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")   # METHODS S1.1
  out <- out[order(out$p_adj, out$p_value), ]
  log_info("xCell group comparison: %d/%d cell types differ at BH-FDR < %.2f",
           sum(out$p_adj < ALPHA, na.rm = TRUE), nrow(out), ALPHA)
  out
}

# ---- correlation of the gene of interest with a deconvolved cell fraction --

# The most direct test of the composition question: does the gene track the
# estimated endothelial fraction across samples?
# Cell types worth correlating against, by tissue context. xCell scores are
# relative enrichment across a mixed-tissue reference, so a type that is absent
# (endothelium in whole blood) or saturated (neutrophils in whole blood) gets
# floored at zero and is not a usable continuous variable.
xcell_types_for <- function(tissue_class) {
  if (identical(tissue_class, "blood")) {
    c("Neutrophils", "Monocytes", "CD4+ T-cells", "CD8+ T-cells", "B-cells",
      "NK cells", "Platelets", "ImmuneScore")
  } else {
    c("Endothelial cells", "mv Endothelial cells", "ly Endothelial cells",
      "Fibroblasts", "Epithelial cells", "StromaScore", "ImmuneScore")
  }
}

# A score with many tied zeros or near-zero spread cannot support a correlation.
# Flag rather than silently produce an impressive-looking r.
.score_quality <- function(v) {
  list(frac_zero = mean(v == 0, na.rm = TRUE),
       sd = stats::sd(v, na.rm = TRUE))
}

correlate_gene_with_xcell <- function(mat_vst, scores, genes = GENES_OF_INTEREST,
                                      cell_types = NULL, tissue_class = NULL,
                                      method = "pearson",
                                      max_frac_zero = 0.25, min_sd = 1e-3) {
  if (is.null(cell_types)) cell_types <- xcell_types_for(tissue_class)
  panel <- resolve_panel(genes, rownames(mat_vst))
  common <- intersect(colnames(mat_vst), colnames(scores))
  if (length(common) < 4) return(NULL)

  rows <- list()
  for (i in seq_len(nrow(panel))) {
    if (!panel$found[i]) next
    expr <- as.numeric(mat_vst[panel$id[i], common])
    for (ct in intersect(cell_types, rownames(scores))) {
      cv <- as.numeric(scores[ct, common])
      q <- .score_quality(cv)
      if (q$sd == 0) next
      usable <- q$frac_zero <= max_frac_zero && q$sd >= min_sd
      if (!usable && i == 1) {
        log_warn("xCell '%s' unreliable here (%.0f%% exact zeros, sd=%.5f) -- reported but flagged",
                 ct, 100 * q$frac_zero, q$sd)
      }
      tt <- stats::cor.test(expr, cv, method = method)
      rows[[length(rows) + 1]] <- data.frame(
        gene = panel$symbol[i], cell_type = ct, n = length(common),
        r = unname(tt$estimate),
        ci_low  = if (!is.null(tt$conf.int)) tt$conf.int[1] else NA_real_,
        ci_high = if (!is.null(tt$conf.int)) tt$conf.int[2] else NA_real_,
        p_value = tt$p.value, method = method,
        frac_zero_score = round(q$frac_zero, 3),
        score_sd = signif(q$sd, 3),
        reliable = usable,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(NULL)
  # adjust only across the reliable comparisons; degenerate scores should not
  # consume multiple-testing budget or be read as findings
  out$p_adj <- NA_real_
  if (any(out$reliable)) {
    out$p_adj[out$reliable] <- stats::p.adjust(out$p_value[out$reliable],
                                               method = "BH")
  }
  out[order(-out$reliable, out$p_value), ]
}
