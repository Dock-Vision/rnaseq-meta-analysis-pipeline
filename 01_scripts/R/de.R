# ---------------------------------------------------------------------------
# R/de.R -- differential expression to the analysis standard standard.
#
#   METHODS S1.1  report raw p AND BH-adjusted padj, method stated
#   METHODS S1.2  shrunken log2FoldChange (apeglm, ashr fallback) -- never raw LFC
#   METHODS S2.3  cell-composition correction QUANTIFIED, not eyeballed:
#             (a) marker-adjusted model re-fit
#             (b) correlation between the gene's change and the marker change,
#                 with a test statistic, CI and p-value
# ---------------------------------------------------------------------------

# ---- core DE --------------------------------------------------------------

run_deseq <- function(dds, contrast_name = "group_Case_vs_Control",
                      shrink_type = LFC_SHRINK, alpha = ALPHA) {
  require_pkg("DESeq2", "differential expression")
  log_step("DESeq2: %s", contrast_name)

  set.seed(SEED)   # METHODS S5
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # unshrunken results carry the p-values; BH adjustment is DESeq2's default
  # but we state it explicitly so the method is never ambiguous
  res <- DESeq2::results(dds, name = contrast_name, alpha = alpha,
                         pAdjustMethod = "BH")

  # METHODS S1.2 -- shrink effect sizes. apeglm needs the contrast to be a model
  # coefficient; fall back to ashr when it is not.
  shrunk <- tryCatch({
    require_pkg(shrink_type, "LFC shrinkage")
    DESeq2::lfcShrink(dds, coef = contrast_name, type = shrink_type, quiet = TRUE)
  }, error = function(e) {
    log_warn("%s shrinkage failed (%s); falling back to ashr",
             shrink_type, conditionMessage(e))
    require_pkg("ashr", "LFC shrinkage fallback")
    shrink_type <<- "ashr"
    DESeq2::lfcShrink(dds, coef = contrast_name, type = "ashr", quiet = TRUE)
  })

  out <- data.frame(
    gene_id        = rownames(res),
    baseMean       = res$baseMean,
    log2FC_raw     = res$log2FoldChange,      # kept for transparency only
    log2FC_shrunk  = shrunk$log2FoldChange,   # <- the reported effect size
    lfcSE_shrunk   = shrunk$lfcSE,
    stat           = res$stat,
    pvalue         = res$pvalue,              # METHODS S1.1 raw p
    padj           = res$padj,                # METHODS S1.1 BH-adjusted
    stringsAsFactors = FALSE
  )
  out$symbol <- ids_to_symbols(out$gene_id)
  out <- out[order(out$padj, out$pvalue), ]

  attr(out, "shrink_type") <- shrink_type
  attr(out, "alpha") <- alpha
  attr(out, "n_sig") <- sum(out$padj < alpha, na.rm = TRUE)
  attr(out, "dds") <- dds

  log_info("%s genes tested; %s significant at BH-FDR < %.2f (shrinkage: %s)",
           fmt_n(nrow(out)), fmt_n(attr(out, "n_sig")), alpha, shrink_type)
  out
}

# ---- limma path for normalized (non-count) matrices -----------------------

# Several GEO series ship only RPKM/FPKM/TPM (e.g. GSE72509). DESeq2 is invalid
# on those, so we use limma-trend, which is the standard approach for
# log-transformed continuous expression.
#
# Note on METHODS S1.2: eBayes moderates the *variance*, not the effect
# size, so limma alone does not satisfy the shrinkage requirement. We therefore
# pass limma's coefficients and standard errors through ashr, which shrinks
# effect sizes given (betahat, sebetahat). The reported log2FC_shrunk is
# comparable in spirit to apeglm output and is safe to pool in the
# meta-analysis.
run_limma <- function(mat_log, meta, coef_name = "groupCase", alpha = ALPHA,
                      min_expr = 1, min_prop = 0.25,
                      design_formula = ~ group) {
  require_pkg("limma", "differential expression on normalized data")
  log_step("limma-trend (normalized input): %s", coef_name)

  mat <- as.matrix(mat_log)
  keep <- rowMeans(mat > min_expr) >= min_prop
  log_info("expression filter: %s -> %s genes (>%g in >=%.0f%% of samples)",
           fmt_n(nrow(mat)), fmt_n(sum(keep)), min_expr, 100 * min_prop)
  mat <- mat[keep, , drop = FALSE]

  design <- stats::model.matrix(design_formula, data = meta)
  set.seed(SEED)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  if (!coef_name %in% colnames(fit$coefficients)) {
    coef_name <- colnames(fit$coefficients)[2]
    log_info("using coefficient '%s'", coef_name)
  }
  tt <- limma::topTable(fit, coef = coef_name, number = Inf,
                        adjust.method = "BH", sort.by = "none")

  # moderated standard error of the coefficient
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name]

  lfc_shrunk <- tt$logFC
  shrink_type <- "none"
  ok <- tryCatch({
    require_pkg("ashr", "effect-size shrinkage for the limma path")
    a <- ashr::ash(betahat = tt$logFC, sebetahat = as.numeric(se),
                   mixcompdist = "normal")
    lfc_shrunk <- ashr::get_pm(a)   # posterior mean effect size
    shrink_type <- "ashr (on limma coefficients)"
    TRUE
  }, error = function(e) {
    log_warn("ashr shrinkage failed (%s) -- reporting unshrunk limma logFC",
             conditionMessage(e)); FALSE
  })
  if (!ok) shrink_type <- "none (ashr unavailable)"

  out <- data.frame(
    gene_id       = rownames(tt),
    baseMean      = tt$AveExpr,
    log2FC_raw    = tt$logFC,
    log2FC_shrunk = lfc_shrunk,
    lfcSE_shrunk  = as.numeric(se),
    stat          = tt$t,
    pvalue        = tt$P.Value,
    padj          = tt$adj.P.Val,
    stringsAsFactors = FALSE)
  out$symbol <- ids_to_symbols(out$gene_id)
  out <- out[order(out$padj, out$pvalue), ]

  attr(out, "shrink_type") <- shrink_type
  attr(out, "method") <- "limma-trend (eBayes, robust), BH FDR"
  attr(out, "alpha") <- alpha
  attr(out, "n_sig") <- sum(out$padj < alpha, na.rm = TRUE)
  attr(out, "fit") <- fit

  log_info("%s genes tested; %s significant at BH-FDR < %.2f (shrinkage: %s)",
           fmt_n(nrow(out)), fmt_n(attr(out, "n_sig")), alpha, shrink_type)
  out
}

# Composition adjustment for the limma path: refit with the marker score as a
# covariate, mirroring adjust_for_composition() for DESeq2.
adjust_for_composition_limma <- function(mat_log, meta, panel_symbols,
                                         genes = GENES_OF_INTEREST,
                                         alpha = ALPHA,
                                         min_expr = 1, min_prop = 0.25) {
  require_pkg("limma", "composition-adjusted model")
  log_step("composition-adjusted model (limma path)")

  score <- marker_score(mat_log, panel_symbols)
  if (is.null(score)) return(NULL)

  mat <- as.matrix(mat_log)
  keep <- rowMeans(mat > min_expr) >= min_prop
  mat <- mat[keep, , drop = FALSE]

  meta$marker_score <- as.numeric(score[colnames(mat)])
  design <- stats::model.matrix(~ marker_score + group, data = meta)
  set.seed(SEED)
  fit <- limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  cf <- grep("^group", colnames(fit$coefficients), value = TRUE)[1]
  tt <- limma::topTable(fit, coef = cf, number = Inf, adjust.method = "BH",
                        sort.by = "none")

  tab <- data.frame(gene_id = rownames(tt),
                    log2FC_adjusted = tt$logFC,
                    pvalue_adjusted_model = tt$P.Value,
                    padj_adjusted_model = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  tab$symbol <- ids_to_symbols(tab$gene_id)
  out <- tab[!is.na(tab$symbol) & tab$symbol %in% genes, ]
  attr(out, "marker_genes") <- attr(score, "genes_used")
  attr(out, "marker_score") <- score
  out
}

# ---- genes of interest ----------------------------------------------------

# A gene of interest can vanish from a DE table for two very different reasons:
# it is genuinely absent from the platform, or it is present but too lowly
# expressed to survive the filter (e.g. an endothelial gene in whole blood).
# Those mean opposite things, so record the expression level of every gene of
# interest regardless of whether it was tested.
detection_report <- function(expr_mat, de_res, symbols = GENES_OF_INTEREST,
                             detect_threshold = 1) {
  panel <- resolve_panel(symbols, rownames(expr_mat))
  rows <- lapply(seq_len(nrow(panel)), function(i) {
    if (!panel$found[i]) {
      return(data.frame(symbol = panel$symbol[i], in_matrix = FALSE,
                        mean_expr = NA_real_, median_expr = NA_real_,
                        frac_detected = NA_real_, tested_in_DE = FALSE,
                        note = "not present on platform / unmapped",
                        stringsAsFactors = FALSE))
    }
    v <- as.numeric(expr_mat[panel$id[i], ])
    tested <- panel$id[i] %in% de_res$gene_id ||
      panel$symbol[i] %in% de_res$symbol
    data.frame(
      symbol = panel$symbol[i], in_matrix = TRUE,
      mean_expr = mean(v, na.rm = TRUE),
      median_expr = stats::median(v, na.rm = TRUE),
      frac_detected = mean(v > detect_threshold, na.rm = TRUE),
      tested_in_DE = tested,
      note = if (tested) "tested"
             else "present but below expression filter -- NOT testable here",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  for (i in seq_len(nrow(out))) {
    if (out$in_matrix[i] && !out$tested_in_DE[i]) {
      log_warn("%s present but not testable: median expression %.3f, detected in %.0f%% of samples",
               out$symbol[i], out$median_expr[i], 100 * out$frac_detected[i])
    }
  }
  out
}

extract_genes <- function(de_res, symbols = c(GENES_OF_INTEREST,
                                              unlist(MARKERS, use.names = FALSE))) {
  sub <- de_res[!is.na(de_res$symbol) & de_res$symbol %in% symbols, , drop = FALSE]
  missing <- setdiff(symbols, sub$symbol)
  if (length(missing)) {
    log_warn("not in DE table (filtered or unmapped): %s",
             paste(missing, collapse = ", "))
  }
  sub[order(match(sub$symbol, symbols)), ]
}

# ---- METHODS S2.3(a) marker-adjusted model --------------------------------------

# Builds a per-sample marker score (mean VST expression of the panel), then
# re-fits DE with that score as a covariate. If the gene's effect survives, the
# signal is not merely a composition shift.
marker_score <- function(mat_vst, panel_symbols, method = c("pc1", "mean"),
                         id_type_rownames = NULL) {
  method <- match.arg(method)
  ids <- resolve_panel(panel_symbols, rownames(mat_vst))
  ids <- ids[ids$found, ]
  if (nrow(ids) == 0) {
    log_warn("no marker genes found in matrix -- marker score unavailable")
    return(NULL)
  }
  sub <- mat_vst[ids$id, , drop = FALSE]
  # scale each marker so one high-expression gene (e.g. VWF, PTPRC) does not
  # dominate purely by magnitude
  z <- t(scale(t(sub)))

  if (method == "mean" || nrow(z) < 2) {
    # simple average: only valid when the panel moves as one block, since
    # markers with opposite behaviour cancel each other out
    score <- colMeans(z, na.rm = TRUE)
    attr(score, "method") <- "mean of z-scores"
  } else {
    # PC1 of the marker panel. Signed loadings mean a panel containing both
    # rising and falling markers (e.g. FCGR3B up, CD3E down in SLE blood)
    # still yields one coherent composition axis instead of cancelling.
    rownames(z) <- ids$symbol   # label the axis by gene symbol, not Ensembl ID
    pc <- stats::prcomp(t(z), center = TRUE, scale. = FALSE)
    score <- pc$x[, 1]
    # orient the axis so it points the same way as the panel's overall level,
    # keeping the sign of the score interpretable across datasets
    if (stats::cor(score, colMeans(z, na.rm = TRUE)) < 0) score <- -score
    var_pct <- 100 * pc$sdev[1]^2 / sum(pc$sdev^2)
    loadings <- pc$rotation[, 1]
    if (stats::cor(pc$x[, 1], colMeans(z, na.rm = TRUE)) < 0) loadings <- -loadings
    attr(score, "method") <- sprintf("PC1 of marker panel (%.1f%% of panel variance)",
                                     var_pct)
    attr(score, "var_pct") <- var_pct
    attr(score, "loadings") <- loadings
    log_info("marker PC1 explains %.1f%% of panel variance; loadings: %s",
             var_pct, paste(sprintf("%s=%+.2f", names(loadings), loadings),
                            collapse = " "))
  }
  names(score) <- colnames(mat_vst)
  log_info("marker score (%s) from %d/%d genes: %s",
           attr(score, "method"), nrow(ids), length(panel_symbols),
           paste(ids$symbol, collapse = ", "))
  attr(score, "genes_used") <- ids$symbol
  score
}

# Re-fit with ~ marker_score + group and report how the gene of interest moves.
adjust_for_composition <- function(dds, meta, mat_vst, panel_symbols,
                                   genes = GENES_OF_INTEREST,
                                   alpha = ALPHA) {
  log_step("composition-adjusted model (METHODS S2)")
  score <- marker_score(mat_vst, panel_symbols)
  if (is.null(score)) return(NULL)

  cd <- SummarizedExperiment::colData(dds)
  cd$marker_score <- as.numeric(score[colnames(dds)])
  SummarizedExperiment::colData(dds) <- cd
  DESeq2::design(dds) <- ~ marker_score + group

  set.seed(SEED)
  dds2 <- DESeq2::DESeq(dds, quiet = TRUE)
  res2 <- DESeq2::results(dds2, name = "group_Case_vs_Control",
                          alpha = alpha, pAdjustMethod = "BH")
  shr2 <- tryCatch(
    DESeq2::lfcShrink(dds2, coef = "group_Case_vs_Control",
                      type = LFC_SHRINK, quiet = TRUE),
    error = function(e) DESeq2::lfcShrink(dds2, coef = "group_Case_vs_Control",
                                          type = "ashr", quiet = TRUE))

  tab <- data.frame(
    gene_id = rownames(res2),
    log2FC_adjusted = shr2$log2FoldChange,
    pvalue_adjusted_model = res2$pvalue,
    padj_adjusted_model = res2$padj,
    stringsAsFactors = FALSE)
  tab$symbol <- ids_to_symbols(tab$gene_id)

  out <- tab[!is.na(tab$symbol) & tab$symbol %in% genes, ]
  # also report how much of the group effect the marker score itself absorbed
  attr(out, "marker_genes") <- attr(score, "genes_used")
  attr(out, "marker_score") <- score
  out
}

# ---- METHODS S2.3(b) gene-change vs marker-change correlation -------------

# Across samples, does the gene's expression track the cell-identity marker
# score? A strong positive correlation means the signal is largely a
# composition readout. Reported with r, 95% CI, and p.
correlate_with_composition <- function(mat_vst, meta, panel_symbols,
                                       genes = GENES_OF_INTEREST,
                                       method = "pearson", de_res = NULL) {
  score <- marker_score(mat_vst, panel_symbols)
  if (is.null(score)) return(NULL)
  panel <- resolve_panel(genes, rownames(mat_vst))
  mpanel <- resolve_panel(panel_symbols, rownames(mat_vst))
  mpanel <- mpanel[mpanel$found, ]

  # The composite score averages z-scored markers. If the markers do not all
  # move in the same direction between groups, that average partially cancels
  # and understates the true association -- so we report the composite AND
  # every marker individually, and warn when the panel is heterogeneous.
  if (!is.null(de_res)) {
    md <- de_res[de_res$symbol %in% mpanel$symbol, ]
    if (nrow(md) > 1) {
      dirs <- sign(md$log2FC_shrunk)
      if (length(unique(dirs[!is.na(dirs)])) > 1) {
        log_warn("marker panel is directionally heterogeneous (%s) -- composite score partially cancels; rely on the per-marker rows",
                 paste(sprintf("%s%s", md$symbol,
                               ifelse(dirs > 0, "+", "-")), collapse = " "))
      }
    }
  }

  cor_one <- function(expr, ref, gene, ref_label) {
    if (stats::sd(expr) == 0 || stats::sd(ref) == 0) return(NULL)
    ct <- stats::cor.test(expr, ref, method = method)
    data.frame(
      gene = gene, reference = ref_label, n = length(expr),
      r = unname(ct$estimate),
      ci_low = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
      ci_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
      statistic = unname(ct$statistic),
      p_value = ct$p.value, method = method,
      stringsAsFactors = FALSE)
  }

  sc <- as.numeric(score[colnames(mat_vst)])
  rows <- list()
  for (i in seq_len(nrow(panel))) {
    if (!panel$found[i]) next
    expr <- as.numeric(mat_vst[panel$id[i], ])
    rows[[length(rows) + 1]] <- cor_one(expr, sc, panel$symbol[i],
                                        "composite marker score")
    for (j in seq_len(nrow(mpanel))) {
      rows[[length(rows) + 1]] <- cor_one(
        expr, as.numeric(mat_vst[mpanel$id[j], ]),
        panel$symbol[i], mpanel$symbol[j])
    }
  }
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) return(NULL)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  attr(out, "marker_genes") <- attr(score, "genes_used")
  out
}

# ---- interpretation helper ------------------------------------------------

# Turns the quantitative composition evidence into the explicit verdict the
# analysis standard asks for, rather than leaving it to prose.
composition_verdict <- function(goi_de, marker_de, adj_res = NULL,
                                alpha = ALPHA) {
  g <- goi_de[goi_de$symbol %in% GENES_OF_INTEREST, ]
  m <- marker_de
  if (nrow(g) == 0 || is.null(m) || nrow(m) == 0) return("insufficient data")

  out <- character(0)
  for (i in seq_len(nrow(g))) {
    gene <- g$symbol[i]
    lfc  <- g$log2FC_shrunk[i]
    sig  <- !is.na(g$padj[i]) && g$padj[i] < alpha

    same_dir <- sign(m$log2FC_shrunk) == sign(lfc)
    m_sig <- !is.na(m$padj) & m$padj < alpha
    n_conc <- sum(same_dir & m_sig, na.rm = TRUE)

    survives <- NA
    if (!is.null(adj_res) && gene %in% adj_res$symbol) {
      a <- adj_res[adj_res$symbol == gene, ][1, ]
      survives <- !is.na(a$padj_adjusted_model) && a$padj_adjusted_model < alpha
    }

    verdict <- if (!sig) {
      "no significant change"
    } else if (n_conc >= ceiling(nrow(m) / 2) && identical(survives, FALSE)) {
      "LIKELY COMPOSITION SHIFT (markers move together; effect lost after adjustment)"
    } else if (n_conc >= ceiling(nrow(m) / 2)) {
      "AMBIGUOUS (markers move with gene; check adjusted model and deconvolution)"
    } else if (identical(survives, TRUE)) {
      "LIKELY GENUINE REGULATION (markers flat; effect survives adjustment)"
    } else {
      "LIKELY GENUINE REGULATION (markers flat)"
    }
    out <- c(out, sprintf("%s: log2FC=%.3f padj=%.3g -> %s",
                          gene, lfc, g$padj[i], verdict))
  }
  paste(out, collapse = " | ")
}
