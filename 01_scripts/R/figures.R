# ---------------------------------------------------------------------------
# R/figures.R -- the seven core figure types (docs/METHODS.md S6).
#
#   1 fig_dotplot        gene-of-interest effect across conditions
#   2 fig_volcano        per-dataset DE, genes of interest highlighted
#   3 fig_forest         meta-analysis pooled effects
#   4 fig_boxplot        genes of interest, case vs control per dataset
#   5 fig_heatmap        genes of interest + markers across conditions
#   6 fig_pca            sample structure / batch, before vs after
#   7 fig_correlation    gene vs comparator across atlas studies
#
# All obey METHODS S6: axis labels with units, n= annotated, one shared
# theme from 00_theme.R. Every function returns a ggplot; save with save_fig().
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

.units_vst  <- "VST-normalized expression (log2 scale)"
.units_lfc  <- expression(log[2]~"fold change (shrunken)")

# ---- 1. dot plot ----------------------------------------------------------

# Effect size per disease/dataset; point size = significance, colour = direction.
fig_dotplot <- function(effect_table, gene_col = "symbol",
                        title = "Genes of interest across conditions") {
  d <- effect_table[!is.na(effect_table$log2FC_shrunk), ]
  d$direction <- ifelse(is.na(d$padj) | d$padj >= ALPHA, "NS",
                        ifelse(d$log2FC_shrunk > 0, "Up", "Down"))
  d$neglog10p <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$label <- sprintf("%s\n(n=%d)", d$study, d$n_samples)

  ggplot(d, aes(x = log2FC_shrunk, y = reorder(label, log2FC_shrunk))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = log2FC_shrunk - 1.96 * lfcSE_shrunk,
                       xmax = log2FC_shrunk + 1.96 * lfcSE_shrunk),
                   height = 0.2, colour = "grey55", linewidth = 0.4) +
    geom_point(aes(size = neglog10p, colour = direction)) +
    facet_wrap(as.formula(paste("~", gene_col)), scales = "free_x") +
    scale_colour_manual(values = PAL_DIRECTION, name = "Direction\n(BH-FDR)") +
    scale_size_continuous(name = expression(-log[10]~"BH-adjusted p"),
                          range = c(1.5, 6)) +
    labs(title = title,
         subtitle = "Points show shrunken effect sizes; bars are 95% CI",
         x = .units_lfc, y = NULL,
         caption = caption_stats()) +
    theme(axis.text.y = element_text(size = rel(0.8)))
}

# ---- 2. volcano -----------------------------------------------------------

fig_volcano <- function(de_res, gse, highlight = GENES_OF_INTEREST,
                        markers = NULL, alpha = ALPHA) {
  d <- de_res[!is.na(de_res$padj) & !is.na(de_res$log2FC_shrunk), ]
  d$direction <- ifelse(d$padj >= alpha, "NS",
                        ifelse(d$log2FC_shrunk > 0, "Up", "Down"))
  d$neglog10p <- -log10(pmax(d$padj, .Machine$double.xmin))

  lab <- d[!is.na(d$symbol) & d$symbol %in% c(highlight, markers), ]
  n_sig <- sum(d$padj < alpha, na.rm = TRUE)

  p <- ggplot(d, aes(x = log2FC_shrunk, y = neglog10p)) +
    geom_point(aes(colour = direction), size = 0.7, alpha = 0.5) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed",
               colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = PAL_DIRECTION, name = NULL) +
    labs(title = sprintf("%s: differential expression", gse),
         subtitle = sprintf("%s genes tested, %s significant at BH-FDR < %.2f",
                            fmt_n(nrow(d)), fmt_n(n_sig), alpha),
         x = .units_lfc,
         y = expression(-log[10]~"BH-adjusted p"),
         caption = caption_stats("Labelled: genes of interest and cell-identity markers."))

  if (nrow(lab) > 0) {
    p <- p +
      geom_point(data = lab, colour = "black", size = 2.2, shape = 21,
                 fill = "white", stroke = 0.7)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, aes(label = symbol), size = 3, min.segment.length = 0,
        box.padding = 0.5, max.overlaps = 30, seed = SEED)
    }
  }
  p
}

# ---- 3. forest ------------------------------------------------------------

# Built in ggplot rather than metafor::forest() so the styling matches the rest
# of the figure set (METHODS S6).
fig_forest <- function(pooled, gene, effect_label = "log2 fold change") {
  d <- pooled$data
  s <- pooled$summary
  fit <- pooled$fit

  is_lfc <- "log2FC_shrunk" %in% colnames(d)
  # correlation tables carry `n`; disease effect tables carry `n_samples`
  n_vec <- if (!is.null(d$n)) d$n else if (!is.null(d$n_samples)) d$n_samples else
    rep(NA_integer_, nrow(d))

  study_df <- data.frame(
    study = d$study,
    est   = if (is_lfc) d$log2FC_shrunk else d$r,
    lo    = if (is_lfc) d$log2FC_shrunk - 1.96 * d$lfcSE_shrunk else d$ci_low,
    hi    = if (is_lfc) d$log2FC_shrunk + 1.96 * d$lfcSE_shrunk else d$ci_high,
    n     = n_vec,
    type  = "Study", stringsAsFactors = FALSE)

  pooled_df <- data.frame(
    study = "Pooled (RE)",
    est = if (is_lfc) s$pooled_log2FC else s$pooled_r,
    lo  = if (is_lfc) s$ci_low else s$ci_low_r,
    hi  = if (is_lfc) s$ci_high else s$ci_high_r,
    n = if (all(is.na(n_vec))) NA_integer_ else sum(n_vec, na.rm = TRUE),
    type = "Pooled", stringsAsFactors = FALSE)

  all_df <- rbind(study_df, pooled_df)
  # METHODS S6 requires n= annotated; fall back gracefully if unavailable
  all_df$label <- ifelse(is.na(all_df$n), all_df$study,
                         sprintf("%s (n=%d)", all_df$study, all_df$n))
  all_df$label <- factor(all_df$label, levels = rev(all_df$label))

  ggplot(all_df, aes(x = est, y = label, colour = type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5) +
    geom_point(aes(shape = type, size = type)) +
    scale_colour_manual(values = c(Study = "grey25", Pooled = "#C0392B"),
                        guide = "none") +
    scale_shape_manual(values = c(Study = 16, Pooled = 18), guide = "none") +
    scale_size_manual(values = c(Study = 2.4, Pooled = 4.5), guide = "none") +
    labs(title = sprintf("%s: random-effects meta-analysis", gene),
         subtitle = heterogeneity_text(pooled),
         x = effect_label, y = NULL,
         caption = "Bars are 95% confidence intervals. Pooled estimate from metafor::rma (REML).")
}

# ---- 4. boxplot -----------------------------------------------------------

fig_boxplot <- function(mat_vst, meta, gse, genes = GENES_OF_INTEREST,
                        group_col = "group", de_res = NULL) {
  panel <- resolve_panel(genes, rownames(mat_vst))
  panel <- panel[panel$found, ]
  if (nrow(panel) == 0) return(NULL)

  long <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i) {
    data.frame(gene = panel$symbol[i],
               expr = as.numeric(mat_vst[panel$id[i], ]),
               group = meta[[group_col]][match(colnames(mat_vst), meta$sample_id)],
               stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$group), ]

  labs_n <- label_with_n(long$group[long$gene == long$gene[1]])

  # annotate each facet with the BH-adjusted p from the DE model
  ann <- NULL
  if (!is.null(de_res)) {
    sub <- de_res[de_res$symbol %in% panel$symbol, ]
    if (nrow(sub) > 0) {
      ann <- data.frame(
        gene = sub$symbol,
        label = sprintf("BH-FDR p = %.3g", sub$padj),
        y = tapply(long$expr, long$gene, max)[sub$symbol],
        stringsAsFactors = FALSE)
    }
  }

  p <- ggplot(long, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.75,
                 linewidth = 0.4) +
    geom_jitter(width = 0.16, size = 1.1, alpha = 0.6, colour = "grey20") +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = PAL_GROUP, guide = "none") +
    scale_x_discrete(labels = labs_n) +
    labs(title = sprintf("%s: genes of interest by group", gse),
         x = NULL, y = .units_vst,
         caption = caption_stats())

  if (!is.null(ann)) {
    p <- p + geom_text(data = ann, aes(x = 1.5, y = y * 1.02, label = label),
                       inherit.aes = FALSE, size = 3, colour = "grey20")
  }
  p
}

# ---- 5. heatmap -----------------------------------------------------------

# Genes of interest + cell-identity markers, z-scored per gene so composition
# and gene-of-interest patterns are directly comparable by eye -- with the statistics
# reported separately, never inferred from the colours alone.
fig_heatmap <- function(mat_vst, meta, gse, genes = NULL, group_col = "group") {
  if (is.null(genes)) {
    genes <- c(GENES_OF_INTEREST, unlist(MARKERS, use.names = FALSE))
  }
  panel <- resolve_panel(genes, rownames(mat_vst))
  panel <- panel[panel$found, ]
  if (nrow(panel) < 2) return(NULL)

  m <- mat_vst[panel$id, , drop = FALSE]
  rownames(m) <- panel$symbol
  z <- t(scale(t(m)))

  ord <- order(meta[[group_col]][match(colnames(z), meta$sample_id)])
  z <- z[, ord, drop = FALSE]
  grp <- meta[[group_col]][match(colnames(z), meta$sample_id)]

  long <- data.frame(
    gene = rep(rownames(z), times = ncol(z)),
    sample = rep(colnames(z), each = nrow(z)),
    z = as.vector(z),
    group = rep(grp, each = nrow(z)),
    stringsAsFactors = FALSE)
  long$gene <- factor(long$gene, levels = rev(panel$symbol))
  long$sample <- factor(long$sample, levels = colnames(z))

  n_txt <- paste(sprintf("%s n=%d", names(table(grp)), table(grp)),
                 collapse = "; ")

  ggplot(long, aes(x = sample, y = gene, fill = z)) +
    geom_tile() +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = PAL_DIVERGING[1], mid = PAL_DIVERGING[2],
                         high = PAL_DIVERGING[3], midpoint = 0,
                         name = "z-score\n(per gene)") +
    labs(title = sprintf("%s: genes of interest and cell-identity markers", gse),
         subtitle = sprintf("VST expression, z-scored per gene. %s", n_txt),
         x = "Sample", y = NULL,
         caption = "Cell-identity markers included so changes in the genes of interest can be read against composition shifts.") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank())
}

# ---- 6. PCA ---------------------------------------------------------------

# METHODS S3: PCA before and after correction, with % variance on axes
# and the variance-explained table reported alongside.
fig_pca <- function(mat_vst, meta, gse, colour_col = "group",
                    shape_col = NULL, subtitle = NULL) {
  pca <- stats::prcomp(t(mat_vst), center = TRUE, scale. = FALSE)
  vp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  idx <- match(rownames(pca$x), meta$sample_id)
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  stringsAsFactors = FALSE)
  d$colour_var <- meta[[colour_col]][idx]
  if (!is.null(shape_col) && !is.null(meta[[shape_col]])) {
    d$shape_var <- as.factor(meta[[shape_col]][idx])
  }

  n_txt <- paste(sprintf("%s n=%d", names(table(d$colour_var)),
                         table(d$colour_var)), collapse = "; ")

  p <- ggplot(d, aes(x = PC1, y = PC2, colour = colour_var)) +
    { if (!is.null(d$shape_var)) geom_point(aes(shape = shape_var), size = 2.6, alpha = 0.9)
      else geom_point(size = 2.6, alpha = 0.9) } +
    stat_ellipse(type = "norm", level = 0.68, linewidth = 0.4,
                 linetype = "dashed", show.legend = FALSE) +
    labs(title = sprintf("%s: sample structure (PCA)", gse),
         subtitle = subtitle %||% sprintf("VST expression, all genes. %s", n_txt),
         x = sprintf("PC1 (%.1f%% of variance)", vp[1]),
         y = sprintf("PC2 (%.1f%% of variance)", vp[2]),
         colour = colour_col, shape = shape_col)

  if (all(names(PAL_GROUP) %in% unique(as.character(d$colour_var)))) {
    p <- p + scale_colour_manual(values = PAL_GROUP)
  }
  p
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Before/after batch correction, side by side.
fig_pca_before_after <- function(mat_before, mat_after, meta, gse,
                                 colour_col = "group", batch_col = NULL) {
  require_pkg("patchwork", "side-by-side PCA (METHODS S3)")
  p1 <- fig_pca(mat_before, meta, gse, colour_col, batch_col,
                subtitle = "Before batch correction")
  p2 <- fig_pca(mat_after, meta, gse, colour_col, batch_col,
                subtitle = "After ComBat (cell type protected)")
  patchwork::wrap_plots(p1, p2, ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("%s: batch correction effect", gse),
      caption = "Quantitative variance decomposition reported in 05_results/.")
}

# Companion bar chart: % of total variance explained by batch vs group.
fig_variance_explained <- function(var_summary, gse,
                                   subtitle = "Sum over the top 5 principal components") {
  d <- var_summary
  d$covariate <- factor(d$covariate,
                        levels = d$covariate[order(d$total_var_explained_pct)])
  ggplot(d, aes(x = total_var_explained_pct, y = covariate)) +
    geom_col(fill = "#4C72B0", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", total_var_explained_pct)),
              hjust = -0.15, size = 3.2, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = sprintf("%s: variance attributable to each factor", gse),
         subtitle = subtitle,
         x = "Total dataset variance explained (%)", y = NULL,
         caption = "PC ~ covariate regression; per-PC variance weighted by that PC's share of total variance.")
}

# ---- 7. correlation -------------------------------------------------------

fig_correlation <- function(mat_vst, meta, study_id,
                            gene_a = GENES_OF_INTEREST[1],
                            gene_b = GENES_OF_INTEREST[2], colour_col = NULL) {
  panel <- resolve_panel(c(gene_a, gene_b), rownames(mat_vst))
  if (!all(panel$found)) return(NULL)

  d <- data.frame(x = as.numeric(mat_vst[panel$id[1], ]),
                  y = as.numeric(mat_vst[panel$id[2], ]),
                  stringsAsFactors = FALSE)
  if (!is.null(colour_col) && !is.null(meta[[colour_col]])) {
    d$grp <- meta[[colour_col]][match(colnames(mat_vst), meta$sample_id)]
  }
  ct <- stats::cor.test(d$x, d$y)

  p <- ggplot(d, aes(x = x, y = y)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#C0392B",
                fill = "#C0392B", alpha = 0.12, linewidth = 0.6) +
    { if (!is.null(d$grp)) geom_point(aes(colour = grp), size = 2.2, alpha = 0.85)
      else geom_point(size = 2.2, alpha = 0.85, colour = "grey20") } +
    labs(title = sprintf("%s: %s vs %s", study_id, gene_a, gene_b),
         subtitle = sprintf("Pearson r = %.3f [%.3f, %.3f], p = %.3g, n = %d",
                            ct$estimate, ct$conf.int[1], ct$conf.int[2],
                            ct$p.value, nrow(d)),
         x = sprintf("%s -- %s", gene_a, .units_vst),
         y = sprintf("%s -- %s", gene_b, .units_vst),
         colour = colour_col,
         caption = "Within-study correlation. Pooling across studies is done by random-effects meta-analysis, not by combining samples.")
  p
}
