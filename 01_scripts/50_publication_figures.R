#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 50_publication_figures.R -- the final figure set -> 07_final_figures/
#
# Writes ONLY to 07_final_figures/ (600 dpi PNG + vector PDF). 06_figures/ is
# the working archive written by the earlier stages and is never touched here,
# so an exploratory plot can never be mistaken for a final figure.
#
#   Figure 1  dot plot      effect size per gene per study, grouped by condition
#   Figure 2  volcano       per-dataset DE with the genes of interest ringed
#   Figure 3  forest        random-effects meta-analysis, overall and by
#                           tissue class, with the moderator test
#   Figure 4  boxplots      expression by group, per dataset
#   Figure 5  heatmap       genes of interest read against cell-identity markers
#   Figure 6  PCA           sample structure before vs after batch correction
#   Figure 7  correlation   gene vs comparator across the atlas studies
#
# Every caption states the test, the multiple-testing method and the shrinkage
# estimator, and every panel carries its n -- a figure that leaves its own
# statistics to the surrounding prose is not finished (docs/METHODS.md S6).
#
# Run the pipeline stages first: 10_run_dataset.R for each dataset, then
# 20_meta_analysis.R, then 30_coexpression_atlas.R.
#
# Usage:
#   Rscript 01_scripts/50_publication_figures.R          # all figures
#   Rscript 01_scripts/50_publication_figures.R 3        # just Figure 3
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(data.table)
})

source("01_scripts/00_config.R")
source("01_scripts/00_theme_publication.R")
for (f in c("utils", "annotate", "meta")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

args <- commandArgs(trailingOnly = TRUE)
WHICH <- if (length(args)) as.integer(args[1]) else 0   # 0 = all

# ---- what these figures are about -----------------------------------------
# Everything below is driven by the config, so the same script produces the
# figure set for whatever genes and conditions you registered.

GENES     <- GENES_OF_INTEREST
G1        <- GENES[1]                                   # primary gene
G2        <- if (length(GENES) > 1) GENES[2] else NA    # comparator / paralog
GENE_TXT  <- paste(GENES, collapse = " and ")

# Human-readable condition labels; override CONDITION_LABELS in 00_config.R to
# get prettier names than the registry codes. Order here is the order the
# facets appear in, so group related tissues together.
DISEASE_LABEL <- CONDITION_LABELS
DISEASE_ORDER <- names(DISEASE_LABEL)

# Figures 2, 4 and 5 show a readable SUBSET of datasets -- a volcano grid of
# twenty panels communicates nothing. Take the PRIMARY dataset of each
# condition, preferring ones that have actually been run, and keep tissue and
# blood both represented so the compartment comparison in Figure 5 is possible.
# Override by setting FEATURED_GSES before sourcing, or edit the rule.
if (!exists("FEATURED_GSES")) {
  FEATURED_GSES <- local({
    d <- DATASETS[DATASETS$arm == "case_control" & DATASETS$de_capable, ]
    has_de <- file.exists(file.path(PATHS$de,
                                    paste0(d$gse, "_effect_for_meta.csv")))
    d <- d[has_de, ]
    if (!nrow(d)) return(character(0))
    d <- d[order(match(d$use, c("PRIMARY", "VALIDATION", "SECONDARY")),
                 match(d$condition, DISEASE_ORDER)), ]
    pick <- unlist(lapply(split(d$gse, d$condition), function(x) x[1]))
    pick <- pick[order(match(DATASETS$condition[match(pick, DATASETS$gse)],
                             DISEASE_ORDER))]
    unname(pick)
  })
}
if (!length(FEATURED_GSES)) {
  stop("No per-dataset results found in ", PATHS$de,
       ".\nRun 01_scripts/10_run_dataset.R on your datasets first.")
}
message(sprintf("[figures] featured datasets: %s",
                paste(FEATURED_GSES, collapse = ", ")))

load_effects <- function() {
  fs <- list.files(PATHS$de, pattern = "_effect_for_meta\\.csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
  d$tissue_class <- DATASETS$tissue_class[match(d$study, DATASETS$gse)]
  d$use <- DATASETS$use[match(d$study, DATASETS$gse)]
  d$disease <- factor(d$disease, levels = DISEASE_ORDER)
  d
}

# ===========================================================================
# FIGURE 1 -- cross-condition dot plot
# ===========================================================================

fig1_dotplot <- function() {
  d <- load_effects()
  d <- d[!is.na(d$log2FC_shrunk) & !is.na(d$lfcSE_shrunk), ]

  d$direction <- ifelse(is.na(d$padj) | d$padj >= ALPHA, "NS",
                        ifelse(d$log2FC_shrunk > 0, "Up", "Down"))
  d$direction <- factor(d$direction, levels = c("Down", "NS", "Up"))
  d$lo <- d$log2FC_shrunk - 1.96 * d$lfcSE_shrunk
  d$hi <- d$log2FC_shrunk + 1.96 * d$lfcSE_shrunk

  # each row is one study; label carries the n per group, per the brief
  d$row_label <- sprintf("%s  (%d vs %d)", d$study, d$n_control, d$n_case)
  # The primary gene leads; alphabetical order would otherwise decide it
  d$symbol <- factor(d$symbol, levels = GENES)
  d$disease_lab <- factor(DISEASE_LABEL[as.character(d$disease)],
                          levels = DISEASE_LABEL[DISEASE_ORDER])

  # order rows within each disease by effect size, keeping diseases grouped
  d <- d[order(d$disease, d$log2FC_shrunk), ]
  lv <- unique(d$row_label[order(d$disease, -d$log2FC_shrunk)])
  d$row_label <- factor(d$row_label, levels = rev(lv))

  # significance stars give a second, non-colour encoding of the same fact
  d$star <- ifelse(is.na(d$padj), "",
            ifelse(d$padj < 0.001, "***",
            ifelse(d$padj < 0.01, "**",
            ifelse(d$padj < ALPHA, "*", ""))))

  n_studies <- length(unique(d$study))
  n_samples <- sum(tapply(d$n_samples, d$study, function(x) x[1]))

  ggplot(d, aes(x = log2FC_shrunk, y = row_label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = INK$muted,
               linetype = "22") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                   linewidth = 0.55, colour = INK$faint) +
    geom_point(aes(fill = direction), shape = 21, size = 2.9,
               stroke = 0.45, colour = INK$surface) +
    geom_text(aes(x = hi, label = star), hjust = -0.35, vjust = 0.78,
              size = 2.5, colour = INK$secondary, family = PUB_FONT) +
    facet_grid(disease_lab ~ symbol, scales = "free_y", space = "free_y",
               switch = "y") +
    scale_fill_manual(values = PUB_DIRECTION, name = NULL,
                      labels = c(Down = "Lower in disease",
                                 NS = "Not significant",
                                 Up = "Higher in disease"),
                      drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.14))) +
    labs(
      title = sprintf("%s across %d %s", GENE_TXT, length(DISEASE_ORDER),
                      ifelse(length(DISEASE_ORDER) == 1, "condition", "conditions")),
      subtitle = sprintf(
        "Shrunken effect sizes with 95%% confidence intervals · %d datasets, %s samples\nRow labels show n control vs n case",
        n_studies, fmt_n(n_samples)),
      x = LAB_LFC, y = NULL,
      caption = paste0(
        "Effect sizes from DESeq2 (apeglm shrinkage) or limma-trend with ashr shrinkage where the series ships normalised values.\n",
        "Significance by Benjamini–Hochberg FDR: * < 0.05, ** < 0.01, *** < 0.001. Bars are 95% CI; absent points were below the expression filter.")
    ) +
    theme_pub(grid = "x") +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold",
                                       colour = INK$primary,
                                       margin = margin(r = 7)),
      strip.text.x = element_text(face = "bold", size = rel(1.05)),
      axis.text.y = element_text(size = rel(0.86), colour = INK$secondary,
                                 hjust = 0),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.spacing.y = unit(3.5, "pt"),
      panel.spacing.x = unit(13, "pt"),
      legend.position = "top",
      legend.justification = c(0, 0.5),
      legend.location = "plot",
      legend.direction = "horizontal",
      legend.margin = margin(b = 4)
    )
}

# ---------------------------------------------------------------------------

if (WHICH %in% c(0, 1)) {
  log_step("Figure 1 -- cross-condition dot plot")
  save_pub(fig1_dotplot(), "Figure1_dotplot_genes_across_conditions",
           width = 8.4, height = 7.2)
}

if (WHICH == 0) log_session_info("50_publication_figures")

# ===========================================================================
# FIGURE 2 -- volcano panels
# ===========================================================================

fig2_volcano <- function() {
  sel <- FEATURED_GSES
  rows <- list()
  for (g in sel) {
    f <- file.path(PATHS$de, paste0(g, "_DE_full.csv"))
    if (!file.exists(f)) next
    d <- data.table::fread(f, data.table = FALSE)
    d <- d[!is.na(d$padj) & !is.na(d$log2FC_shrunk), ]
    e <- utils::read.csv(file.path(PATHS$de, paste0(g, "_effect_for_meta.csv")))
    d$study <- g
    d$n_control <- e$n_control[1]; d$n_case <- e$n_case[1]
    d$disease <- DISEASE_LABEL[as.character(e$disease[1])]
    rows[[g]] <- d[, c("study","disease","n_control","n_case","symbol",
                       "log2FC_shrunk","padj")]
  }
  d <- do.call(rbind, rows)
  d$direction <- factor(ifelse(d$padj >= ALPHA, "NS",
                        ifelse(d$log2FC_shrunk > 0, "Up", "Down")),
                        levels = c("Down","NS","Up"))
  d$negLog <- -log10(pmax(d$padj, .Machine$double.xmin))
  # cap the y axis so a handful of extreme genes do not flatten every panel
  ycap <- stats::quantile(d$negLog, 0.999, na.rm = TRUE)
  d$negLog_c <- pmin(d$negLog, ycap)
  d$capped <- d$negLog > ycap

  lab <- d[d$symbol %in% GENES_OF_INTEREST, ]
  hdr <- unique(d[, c("study","disease","n_control","n_case")])
  hdr$panel <- sprintf("%s · %s\n%d control vs %d case",
                       hdr$disease, hdr$study, hdr$n_control, hdr$n_case)
  d$panel  <- hdr$panel[match(d$study, hdr$study)]
  lab$panel <- hdr$panel[match(lab$study, hdr$study)]
  ord <- hdr$panel[match(sel, hdr$study)]
  d$panel <- factor(d$panel, levels = ord); lab$panel <- factor(lab$panel, levels = ord)

  p <- ggplot(d, aes(log2FC_shrunk, negLog_c)) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "22",
               colour = INK$muted, linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "22", colour = INK$faint,
               linewidth = 0.3) +
    geom_point(aes(colour = direction), size = 0.30, alpha = 0.42, shape = 16) +
    # Highlight marker is deliberately NEUTRAL, not direction-coloured. The
    # background cloud is already red/blue/grey by direction, so a
    # direction-coloured highlight would camouflage the very genes it is meant
    # to make findable -- and the point's position on the x-axis already states
    # the direction unambiguously. A white halo lifts the marker off the cloud;
    # the near-black ring carries "gene of interest" without competing with the
    # up/down palette used in Figures 1, 3, 4 and 5.
    geom_point(data = lab, shape = 21, size = 4.6, stroke = 0,
               fill = INK$surface, colour = NA, alpha = 0.95) +
    geom_point(data = lab, shape = 21, size = 2.9, stroke = 1.05,
               fill = INK$surface, colour = INK$primary) +
    facet_wrap(~ panel, ncol = 3, scales = "free") +
    scale_colour_manual(values = PUB_DIRECTION, name = NULL,
                        labels = c(Down="Lower in disease", NS="Not significant",
                                   Up="Higher in disease"),
                        guide = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
    labs(title = sprintf("Differential expression in %d representative datasets",
                         length(unique(d$study))),
         subtitle = sprintf("%s highlighted against the transcriptome-wide background",
                            GENE_TXT),
         x = LAB_LFC, y = LAB_PADJ,
         caption = paste0(
           "Benjamini\u2013Hochberg FDR; dashed horizontal line at FDR = 0.05. y-axis capped at the 99.9th percentile so single extreme genes do not flatten the panels.\n",
           sprintf("%s are ringed in black; their position on the x-axis, not the ring colour, carries the direction of change. A gene absent from a panel fell below\n", GENE_TXT),
           "that dataset's expression filter. Where a gene sits close to the detection floor in a tissue, its effect size is less well determined than the interval\n",
           "alone suggests -- 40_sensitivity_tests.R quantifies that directly.")) +
    theme_pub(grid = "both") +
    theme(strip.text = element_text(size = rel(0.86), lineheight = 1.15),
          legend.position = "top", legend.justification = c(0, 0.5),
          legend.location = "plot", legend.direction = "horizontal",
          panel.grid.major = element_line(colour = INK$rule, linewidth = 0.2))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = lab, aes(label = symbol), size = 2.6, family = PUB_FONT,
      fontface = "bold", colour = INK$primary, seed = SEED,
      min.segment.length = 0, box.padding = 1.15, point.padding = 0.55,
      force = 12, force_pull = 0.4, nudge_y = 1.2,
      segment.colour = INK$secondary, segment.size = 0.35,
      segment.alpha = 0.9, max.overlaps = 60,
      bg.color = INK$surface, bg.r = 0.16)
  }
  p
}

# ===========================================================================
# FIGURE 3 -- forest plots (overall + tissue-stratified)
# ===========================================================================

fig3_forest <- function() {
  eff <- load_effects()
  sg  <- utils::read.csv(file.path(PATHS$results, "meta_subgroup_by_tissue.csv"))
  pl  <- utils::read.csv(file.path(PATHS$results, "meta_pooled_effect_sizes.csv"))

  # ONE row order shared by both gene panels. Ordering each panel by its own
  # effect size put the same GSE on different rows, making the two genes for
  # one dataset impossible to compare across panels. Order is fixed by tissue
  # class, then by the primary gene's effect, and reused for the comparator.
  ref <- eff[eff$symbol == G1 & !is.na(eff$lfcSE_shrunk) & eff$lfcSE_shrunk > 0, ]
  ref <- ref[order(ref$tissue_class, ref$log2FC_shrunk), ]
  study_lab <- setNames(sprintf("%s  (%d vs %d)", ref$study, ref$n_control, ref$n_case),
                        ref$study)
  study_tc  <- setNames(ref$tissue_class, ref$study)

  rows <- list()
  for (gene in GENES) {
    d <- eff[eff$symbol == gene & !is.na(eff$lfcSE_shrunk) & eff$lfcSE_shrunk > 0, ]
    d <- d[match(names(study_lab), d$study), ]      # SAME order in both panels
    d$study <- names(study_lab)
    rows[[length(rows)+1]] <- data.frame(
      gene = gene, label = unname(study_lab),
      block = ifelse(unname(study_tc) == "blood", "Blood", "Solid tissue"),
      est = d$log2FC_shrunk,
      lo = d$log2FC_shrunk - 1.96*d$lfcSE_shrunk,
      hi = d$log2FC_shrunk + 1.96*d$lfcSE_shrunk,
      kind = ifelse(unname(study_tc) == "blood", "Blood study", "Tissue study"),
      klab = NA_character_, stringsAsFactors = FALSE)
    for (tc in c("tissue","blood")) {
      r <- sg[sg$gene == gene & sg$tissue_class == tc, ]
      if (!nrow(r)) next
      rows[[length(rows)+1]] <- data.frame(
        gene = gene,
        label = sprintf("POOLED \u2014 %s", ifelse(tc=="blood","blood","solid tissue")),
        block = "Pooled", est = r$pooled_log2FC, lo = r$ci_low, hi = r$ci_high,
        kind = ifelse(tc=="blood","Pooled blood","Pooled tissue"),
        klab = sprintf("k = %d", r$k), stringsAsFactors = FALSE)
    }
    r <- pl[pl$gene == gene, ]
    if (nrow(r)) rows[[length(rows)+1]] <- data.frame(
      gene = gene, label = "POOLED \u2014 all studies",
      block = "Pooled", est = r$pooled_log2FC, lo = r$ci_low, hi = r$ci_high,
      kind = "Pooled overall", klab = sprintf("k = %d", r$k),
      stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  d <- d[!is.na(d$est), ]
  d$gene <- factor(d$gene, levels = GENES)
  d$block <- factor(d$block, levels = c("Solid tissue","Blood","Pooled"))
  d$kind <- factor(d$kind, levels = c("Tissue study","Blood study",
                                      "Pooled tissue","Pooled blood","Pooled overall"))
  ord <- unique(d$label[order(d$block)])
  d$label <- factor(d$label, levels = rev(ord))
  d$is_pool <- d$block == "Pooled"

  het <- sapply(GENES, function(g) {
    r <- pl[pl$gene == g, ]
    sprintf("%s I\u00b2 = %.0f%%", g, r$I2) })
  mod <- sapply(GENES, function(g) {
    r <- sg[sg$gene == g & sg$tissue_class == "MODERATOR TEST", ]
    sprintf("%s p = %s", g, fmt_p(r$p_value)) })

  # shaded bands behind each block, so tissue / blood / pooled read as sections
  bands <- do.call(rbind, lapply(levels(d$block), function(b) {
    idx <- which(rev(ord) %in% d$label[d$block == b])
    data.frame(block = b, ymin = min(idx) - 0.5, ymax = max(idx) + 0.5,
               stringsAsFactors = FALSE) }))
  bands$fill <- c("Solid tissue" = "#eef6f2", "Blood" = "#f0eef7",
                  "Pooled" = "#fbeceb")[bands$block]

  ggplot(d, aes(est, label)) +
    geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
                                fill = I(fill)), inherit.aes = FALSE) +
    # hard rule marking the break before the pooled summaries
    geom_hline(yintercept = min(which(rev(ord) %in% d$label[d$block == "Pooled"])) - 0.5,
               colour = INK$secondary, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "22", colour = INK$muted, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = kind, linewidth = is_pool),
                   height = 0) +
    geom_point(aes(colour = kind, shape = kind, size = kind)) +
    geom_text(data = d[!is.na(d$klab), ], aes(x = hi, label = klab),
              hjust = -0.28, size = 2.3, family = PUB_FONT,
              colour = INK$secondary) +
    facet_wrap(~ gene, ncol = 2) +
    scale_colour_manual(
      values = c("Tissue study"   = PUB_TISSUE[["tissue"]],
                 "Blood study"    = PUB_TISSUE[["blood"]],
                 "Pooled tissue"  = "#0f7a52",
                 "Pooled blood"   = "#2f2470",
                 "Pooled overall" = "#b32726"), name = NULL) +
    scale_shape_manual(values = c("Tissue study"=16,"Blood study"=16,
                                  "Pooled tissue"=23,"Pooled blood"=23,
                                  "Pooled overall"=23), name = NULL) +
    scale_size_manual(values = c("Tissue study"=2.1,"Blood study"=2.1,
                                 "Pooled tissue"=5.2,"Pooled blood"=5.2,
                                 "Pooled overall"=6.0), name = NULL) +
    scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.0), guide = "none") +
    labs(title = sprintf("Random-effects meta-analysis of %s", GENE_TXT),
         subtitle = sprintf("Pooled with metafor (REML, Fisher-z for correlations). Heterogeneity: %s. Tissue-class moderator test: %s.\nStudies grouped by compartment; diamonds are pooled estimates. Row order is identical in both panels.",
                            paste(het, collapse = ", "), paste(mod, collapse = ", ")),
         x = LAB_LFC, y = NULL,
         caption = paste0(
           "Bars are 95% confidence intervals; study labels show n control vs n case. Diamonds below the rule are random-effects pooled estimates, not observations.\n",
           "A gene pooled over fewer studies than another fell below the expression filter in the remaining ones: that is a measurability difference, not a missing\n",
           "result. Read the pooled estimates against the cell-identity markers in Figure 5 before interpreting any of them as regulation.")) +
    theme_pub(grid = "x") +
    theme(axis.text.y = element_text(size = rel(0.82), hjust = 0),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = rel(1.15), face = "bold"),
          strip.background = element_rect(fill = "#eeeef2", colour = NA),
          legend.position = "top", legend.justification = c(0, 0.5),
          legend.location = "plot", legend.direction = "horizontal",
          legend.key.size = unit(10, "pt"),
          plot.margin = margin(13, 16, 11, 13))
}

# ===========================================================================
# FIGURE 4 -- boxplots by group
# ===========================================================================

get_vst_group <- function(gse) {
  rds <- file.path(PATHS$processed, paste0(gse, "_pipeline_result.rds"))
  if (!file.exists(rds)) return(NULL)
  r <- readRDS(rds)
  # prefer the grouping stored by the pipeline -- always correct, including for
  # replicate-collapsed datasets
  if (!is.null(r$group)) {
    grp <- factor(r$group, levels = c("Control", "Case"))
    names(grp) <- names(r$group)
    common <- intersect(colnames(r$vst), names(grp))
    return(list(vst = r$vst[, common, drop = FALSE],
                grp = droplevels(factor(grp[common]))))
  }
  dds <- attr(r$de, "dds")
  if (!is.null(dds)) {
    grp <- SummarizedExperiment::colData(dds)$group; names(grp) <- colnames(dds)
  } else {
    md <- utils::read.csv(file.path(PATHS$metadata, paste0(gse, "_metadata.csv")),
                          stringsAsFactors = FALSE)
    if (is.null(md$sample_id)) md$sample_id <- md$gsm_id
    spec <- CONTRAST_SPEC[[gse]]
    if (!is.null(spec) && !is.null(md[[spec$col]])) {
      v <- trimws(as.character(md[[spec$col]]))
      md <- md[v %in% c(spec$case, spec$control), ]
      grp <- factor(ifelse(trimws(as.character(md[[spec$col]])) == spec$control,
                           "Control","Case"), levels = c("Control","Case"))
    } else {
      grp <- assign_groups(md, r$group_col, verbose = FALSE)
    }
    names(grp) <- md$sample_id
  }
  common <- intersect(colnames(r$vst), names(grp))
  list(vst = r$vst[, common, drop = FALSE], grp = droplevels(factor(grp[common])))
}

fig4_boxplots <- function() {
  sel <- FEATURED_GSES
  eff <- load_effects()
  rows <- list()
  for (g in sel) {
    o <- get_vst_group(g); if (is.null(o)) next
    panel <- resolve_panel(GENES_OF_INTEREST, rownames(o$vst))
    for (i in seq_len(nrow(panel))) {
      if (!panel$found[i]) next
      e <- eff[eff$study == g & eff$symbol == panel$symbol[i], ]
      rows[[length(rows)+1]] <- data.frame(
        study = g, gene = panel$symbol[i],
        expr = as.numeric(o$vst[panel$id[i], ]),
        group = as.character(o$grp),
        padj = if (nrow(e)) e$padj[1] else NA_real_,
        lfc  = if (nrow(e)) e$log2FC_shrunk[1] else NA_real_,
        disease = DISEASE_LABEL[as.character(eff$disease[eff$study == g][1])],
        stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, rows)

  d$group <- factor(d$group, levels = c("Control", "Case"))
  d$gene  <- factor(d$gene, levels = GENES)

  # Direction, taken straight from the differential-expression result so the
  # fill can never disagree with the reported statistic.
  d$direction <- ifelse(is.na(d$padj) | d$padj >= ALPHA, "NS",
                        ifelse(d$lfc > 0, "Up", "Down"))
  # Control boxes stay neutral; only the CASE box carries the finding.
  d$fill_key <- ifelse(d$group == "Control", "Control", d$direction)
  d$fill_key <- factor(d$fill_key, levels = c("Control", "Down", "NS", "Up"))

  nn <- aggregate(expr ~ study + gene + group, d, length); names(nn)[4] <- "n"
  d <- merge(d, nn, by = c("study","gene","group"))
  d$xlab <- factor(sprintf("%s\n(n = %d)", d$group, d$n),
                   levels = unique(sprintf("%s\n(n = %d)", d$group, d$n)[
                     order(d$group, d$study)]))

  d$strip <- sprintf("%s  \u00b7  %s\n%s", d$gene, d$disease, d$study)
  lev <- unlist(lapply(GENES, function(gg)
    unique(d$strip[d$gene == gg])[order(match(
      d$study[d$gene == gg][!duplicated(d$strip[d$gene == gg])], sel))]))
  d$strip <- factor(d$strip, levels = lev)

  # p-value annotation: stars alongside, and bold when significant
  ann <- unique(d[, c("strip","padj","direction")])
  ann$star <- ifelse(is.na(ann$padj), "",
              ifelse(ann$padj < 0.001, "***",
              ifelse(ann$padj < 0.01,  "**",
              ifelse(ann$padj < ALPHA, "*", ""))))
  ann$lab <- ifelse(is.na(ann$padj), "not tested",
                    trimws(paste0("FDR ", fmt_p(ann$padj), "  ", ann$star)))
  ann$face <- ifelse(!is.na(ann$padj) & ann$padj < ALPHA, "bold", "plain")
  ann$col  <- ifelse(!is.na(ann$padj) & ann$padj < ALPHA,
                     INK$primary, INK$muted)
  ann$y <- sapply(ann$strip, function(s) max(d$expr[d$strip == s], na.rm = TRUE))

  ggplot(d, aes(xlab, expr)) +
    geom_boxplot(aes(fill = fill_key), width = 0.52, outlier.shape = NA,
                 linewidth = 0.42, colour = INK$secondary) +
    geom_jitter(aes(colour = PUB_BOX_PT[as.character(fill_key)]), width = 0.13, size = 1.15,
                alpha = 0.45, shape = 16, show.legend = FALSE) +
    geom_text(data = ann,
              aes(x = 1.5, y = y, label = lab, fontface = face, colour = col),
              inherit.aes = FALSE, size = 2.5, family = PUB_FONT,
              vjust = -1.15, show.legend = FALSE) +
    scale_colour_identity() +
    facet_wrap(~ strip, nrow = 2, scales = "free") +
    scale_fill_manual(
      values = PUB_BOX, name = NULL, drop = FALSE,
      labels = c(Control = "Control",
                 Down    = "Case — lower in disease",
                 NS      = "Case — not significant",
                 Up      = "Case — higher in disease"),
      guide = guide_legend(override.aes = list(colour = INK$secondary,
                                               linewidth = 0.42))) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5),
                       expand = expansion(mult = c(0.10, 0.26))) +
    labs(title = sprintf("%s expression by group", GENE_TXT),
         subtitle = sprintf("%d representative datasets \u00b7 case boxes are coloured by direction of change, matching Figures 1-3 and 5",
                            length(unique(d$study))),
         x = NULL, y = LAB_VST,
         caption = paste0(
           "Boxes show median and interquartile range; whiskers 1.5 x IQR; every sample plotted. Significance by Benjamini-Hochberg FDR: * < 0.05, ** < 0.01, *** < 0.001;\n",
           "significant values are set in bold. Colour of the case box is taken directly from the differential-expression result (DESeq2 with apeglm, or limma-trend with\n",
           "ashr for normalised series). Panels use independent y-scales, so box heights are comparable within a panel but not across panels.")) +
    theme_pub(grid = "y") +
    theme(strip.text = element_text(size = rel(0.82), lineheight = 1.22,
                                    face = "plain", colour = INK$primary,
                                    margin = margin(b = 5)),
          axis.text.x = element_text(size = rel(0.78), lineheight = 1.10),
          axis.text.y = element_text(size = rel(0.80)),
          panel.spacing.x = unit(13, "pt"),
          panel.spacing.y = unit(19, "pt"),
          legend.position = "top", legend.justification = c(0, 0.5),
          legend.location = "plot", legend.direction = "horizontal",
          legend.key.size = unit(10, "pt"),
          legend.margin = margin(b = 6),
          plot.margin = margin(13, 18, 11, 13))
}

# ===========================================================================
# FIGURE 5 -- marker heatmap
# ===========================================================================

fig5_heatmap <- function() {
  sel <- FEATURED_GSES
  # genes of interest first, then every cell-identity panel in MARKERS
  gene_order <- c(GENES, unlist(unname(MARKERS)))
  gene_order <- gene_order[!duplicated(gene_order)]
  eff <- load_effects()

  # Read the FULL DE tables, not the tissue-class marker subset. The subset
  # only stores the panel matching each dataset's compartment, which left the
  # endothelial rows blank for blood datasets even though those genes were
  # measured -- and those cells are exactly what the tissue/blood comparison
  # needs.
  rows <- list()
  for (g in sel) {
    f <- file.path(PATHS$de, paste0(g, "_DE_full.csv"))
    if (!file.exists(f)) next
    d <- data.table::fread(f, data.table = FALSE)
    d <- d[d$symbol %in% gene_order, c("symbol","log2FC_shrunk","padj")]
    d <- d[!duplicated(d$symbol), ]
    d$study <- g
    d$disease <- DISEASE_LABEL[as.character(eff$disease[eff$study == g][1])]
    d$tissue_class <- DATASETS$tissue_class[DATASETS$gse == g]
    rows[[g]] <- d
  }
  d <- do.call(rbind, rows)

  d$symbol <- factor(d$symbol, levels = rev(gene_order))
  panel_of <- function(sym) {
    for (nm in names(MARKERS)) if (sym %in% MARKERS[[nm]]) {
      return(sprintf("%s markers", tools::toTitleCase(nm)))
    }
    "Genes of interest"
  }
  d$grp <- vapply(as.character(d$symbol), panel_of, character(1))
  d$grp <- factor(d$grp, levels = c("Genes of interest",
                                    sprintf("%s markers",
                                            tools::toTitleCase(names(MARKERS)))))
  # The compartment header carries the insight, so it can never be clipped by a
  # panel edge and sits directly above the columns it describes.
  # Column headers name the compartment and the datasets in it, so the reader
  # can see at a glance which columns are being compared.
  lab_for <- function(tc) {
    st <- unique(d$study[d$tissue_class == tc])
    cond <- unique(as.character(d$disease[d$tissue_class == tc]))
    sprintf("%s   %s", toupper(tc), paste(cond, collapse = " \u00b7 "))
  }
  LAB_TISSUE <- lab_for("tissue")
  LAB_BLOOD  <- lab_for("blood")
  d$compartment <- factor(ifelse(d$tissue_class == "blood", LAB_BLOOD, LAB_TISSUE),
                          levels = c(LAB_TISSUE, LAB_BLOOD))

  n_by <- unique(eff[, c("study","n_control","n_case")])
  d <- merge(d, n_by, by = "study")
  d$col_lab <- sprintf("%s\n%s\n%d v %d", d$disease, d$study, d$n_control, d$n_case)
  d$col_lab <- factor(d$col_lab, levels = unique(d$col_lab[order(match(d$study, sel))]))

  d$star <- ifelse(is.na(d$padj), "",
            ifelse(d$padj < 0.001,"***", ifelse(d$padj < 0.01,"**",
            ifelse(d$padj < ALPHA,"*",""))))

  # FIXED colour cap. A data-driven cap (previously the 90th percentile of
  # |log2FC|) drifts whenever the set of plotted cells changes -- adding the
  # off-panel marker cells moved it from 1.3 to 1.2 with no change in method.
  # Pinning it keeps the scale stable and comparable across revisions.
  lim <- 1.5
  n_sq <- sum(abs(d$log2FC_shrunk) > lim, na.rm = TRUE)

  ggplot(d, aes(col_lab, symbol, fill = log2FC_shrunk)) +
    geom_tile(colour = INK$surface, linewidth = 1.6) +
    geom_text(aes(label = star), size = 2.9, colour = INK$primary,
              fontface = "bold", vjust = 0.5, family = PUB_FONT) +
    facet_grid(grp ~ compartment, scales = "free", space = "free", switch = "y") +
    # emphasise the two genes of interest with taller cells than the marker rows
    ggh4x::force_panelsizes(rows = c(1.40, 1, 1), respect = FALSE) +
    scale_fill_gradient2(low = PUB_DIVERGING[["low"]], mid = PUB_DIVERGING[["mid"]],
                         high = PUB_DIVERGING[["high"]], midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish,
                         breaks = c(-lim, -lim/2, 0, lim/2, lim),
                         labels = c(paste0("\u2264 -", lim), sprintf("%.2g", -lim/2),
                                    "0", sprintf("%.2g", lim/2), paste0("\u2265 ", lim)),
                         name = expression(atop("Shrunken", log[2]*" fold change")),
                         guide = guide_colourbar(barwidth = unit(9, "pt"),
                                                 barheight = unit(92, "pt"),
                                                 title.position = "top",
                                                 ticks.colour = INK$surface,
                                                 ticks.linewidth = 1.1,
                                                 frame.colour = INK$faint,
                                                 frame.linewidth = 0.4)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
    labs(title = sprintf("%s read against cell-identity markers", GENE_TXT),
         subtitle = paste0(
           "Case-vs-control effect size per gene. A gene of interest moving TOGETHER with the markers below it indicates a shift in cell composition;\n",
           "the same gene moving while the markers stay flat indicates regulation within the cells. This panel is what licenses that distinction."),
         x = NULL, y = NULL,
         caption = paste0(
           sprintf("Cell colour encodes the shrunken log\u2082 fold change for disease versus control \u2014 blue lower in disease, red higher, grey unchanged \u2014 on a scale capped at \u00b1%.1f;\n", lim),
           sprintf("%s exceed the cap and %s drawn at the extreme. Stars denote Benjamini\u2013Hochberg FDR (* < 0.05, ** < 0.01, *** < 0.001), and blank cells mark genes falling\n",
                   ifelse(n_sq == 1, "One cell", sprintf("%d cells", n_sq)),
                   ifelse(n_sq == 1, "is", "are")),
           "below that dataset's expression filter. Column labels give n control v n case.")) +
    theme_pub(grid = "none") +
    theme(
      axis.text.x = element_text(size = rel(0.74), lineheight = 1.10,
                                 colour = INK$secondary),
      axis.text.y = element_text(face = "bold", size = rel(0.94),
                                 colour = INK$primary),
      axis.line = element_blank(), axis.ticks = element_blank(),
      # row-group bands: a filled strip on the left names each block
      strip.text.y.left = element_text(angle = 90, hjust = 0.5, face = "bold",
                                       size = rel(0.80), colour = INK$primary,
                                       margin = margin(r = 5, l = 3)),
      strip.background.y = element_rect(fill = "#f0f0f3", colour = NA),
      # compartment header band across the columns
      strip.text.x = element_text(face = "bold", size = rel(0.82),
                                  colour = INK$primary, lineheight = 1.5,
                                  margin = margin(t = 6, b = 6)),
      strip.background.x = element_rect(fill = "#e8e8ee", colour = NA),
      strip.placement = "outside",
      panel.spacing.y = unit(7, "pt"),
      panel.spacing.x = unit(9, "pt"),
      legend.position = "right", legend.direction = "vertical",
      legend.title = element_text(size = rel(0.82), colour = INK$secondary,
                                  lineheight = 1.1),
      plot.margin = margin(13, 16, 11, 13))
}

# ===========================================================================
# FIGURE 6 -- PCA / batch structure
# ===========================================================================

# Ellipse layer for Figure 6, switchable so the same figure can be rendered with
# shaded hulls, outline-only hulls, none at all, or hulls on the BEFORE panel
# only (where they carry the "studies separate" point without cluttering the
# corrected panel, whose whole message is that the hulls now overlap).
ell_layer <- function(d, mode) {
  if (mode == "none") return(NULL)
  base <- if (mode == "before_only") d[d$panel == levels(d$panel)[1], ] else d
  if (mode == "outline") {
    stat_ellipse(data = base, aes(colour = study), geom = "path",
                 type = "norm", level = 0.75, linewidth = 0.6)
  } else {
    stat_ellipse(data = base, aes(fill = study, colour = study),
                 geom = "polygon", type = "norm", level = 0.75,
                 alpha = 0.11, linewidth = 0.45)
  }
}

ell_note <- function(mode) switch(mode,
  filled      = "Shaded hulls are 75% normal-probability ellipses per study. ",
  outline     = "Outlines are 75% normal-probability ellipses per study. ",
  none        = "",
  before_only = "Shaded hulls in the left panel are 75% normal-probability ellipses per study; the corrected panel is shown as points alone, since its message is the absence of study separation. ")

fig6_pca <- function(ellipse = c("filled","outline","none","before_only")) {
  ellipse <- match.arg(ellipse)
  atlas <- file.path(PATHS$processed, "coexpression_atlas.rds")
  if (!file.exists(atlas)) return(NULL)
  a <- readRDS(atlas)
  if (is.null(a$corrected)) return(NULL)

  # A COMMON BASIS for both panels. Fitting a separate PCA to each matrix gives
  # two different latent spaces, so "PC1 before" and "PC1 after" are not the
  # same axis and the panels cannot be compared directly even with identical
  # limits. Instead the basis is fitted once on the UNCORRECTED data and the
  # corrected data is projected onto it, so both panels share one coordinate
  # system and the collapse of between-study spread is read off directly.
  pc <- stats::prcomp(t(a$combined), center = TRUE, scale. = FALSE)
  vp <- 100 * pc$sdev^2 / sum(pc$sdev^2)

  proj <- function(m) {
    z <- scale(t(m), center = pc$center, scale = FALSE)
    as.data.frame(z %*% pc$rotation[, 1:2])
  }
  b <- proj(a$combined); af <- proj(a$corrected)
  names(b) <- names(af) <- c("PC1", "PC2")

  d <- rbind(
    data.frame(b,  panel = "BEFORE correction", a$meta[, c("study","cell_type")]),
    data.frame(af, panel = "AFTER ComBat  (cell type protected)",
               a$meta[, c("study","cell_type")]))
  d$panel <- factor(d$panel, levels = c("BEFORE correction",
                                        "AFTER ComBat  (cell type protected)"))

  # identical limits on both panels, padded so ellipses are not clipped
  xr <- range(d$PC1); yr <- range(d$PC2)
  pad <- function(r) r + c(-1, 1) * diff(r) * 0.10

  vsum <- utils::read.csv(file.path(PATHS$coexpr,
                                    "atlas_variance_before_after_combat.csv"))
  gv <- function(cov, when) vsum[[paste0("total_var_explained_pct_", when)]][vsum$covariate == cov]

  # Both categorical channels are sized from the DATA, never hardcoded: the
  # number of studies and cell types comes from the user's registry, and a
  # fixed-length manual scale errors out the moment a fifth level appears.
  n_study <- length(unique(d$study))
  n_ctype <- length(unique(d$cell_type))
  sh <- pub_cat_shapes(n_ctype)

  ggplot(d, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = INK$rule, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = INK$rule, linewidth = 0.3) +
    ell_layer(d, ellipse) +
    (if (is.null(sh))
       geom_point(aes(colour = study), shape = 16, size = 3.0,
                  alpha = 0.85, stroke = 0.8)
     else
       geom_point(aes(colour = study, shape = cell_type), size = 3.0,
                  alpha = 0.85, stroke = 0.8)) +
    facet_wrap(~ panel, ncol = 2) +
    coord_cartesian(xlim = pad(xr), ylim = pad(yr)) +
    scale_colour_manual(values = pub_cat_colours(n_study),
                        name = "Study", aesthetics = c("colour","fill")) +
    (if (is.null(sh)) NULL else scale_shape_manual(values = sh, name = "Cell type")) +
    labs(
      title = "Batch correction of the baseline co-expression atlas",
      subtitle = sprintf(
        "%d samples · %d studies · %d cell types. Both panels share one PCA basis, fitted on the uncorrected data, so the axes and limits are identical.\nBy PC regression across the top 5 components, study batch accounts for %.1f%% of total variance before correction and %.1f%% after; cell type %.1f%% and %.1f%%.",
        nrow(a$meta), length(unique(a$meta$study)), length(unique(a$meta$cell_type)),
        gv("study","before"), gv("study","after"),
        gv("cell_type","before"), gv("cell_type","after")),
      x = sprintf("PC1 of the uncorrected basis (%.1f%% of variance in that fit)", vp[1]),
      y = sprintf("PC2 (%.1f%%)", vp[2]),
      caption = paste0(
        ell_note(ellipse),
        "Axis percentages describe the shared PCA basis; the subtitle percentages are a different quantity —\n",
        "variance attributable to each factor by PC regression — so the two sets of numbers are not expected to match. ComBat was applied with cell type protected via the\n",
        "model matrix, restricted to cell types present in >=2 studies and studies carrying >=2 cell types, because cell type is otherwise nested within study and the\n",
        "correction is not estimable. Cell-type variance falls alongside batch variance, which is stated as a limitation. Samples that stay separated in the corrected\n",
        "panel are inspected individually rather than trimmed, and are retained when they pass expression QC.")) +
    theme_pub(grid = "both") +
    theme(strip.text = element_text(size = rel(1.0), face = "bold"),
          strip.background = element_rect(fill = "#eeeef2", colour = NA),
          legend.position = "right", legend.box = "vertical",
          legend.key.size = unit(11, "pt"),
          panel.grid.major = element_line(colour = INK$rule, linewidth = 0.2),
          panel.spacing.x = unit(13, "pt"),
          plot.margin = margin(13, 16, 11, 13))
}

# ===========================================================================
# FIGURE 7 -- gene vs comparator correlation across the atlas studies
# ===========================================================================

fig7_correlation <- function() {
  cors <- utils::read.csv(file.path(PATHS$coexpr,
                                    "atlas_within_study_correlations.csv"))
  cors <- cors[cors$subgroup == "all", ]
  pool <- utils::read.csv(file.path(PATHS$coexpr, "atlas_pooled_correlation.csv"))

  atlas <- file.path(PATHS$processed, "coexpression_atlas.rds")
  pts <- NULL
  if (file.exists(atlas)) {
    a <- readRDS(atlas)
    m <- a$combined
    pan <- resolve_panel(c(G1, G2), rownames(m))
    if (all(pan$found)) {
      pts <- data.frame(gene_a = as.numeric(m[pan$id[1], ]),
                        gene_b = as.numeric(m[pan$id[2], ]),
                        study = a$meta$study, cell_type = a$meta$cell_type,
                        stringsAsFactors = FALSE)
    }
  }

  # left: per-study forest of r; right: pooled scatter across vascular beds
  # Ordered by SAMPLE SIZE, largest first. This is not cosmetic: it groups the
  # wide, uninformative intervals at the bottom where n is smallest, so the
  # relationship between precision and n is legible instead of looking like
  # disagreement between studies.
  f <- cors[order(cors$n), ]
  f$label <- sprintf("%s  (n = %d)", f$study, f$n)
  fr <- rbind(
    data.frame(label = f$label, est = f$r, lo = f$ci_low, hi = f$ci_high,
               kind = "Study", stringsAsFactors = FALSE),
    data.frame(label = sprintf("Pooled random-effects  (k = %d, n = %d)",
                               pool$k, sum(f$n)),
               est = pool$pooled_r, lo = pool$ci_low_r, hi = pool$ci_high_r,
               kind = "Pooled", stringsAsFactors = FALSE))
  fr$label <- factor(fr$label, levels = rev(fr$label))

  pA <- ggplot(fr, aes(est, label)) +
    geom_vline(xintercept = 0, linetype = "22", colour = INK$muted, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = kind), height = 0, linewidth = 0.6) +
    geom_point(aes(colour = kind, shape = kind, size = kind)) +
    scale_colour_manual(values = c(Study = PUB_NEUTRAL, Pooled = PUB_ACCENT), guide = "none") +
    scale_shape_manual(values = c(Study = 16, Pooled = 18), guide = "none") +
    scale_size_manual(values = c(Study = 2.4, Pooled = 4.4), guide = "none") +
    scale_x_continuous(limits = c(-0.5, 1.02), breaks = seq(-0.5, 1, 0.5)) +
    labs(subtitle = sprintf("PRIMARY RESULT \u2014 correlated within each study, then pooled\nr = %.3f [%.3f, %.3f], k = %d, total n = %d, I\u00b2 = %.0f%%",
                            pool$pooled_r, pool$ci_low_r, pool$ci_high_r,
                            pool$k, sum(f$n), pool$I2),
         x = LAB_R, y = NULL) +
    theme_pub(grid = "x") +
    theme(axis.text.y = element_text(size = rel(0.84), hjust = 0),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  if (is.null(pts)) return(pA)

  ct <- stats::cor.test(pts$gene_a, pts$gene_b)
  # Same rule as Figure 6: both categorical scales are sized from the data.
  n_ctype <- length(unique(pts$cell_type))
  n_study <- length(unique(pts$study))
  sh <- pub_cat_shapes(n_study)

  pB <- ggplot(pts, aes(gene_a, gene_b)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = PUB_ACCENT,
                fill = PUB_ACCENT, alpha = 0.10, linewidth = 0.65) +
    (if (is.null(sh))
       geom_point(aes(colour = cell_type), shape = 16, size = 2.3, alpha = 0.9)
     else
       geom_point(aes(colour = cell_type, shape = study), size = 2.3, alpha = 0.9)) +
    scale_colour_manual(values = pub_cat_colours(n_ctype), name = "Cell type") +
    (if (is.null(sh)) NULL else scale_shape_manual(values = sh, name = "Study")) +
    labs(subtitle = sprintf("SUPPORTING VIEW \u2014 batch-corrected integrated atlas\nr = %.3f [%.3f, %.3f], n = %d samples from %d studies",
                            ct$estimate, ct$conf.int[1], ct$conf.int[2],
                            nrow(pts), length(unique(pts$study))),
         x = paste0(G1, " \u2014 ", LAB_VST),
         y = paste0(G2, " \u2014 ", LAB_VST)) +
    theme_pub(grid = "both") +
    theme(legend.position = "right", legend.box = "vertical",
          panel.grid.major = element_line(colour = INK$rule, linewidth = 0.2))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.15)) +
      patchwork::plot_annotation(
        title = sprintf("%s co-expression across baseline samples", GENE_TXT),
        # Every number here is read from the data, never written in by hand --
        # a caption that hardcodes a result silently lies the first time the
        # analysis is re-run on anything else.
        caption = paste0(
          sprintf("THE HEADLINE ESTIMATE IS THE LEFT PANEL: r = %.3f [%.3f, %.3f], pooled across all %d studies by random-effects meta-analysis (metafor, REML, Fisher z).\n",
                  pool$pooled_r, pool$ci_low_r, pool$ci_high_r, pool$k),
          sprintf("Correlations are computed WITHIN each study and only then pooled, so no between-study batch structure can inflate them. The right panel (r = %.3f) pools\n",
                  unname(ct$estimate)),
          sprintf("samples across the %d studies that survived the ComBat restriction; it shows the corrected atlas and is a supporting visualisation, NOT the claimed estimate\n",
                  length(unique(pts$study))),
          sprintf("\u2014 pooling samples before correlating is exactly what the within-study design avoids. Small studies carry wide intervals (smallest here n = %d); the\n",
                  min(f$n)),
          "random-effects model weights them accordingly. Baseline (untreated) samples only \u2014 stimulated, shear-stressed and cytokine-treated samples were excluded.\n",
          "Shaded band is the 95% confidence interval of the linear fit.\n",
          "Cell-type colours are a categorical scale; this figure contains no case/control contrast, so the direction scheme of Figures 1\u20135 does not apply here."),
        theme = theme_pub())
  } else pA
}

# ---------------------------------------------------------------------------

if (WHICH %in% c(0, 2)) { log_step("Figure 2 -- volcano")
  save_pub(fig2_volcano(), "Figure2_volcano_panels", width = 9.6, height = 6.6) }
if (WHICH %in% c(0, 3)) { log_step("Figure 3 -- forest")
  save_pub(fig3_forest(), "Figure3_forest_meta_analysis", width = 10.0, height = 6.4) }
if (WHICH %in% c(0, 4)) { log_step("Figure 4 -- boxplots")
  save_pub(fig4_boxplots(), "Figure4_boxplots_by_group", width = 10.2, height = 5.8) }
if (WHICH %in% c(0, 5)) { log_step("Figure 5 -- heatmap")
  save_pub(fig5_heatmap(), "Figure5_heatmap_markers", width = 8.8, height = 5.6) }
if (WHICH %in% c(0, 6)) {
  log_step("Figure 6 -- PCA")
  # Approved variant: shaded study hulls on the uncorrected panel, where they
  # make the disjoint study territories unmissable; bare points on the corrected
  # panel, whose message is the ABSENCE of study separation -- drawing four
  # overlapping hulls to depict "no structure" adds ink without adding meaning.
  # The asymmetry is disclosed in the caption rather than left implicit.
  p <- fig6_pca("before_only")
  if (!is.null(p)) save_pub(p, "Figure6_PCA_batch_correction",
                            width = 10.0, height = 5.0)
}
if (WHICH %in% c(0, 7)) { log_step("Figure 7 -- correlation")
  save_pub(fig7_correlation(), "Figure7_gene_pair_correlation_atlas",
           width = 11.0, height = 5.4) }
