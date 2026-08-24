#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 61_celltype_atlas_figure.R -- Figures 8+ : "which cell populations does the
# gene track with?"
#
#   Panel A (bars)  pooled random-effects correlation per cell type
#   Panel B (grid)  cohort x cell type, the per-cohort correlations behind it
#
# One figure per gene in GENES_OF_INTEREST, numbered from 8 so it continues the
# set produced by 50_publication_figures.R.
#
# Reads only 60_celltype_atlas.R's output; does no statistics of its own, so
# the figure and the tables can never disagree.
#
# ENCODING -- deliberately chosen so nothing rests on colour alone
#   colour        cell lineage family, also written on the x axis
#   grey          tested but not significant after BH  (both panels)
#   solid bar     pooled association significant at BH-FDR < ALPHA
#   circle / down-triangle  positive / negative correlation in that cohort
#   point area    |r|  (area, not radius -- radius exaggerates by the square)
#   nothing drawn xCell could not score that cell type in that cohort, which
#                 means "not measurable here", never "no association"
#
# Run from the project root:
#   Rscript 01_scripts/61_celltype_atlas_figure.R 2>&1 \
#     | tee 08_logs/61_celltype_atlas_figure.log
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(grid)
})

source("01_scripts/00_config.R")
source("01_scripts/00_theme_publication.R")
for (f in c("utils")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

SCRIPT <- "61_celltype_atlas_figure"
set.seed(SEED)

# save_pub() writes 600 dpi PNG + vector PDF into 07_final_figures/, the same
# saver the rest of the final figure set uses.

# Must match MIN_K in 60_celltype_atlas.R -- it is quoted in the caption.
MIN_K <- 3

# ---- palette ---------------------------------------------------------------
# Okabe-Ito five. Validated with the project's palette checker at --pairs all:
# lightness band PASS, chroma floor PASS, normal-vision floor PASS (worst
# adjacent dE 15.6), CVD worst-pair dE 7.6 -- inside the 6-8 band, which is
# permitted only alongside secondary encoding, and here every cell type is
# named on the axis and the families sit in contiguous labelled blocks.
FAMILY_COLS <- c(
  "Structural & parenchymal"       = "#0072B2",
  "Endothelial"                    = "#CC79A7",
  "Myeloid"                        = "#D55E00",
  "Lymphoid"                       = "#009E73",
  "HSC, progenitor & Mk/erythroid" = "#E69F00",
  # Only appears if the deconvolution vocabulary grew past the family map in
  # 60_celltype_atlas.R. The audit above covers the validated five; a sixth hue
  # is a fallback, not part of the checked palette.
  "Other"                          = "#56B4E9")
NS_GREY <- "#b9b9c0"
NS_LAB  <- "Tested, not significant"

# ---- data ------------------------------------------------------------------

STUDY_CSV  <- file.path(PATHS$cellcomp, "celltype_atlas_per_study.csv")
POOLED_CSV <- file.path(PATHS$results,  "meta_celltype_atlas_pooled.csv")
for (f in c(STUDY_CSV, POOLED_CSV)) {
  if (!file.exists(f)) {
    stop("Missing ", f, ".\nRun 01_scripts/60_celltype_atlas.R first.")
  }
}
study <- fread(STUDY_CSV)
pool  <- fread(POOLED_CSV)

# Family blocks contiguous, cell types alphabetical inside each block. Both the
# families and the cell types come from whatever 60_celltype_atlas.R actually
# scored, so a different deconvolution vocabulary changes the axis rather than
# breaking the script.
FAM_LEVELS  <- intersect(names(FAMILY_COLS), unique(study$family))
CELL_ORDER  <- unique(study[order(match(family, FAM_LEVELS), cell_type)]$cell_type)
stopifnot(length(CELL_ORDER) > 0, !anyDuplicated(CELL_ORDER))

# Cohort rows: tissue cohorts first, then blood, matching how Figures 1-5 order
# the set and keeping the two tissue classes readable as blocks.
meta_rows <- unique(study[, .(study, condition, tissue_class, n)])
setorder(meta_rows, -tissue_class, condition, study)
meta_rows[, label := sprintf("%s  ·  %s  (n = %d)", condition, study, n)]
ROW_ORDER <- rev(meta_rows$label)   # ggplot draws discrete y bottom-up

study <- merge(study, meta_rows[, .(study, label)], by = "study")

SIG <- ALPHA

# ---- one gene --------------------------------------------------------------

build_gene <- function(the_gene, fig_no) {

  gene <- the_gene   # plain local: `..gene` is for column selection, not i-filters

  p <- pool[gene == the_gene]
  p[, cell_type := factor(cell_type, levels = CELL_ORDER)]
  p[, family    := factor(family,    levels = FAM_LEVELS)]
  p[, sig := !is.na(p_adj) & p_adj < SIG]
  # Bars are filled by family when the pooled association survives BH, and grey
  # when it does not, so "grey = not significant" reads the same in both panels.
  p[, bar_fill := fifelse(sig, as.character(family), NS_LAB)]

  s <- study[gene == the_gene & reliable == TRUE & !is.na(r)]
  s[, cell_type := factor(cell_type, levels = CELL_ORDER)]
  s[, label     := factor(label,     levels = ROW_ORDER)]
  s[, sig := !is.na(p_adj) & p_adj < SIG]
  s[, fam_fill := fifelse(sig, as.character(family), NS_LAB)]
  s[, fam_fill := factor(fam_fill, levels = c(FAM_LEVELS, NS_LAB))]
  s[, absr := abs(r)]
  s[, direction := fifelse(r >= 0, "Positive (r > 0)", "Negative (r < 0)")]
  s[, direction := factor(direction,
                          levels = c("Positive (r > 0)", "Negative (r < 0)"))]

  stopifnot(!anyNA(s$cell_type), !anyNA(s$label), !anyNA(s$fam_fill))

  # A cohort where the gene never passed the expression filter contributes no
  # mark at all. Keeping the empty row reads as a rendering fault and
  # contradicts the cohort count in the title, so the row is dropped and the
  # omission is counted in the caption instead.
  rows_here <- meta_rows[label %in% unique(s$label)]
  dropped   <- meta_rows[!label %in% unique(s$label)]
  row_levels <- rev(rows_here$label)
  s[, label := factor(as.character(label), levels = row_levels)]

  n_studies <- uniqueN(s$study)
  k_rng <- range(p$k)
  message(sprintf("[%s] %d cohorts; panel A %d cell types pooled (k %d-%d), %d significant; panel B %d marks (%d significant, %d of them negative)",
                  gene, n_studies, nrow(p), k_rng[1], k_rng[2], sum(p$sig),
                  nrow(s), sum(s$sig), sum(s$sig & s$r < 0)))

  fill_vals <- c(FAMILY_COLS, setNames(NS_GREY, NS_LAB))

  # Boundaries between family blocks, drawn as faint rules in both panels so the
  # blocks are readable by position and not only by hue.
  fam_of_x <- factor(sapply(CELL_ORDER, function(ct)
                       as.character(unique(study$family[study$cell_type == ct]))),
                     levels = FAM_LEVELS)
  breaks_x <- which(diff(as.integer(fam_of_x)) != 0) + 0.5

  # ---- test counts, derived rather than asserted ---------------------------
  # The x axis carries every cell type scored, but only those pooled from at
  # least MIN_K cohorts were tested, and BH ran across exactly those. Quoting
  # the axis count in the caption would misstate the correction denominator;
  # quoting the tested count without saying the rest are untested would
  # misstate the axis. Both numbers are reported, and the denominator is
  # re-derived here so the caption can never drift from the analysis.
  n_axis   <- length(CELL_ORDER)
  n_tested <- nrow(p)
  n_untested <- n_axis - n_tested
  stopifnot(n_tested > 0, n_untested >= 0)
  # BH must have been applied across the tested cell types for this gene only.
  bh_check <- p.adjust(p$p_value, method = "BH")
  stopifnot(isTRUE(all.equal(bh_check, p$p_adj, tolerance = 1e-10)))
  message(sprintf("[%s] BH denominator verified: %d tested of %d plotted (%d untested, k < %d)",
                  the_gene, n_tested, n_axis, n_untested, MIN_K))

  # One direct label per family, not the global top n: the strongest hits are
  # typically near-duplicate members of the SAME family (the lineage groups
  # below are deliberately close together), so a global top n stacks two labels
  # on the same spot and leaves the other families unlabelled.
  top <- p[sig == TRUE][order(-abs(pooled_r))][, .SD[1], by = family]
  # Park each label just beyond the far end of its whisker, on the side the bar
  # points, so it never sits on top of the interval it annotates.
  LAB_PAD <- 0.045
  top[, lab_y := fifelse(pooled_r >= 0, ci_high_r + LAB_PAD, ci_low_r - LAB_PAD)]
  # Axis limits are then taken from the LABELS, not the whiskers, so a label on
  # the tallest bar cannot be clipped at the panel edge.
  y_lo <- min(c(p$ci_low_r, top$lab_y), na.rm = TRUE) - 0.045
  y_hi <- max(c(p$ci_high_r, top$lab_y), na.rm = TRUE) + 0.045

  # ---- panel A -------------------------------------------------------------
  pA <- ggplot(p, aes(cell_type, pooled_r)) +
    geom_vline(xintercept = breaks_x, colour = INK$rule, linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = INK$secondary, linewidth = 0.35) +
    geom_col(aes(fill = bar_fill), width = 0.66,
             colour = INK$secondary, linewidth = 0.18) +
    geom_linerange(aes(ymin = ci_low_r, ymax = ci_high_r),
                   colour = INK$secondary, linewidth = 0.28) +
    geom_label(data = top, aes(y = lab_y, label = sprintf("%.2f", pooled_r)),
               size = 2.3, fontface = "bold", colour = INK$primary,
               family = PUB_FONT, fill = INK$surface, alpha = 0.9,
               linewidth = 0, label.padding = unit(0.9, "pt"),
               label.r = unit(1, "pt"), vjust = 0.5) +
    scale_fill_manual(values = fill_vals, limits = c(FAM_LEVELS, NS_LAB),
                      drop = FALSE, guide = "none") +
    scale_x_discrete(limits = CELL_ORDER, drop = FALSE) +
    scale_y_continuous(
      name = sprintf("Pooled correlation with %s\n(random effects, Fisher z, r)", gene),
      limits = c(y_lo, y_hi),
      breaks = seq(-1, 1, 0.25)) +
    theme_pub(base_size = 8.4, grid = "y") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin  = margin(2, 2, 1, 2)) +
    labs(
      title = sprintf("Figure %d.  %s expression tracks these cell populations across %d cohorts",
                      fig_no, gene, n_studies),
      subtitle = sprintf(
        "xCell deconvolution of every bulk cohort, correlated with %s expression within each cohort, then pooled (metafor, Fisher z, REML).\nAssociation across samples — NOT expression measured in each cell type. Bars grey where the pooled association is not significant at BH-FDR < %.2f.",
        gene, SIG))

  # ---- panel B -------------------------------------------------------------
  pB <- ggplot(s, aes(cell_type, label)) +
    geom_vline(xintercept = breaks_x, colour = INK$rule, linewidth = 0.3) +
    geom_point(aes(fill = fam_fill, size = absr, shape = direction),
               stroke = 0.16, colour = INK$secondary) +
    scale_fill_manual(values = fill_vals, limits = c(FAM_LEVELS, NS_LAB),
                      drop = FALSE, name = "Cell lineage family",
                      guide = guide_legend(
                        order = 1, ncol = 1,
                        override.aes = list(shape = 21, size = 2.6,
                                            linetype = 0, stroke = 0.16,
                                            colour = INK$secondary))) +
    # Circle up / triangle down: the sign of the correlation never depends on
    # telling two colours apart.
    scale_shape_manual(values = c("Positive (r > 0)" = 21,
                                  "Negative (r < 0)" = 25),
                       name = "Direction",
                       guide = guide_legend(
                         order = 3, ncol = 1,
                         override.aes = list(size = 2.6, fill = INK$muted,
                                             colour = INK$secondary))) +
    scale_size(name = "|r| within cohort", range = c(0.8, 5.6),
               limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8),
               guide = guide_legend(
                 order = 2, ncol = 1,
                 override.aes = list(shape = 21, fill = INK$muted,
                                     colour = INK$secondary, stroke = 0.16))) +
    scale_x_discrete(limits = CELL_ORDER, drop = FALSE) +
    scale_y_discrete(limits = row_levels, drop = FALSE) +
    theme_pub(base_size = 8.4, grid = "none") +
    theme(panel.grid.major = element_line(colour = "#eeeef1", linewidth = 0.2),
          axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.4),
          axis.text.y  = element_text(size = 6.8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line    = element_blank(),
          axis.ticks   = element_blank(),
          plot.margin  = margin(1, 2, 2, 2)) +
    labs(caption = paste(
      # Three lines, identical structure for both genes. Any cohort exclusion is
      # appended to line 1 -- the same slot in both figures -- rather than
      # appearing further down in one and not the other.
      sprintf("%d of %d cohorts, n = %d samples%s. Per-cohort Pearson r (%s expression vs. xCell cell-type score), pooled by random-effects meta-analysis (metafor, Fisher z, REML; k = %d-%d cohorts per cell type).",
              nrow(rows_here), nrow(meta_rows), sum(rows_here$n),
              if (nrow(dropped))
                sprintf("; %d cohort(s) excluded, %s below the expression filter there",
                        nrow(dropped), gene)
              else "",
              gene, k_rng[1], k_rng[2]),
      sprintf("%d cell types are plotted; %d had k >= 3 cohorts and were tested, and BH-FDR < %.2f was applied across those %d - the remaining %d are plotted but untested and carry no bar. Grey = not significant; dot area = within-cohort |r|; whiskers = 95%% CI; blank = xCell could not score that cell type there.",
              n_axis, n_tested, SIG, n_tested, n_untested),
      sprintf("xCell scores are bulk-deconvolution estimates, so r shows what expression co-varies with, not expression measured inside a cell type. Full statistics (tau^2, I^2, Q) in %s.",
              basename(POOLED_CSV)),
      sep = "\n"))

  # Echo the footnote to the log so its wording can be reviewed without
  # opening the image.
  message(sprintf("\n---- Figure %d (%s) footnote ----\n%s\n", fig_no, the_gene,
                  pB$labels$caption))

  pA / pB +
    plot_layout(heights = c(1.25, 1.5), guides = "collect") &
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.spacing.y = unit(8, "pt"),
          legend.justification = "top")
}

# ---- render ----------------------------------------------------------------
# Numbered from 8 so these continue the set from 50_publication_figures.R.
FIRST_FIG_NO <- 8

for (i in seq_along(GENES_OF_INTEREST)) {
  g <- GENES_OF_INTEREST[i]
  if (!nrow(pool[gene == g])) {
    log_warn("%s: nothing pooled -- no figure", g); next
  }
  fig_no <- FIRST_FIG_NO + i - 1
  log_step("Figure %d -- cell-type association atlas (%s)", fig_no, g)
  save_pub(build_gene(g, fig_no),
           sprintf("Figure%d_celltype_association_%s", fig_no, g),
           width = 15.0, height = 9.4)
}

log_session_info(SCRIPT)
log_step("COMPLETE: cell-type association figures")
