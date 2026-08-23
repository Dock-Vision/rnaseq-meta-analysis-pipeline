#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20_meta_analysis.R -- pools the per-dataset case/control results.
#
# Reads every *_effect_for_meta.csv written by 10_run_dataset.R, pools each
# gene of interest across studies with metafor random-effects (METHODS S4),
# and produces the cross-condition dot plot and the forest plots.
#
# Pooling is random-effects, never a simple average: studies differ in tissue,
# platform and cohort, so tau^2, I^2 and Q are part of the result rather than
# diagnostics to be skipped.
#
# Usage:  Rscript 01_scripts/20_meta_analysis.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2) })

source("01_scripts/00_config.R")
source("01_scripts/00_theme.R")
for (f in c("utils", "annotate", "meta", "figures")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

SCRIPT <- "20_meta_analysis"

# ---- gather per-study effects --------------------------------------------

log_step("collecting per-study effect sizes")
files <- list.files(PATHS$de, pattern = "_effect_for_meta\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No *_effect_for_meta.csv found in ", PATHS$de,
       ".\nRun 01_scripts/10_run_dataset.R on at least two datasets first.")
}
log_info("found %d per-study effect files", length(files))

effects <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
log_info("%d rows across %d studies, %d genes",
         nrow(effects), length(unique(effects$study)),
         length(unique(effects$symbol)))
save_result(effects, PATHS$results, "meta_all_study_effects")

# ---- pool each gene -------------------------------------------------------

pooled_summaries <- list()
for (gene in GENES_OF_INTEREST) {
  pooled <- pool_effect_sizes(effects, gene)
  if (is.null(pooled)) next
  pooled_summaries[[gene]] <- pooled$summary

  p <- fig_forest(pooled, gene, effect_label = "log2 fold change (shrunken, disease vs control)")
  save_fig(p, sprintf("meta_forest_%s", gene), width = 7, height = 4.5)

  append_progress("ALL (meta)", sprintf("Random-effects meta-analysis: %s", gene),
                  "DONE",
                  sprintf("k=%d studies; pooled log2FC=%.3f [%.3f, %.3f], p=%.3g; tau^2=%.4f, I^2=%.1f%%, Q=%.2f (p=%.3g).",
                          pooled$summary$k, pooled$summary$pooled_log2FC,
                          pooled$summary$ci_low, pooled$summary$ci_high,
                          pooled$summary$p_value, pooled$summary$tau2,
                          pooled$summary$I2, pooled$summary$Q, pooled$summary$Q_p),
                  sprintf("06_figures/png/meta_forest_%s.png", gene))
}

if (length(pooled_summaries)) {
  save_result(do.call(rbind, pooled_summaries), PATHS$results,
              "meta_pooled_effect_sizes")
}

# ---- cross-condition dot plot ---------------------------------------------

log_step("cross-condition dot plot")
p <- fig_dotplot(effects)
DOT_NAME <- sprintf("meta_dotplot_%s", paste(GENES_OF_INTEREST, collapse = "_"))
save_fig(p, DOT_NAME, width = 9, height = 6)

append_progress("ALL (meta)", "Cross-condition dot plot", "DONE",
                sprintf("Effect sizes for %s across %d studies.",
                        paste(GENES_OF_INTEREST, collapse = "/"),
                        length(unique(effects$study))),
                sprintf("06_figures/png/%s.png", DOT_NAME))

# ---- tissue-stratified subgroup analysis ---------------------------------

# When a gene moves in opposite directions in solid tissue and in blood,
# overall heterogeneity is extreme (I^2 above ~75-90%) and a single pooled
# estimate averages the two effects away into something uninterpretable. So
# always stratify by tissue class as well, and run a formal moderator test of
# whether the subgroups genuinely differ (METHODS S4) rather than asserting it
# from the two subgroup estimates by eye.

log_step("tissue-stratified subgroup meta-analysis")

effects$tissue_class <- DATASETS$tissue_class[match(effects$study, DATASETS$gse)]
log_info("tissue classes: %s",
         paste(sprintf("%s=%d", names(table(effects$tissue_class)),
                       table(effects$tissue_class)), collapse = ", "))

subgroup_rows <- list()
for (gene in GENES_OF_INTEREST) {
  sub <- effects[effects$symbol == gene & !is.na(effects$log2FC_shrunk) &
                   !is.na(effects$lfcSE_shrunk) & effects$lfcSE_shrunk > 0, ]
  if (nrow(sub) < 4) next

  for (tc in unique(sub$tissue_class)) {
    d <- sub[sub$tissue_class == tc, ]
    if (nrow(d) < 2) {
      log_info("%s / %s: only %d study, not pooled", gene, tc, nrow(d))
      next
    }
    fit <- metafor::rma(yi = d$log2FC_shrunk, sei = d$lfcSE_shrunk,
                        method = "REML", slab = d$study)
    subgroup_rows[[length(subgroup_rows) + 1]] <- data.frame(
      gene = gene, tissue_class = tc, k = fit$k,
      pooled_log2FC = as.numeric(fit$beta),
      ci_low = fit$ci.lb, ci_high = fit$ci.ub,
      p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
      Q = fit$QE, Q_p = fit$QEp, stringsAsFactors = FALSE)
    log_info("%s / %-6s: k=%d, pooled log2FC = %+.3f [%.3f, %.3f], p = %.3g, I^2 = %.1f%%",
             gene, tc, fit$k, as.numeric(fit$beta), fit$ci.lb, fit$ci.ub,
             fit$pval, fit$I2)
  }

  # formal test of whether tissue class explains the heterogeneity
  if (length(unique(sub$tissue_class)) > 1) {
    mod <- metafor::rma(yi = log2FC_shrunk, sei = lfcSE_shrunk,
                        mods = ~ factor(tissue_class), data = sub,
                        method = "REML")
    log_info("%s: moderator test for tissue class -- QM(df=%d) = %.2f, p = %.4g; residual I^2 = %.1f%%",
             gene, mod$m, mod$QM, mod$QMp, mod$I2)
    subgroup_rows[[length(subgroup_rows) + 1]] <- data.frame(
      gene = gene, tissue_class = "MODERATOR TEST", k = mod$k,
      pooled_log2FC = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = mod$QMp, tau2 = mod$tau2, I2 = mod$I2,
      Q = mod$QM, Q_p = mod$QMp, stringsAsFactors = FALSE)
  }
}

if (length(subgroup_rows)) {
  sg <- do.call(rbind, subgroup_rows)
  save_result(sg, PATHS$results, "meta_subgroup_by_tissue")

  # Report the stratified estimates, plus the overall I^2 actually observed --
  # never a remembered number.
  g1_sg <- sg[sg$gene == GENES_OF_INTEREST[1] &
              sg$tissue_class != "MODERATOR TEST", ]
  i2_txt <- if (length(pooled_summaries)) {
    ov <- do.call(rbind, pooled_summaries)
    paste(sprintf("%s I^2=%.1f%%", ov$gene, ov$I2), collapse = ", ")
  } else "not computed"
  append_progress("ALL (meta)", "Tissue-stratified subgroup meta-analysis", "DONE",
                  sprintf("Overall heterogeneity: %s. Stratified by tissue: %s. Moderator test p-values in the CSV.",
                          i2_txt,
                          paste(sprintf("%s %s k=%d %+.3f (p=%.3g)",
                                        g1_sg$gene, g1_sg$tissue_class, g1_sg$k,
                                        g1_sg$pooled_log2FC, g1_sg$p_value),
                                collapse = "; ")),
                  "05_results/meta_subgroup_by_tissue.csv")
}

log_session_info(SCRIPT)
log_step("COMPLETE: meta-analysis")
