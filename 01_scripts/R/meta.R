# ---------------------------------------------------------------------------
# R/meta.R -- formal random-effects meta-analysis (METHODS S4).
#
# Two pooling jobs:
#   1. within-study gene-pair correlations across the atlas studies
#      -> Fisher z transform, rma(), back-transform to r
#   2. case/control effect sizes (shrunken log2FC + SE) across studies
#      -> rma() on the log2FC scale, feeding the forest plot
#
# Both report tau^2, I^2, Q and the pooled estimate with CI. Simple averaging
# is explicitly not acceptable.
# ---------------------------------------------------------------------------

# ---- within-study correlation ---------------------------------------------

# Computes the gene-pair correlation WITHIN one study. Pooling happens later --
# never correlate across pooled raw samples, which would manufacture
# correlation out of between-study batch structure (docs/METHODS.md S8.3).
study_correlation <- function(mat_vst, study_id,
                              gene_a = GENES_OF_INTEREST[1],
                              gene_b = GENES_OF_INTEREST[2], method = "pearson",
                              group_by = NULL, meta = NULL) {
  panel <- resolve_panel(c(gene_a, gene_b), rownames(mat_vst))
  if (!all(panel$found)) {
    log_warn("[%s] %s missing -- cannot correlate", study_id,
             paste(panel$symbol[!panel$found], collapse = ", "))
    return(NULL)
  }
  a <- as.numeric(mat_vst[panel$id[1], ])
  b <- as.numeric(mat_vst[panel$id[2], ])

  make_row <- function(x, y, label, n_lab) {
    if (length(x) < 4 || stats::sd(x) == 0 || stats::sd(y) == 0) return(NULL)
    ct <- stats::cor.test(x, y, method = method)
    data.frame(
      study = study_id, subgroup = label, n = length(x),
      r = unname(ct$estimate),
      ci_low  = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
      ci_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
      p_value = ct$p.value, method = method,
      gene_a = gene_a, gene_b = gene_b,
      stringsAsFactors = FALSE)
  }

  if (is.null(group_by) || is.null(meta) || is.null(meta[[group_by]])) {
    return(make_row(a, b, "all", length(a)))
  }
  # per-cell-type correlations where the study covers several vascular beds
  g <- as.character(meta[[group_by]])
  rows <- c(list(make_row(a, b, "all", length(a))),
            lapply(unique(g), function(lv) {
              idx <- which(g == lv)
              make_row(a[idx], b[idx], lv, length(idx))
            }))
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ---- pooling correlations -------------------------------------------------

# Fisher z-transform, random-effects REML, back-transform. Returns both the
# model object (for forest plots) and a tidy summary row.
pool_correlations <- function(cor_table, subgroup = "all") {
  require_pkg("metafor", "random-effects meta-analysis (METHODS S4)")

  dat <- cor_table[cor_table$subgroup == subgroup & !is.na(cor_table$r), ]
  dat <- dat[dat$n > 3, ]
  if (nrow(dat) < 2) {
    log_warn("only %d study/studies -- cannot pool", nrow(dat))
    return(NULL)
  }
  log_step("pooling %d studies (Fisher z, RE/REML)", nrow(dat))

  # escalc: ZCOR = Fisher r-to-z, with variance 1/(n-3)
  es <- metafor::escalc(measure = "ZCOR", ri = dat$r, ni = dat$n,
                        slab = dat$study)
  fit <- metafor::rma(yi = es$yi, vi = es$vi, method = "REML", slab = dat$study)

  back <- function(z) (exp(2 * z) - 1) / (exp(2 * z) + 1)   # z -> r
  summary_row <- data.frame(
    k = fit$k,
    pooled_z = as.numeric(fit$beta),
    pooled_r = back(as.numeric(fit$beta)),
    ci_low_r = back(fit$ci.lb),
    ci_high_r = back(fit$ci.ub),
    se_z = fit$se,
    p_value = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2,
    H2 = fit$H2,
    Q = fit$QE,
    Q_p = fit$QEp,
    model = "random-effects (REML), Fisher z",
    stringsAsFactors = FALSE)

  log_info("pooled r = %.3f [%.3f, %.3f], p = %.3g",
           summary_row$pooled_r, summary_row$ci_low_r,
           summary_row$ci_high_r, summary_row$p_value)
  log_info("heterogeneity: tau^2 = %.4f, I^2 = %.1f%%, Q(%d) = %.2f, p = %.3g",
           fit$tau2, fit$I2, fit$k - 1, fit$QE, fit$QEp)

  list(fit = fit, summary = summary_row, data = dat, escalc = es)
}

# ---- pooling effect sizes -------------------------------------------------

# Pools shrunken log2FC across studies for one gene. The SE used is the
# shrunken lfcSE, so the pooling is consistent with the reported effect size.
pool_effect_sizes <- function(effect_table, gene) {
  require_pkg("metafor", "random-effects meta-analysis (METHODS S4)")

  dat <- effect_table[effect_table$symbol == gene &
                        !is.na(effect_table$log2FC_shrunk) &
                        !is.na(effect_table$lfcSE_shrunk) &
                        effect_table$lfcSE_shrunk > 0, ]
  if (nrow(dat) < 2) {
    log_warn("gene %s present in only %d study/studies -- cannot pool",
             gene, nrow(dat))
    return(NULL)
  }
  log_step("pooling %s across %d studies (RE/REML)", gene, nrow(dat))

  fit <- metafor::rma(yi = dat$log2FC_shrunk, sei = dat$lfcSE_shrunk,
                      method = "REML", slab = dat$study)

  summary_row <- data.frame(
    gene = gene, k = fit$k,
    pooled_log2FC = as.numeric(fit$beta),
    ci_low = fit$ci.lb, ci_high = fit$ci.ub,
    se = fit$se, z = as.numeric(fit$zval), p_value = fit$pval,
    tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE, Q_p = fit$QEp,
    model = "random-effects (REML) on shrunken log2FC",
    stringsAsFactors = FALSE)

  log_info("%s: pooled log2FC = %.3f [%.3f, %.3f], p = %.3g; I^2 = %.1f%%",
           gene, summary_row$pooled_log2FC, fit$ci.lb, fit$ci.ub,
           fit$pval, fit$I2)

  list(fit = fit, summary = summary_row, data = dat)
}

# ---- reporting ------------------------------------------------------------

# One-line heterogeneity statement for figure captions and the written report.
heterogeneity_text <- function(pooled) {
  s <- pooled$summary
  sprintf("Random-effects (REML), k = %d; tau^2 = %.4f, I^2 = %.1f%%, Q = %.2f (p = %.3g)",
          s$k, s$tau2, s$I2, s$Q, s$Q_p)
}
