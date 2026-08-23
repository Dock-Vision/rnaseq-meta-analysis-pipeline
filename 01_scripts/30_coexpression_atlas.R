#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 30_coexpression_atlas.R -- baseline co-expression atlas for the gene pair.
#
# Answers a different question from 10_run_dataset.R: not "does this gene
# change in disease" but "is its normal relationship to its comparator
# preserved". Implements docs/METHODS.md S8:
#
#   1. Untreated / baseline samples ONLY -- every stimulated condition is
#      excluded, and the exclusion is recorded explicitly rather than assumed.
#   2. ComBat batch correction with CELL TYPE PROTECTED via `mod`, so real
#      biological differences between cell types are not regressed away along
#      with the study batch.
#   3. The two genes correlated WITHIN each study first, then pooled with
#      metafor random-effects (METHODS S4). Pooling raw samples across studies
#      and correlating once manufactures correlation out of between-study batch
#      structure.
#
# Plus METHODS S3: PCA before/after correction and the % variance
# attributable to study batch vs cell type.
#
# Usage:  Rscript 01_scripts/30_coexpression_atlas.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery); library(DESeq2); library(sva); library(ggplot2)
  library(data.table); library(SummarizedExperiment); library(Biobase)
})

source("01_scripts/00_config.R")
source("01_scripts/00_theme.R")
for (f in c("utils", "annotate", "download", "prepare", "meta", "figures")) {
  source(file.path("01_scripts", "R", paste0(f, ".R")))
}

SCRIPT <- "30_coexpression_atlas"

# ---- baseline sample selection (METHODS S8.1) -------------------------------------

# Terms that mark a STIMULATED / PERTURBED sample. Any sample whose metadata
# matches is excluded from the atlas. Deliberately broad -- for this analysis a
# false exclusion is much cheaper than a treated sample contaminating a
# "healthy baseline" correlation.
STIMULATION_TERMS <- c(
  "tnf", "il-?1", "il1b", "il-?6", "tgf-?b", "tgfb", "ifn", "lps",
  "vegf", "hypoxi", "shear", "flow", "stretch", "oxldl", "ox-ldl",
  "sirna", "shrna", "knockdown", "knock-down", "overexpress", "transfect",
  "crispr", "treated", "treatment", "stimulat", "activat", "induc",
  "drug", "inhibitor", "agonist", "antagonist", "dose", "\\bhg\\b",
  "high glucose", "senescen", "irradiat", "infect"
)

# Terms that positively mark a baseline sample.
BASELINE_TERMS <- c("untreated", "control", "vehicle", "unstimulated",
                    "baseline", "static", "normoxi", "resting", "\\bnt\\b",
                    "scramble", "mock", "\\bdmso\\b",
                    # vehicle values that mark a control sample
                    "\\bh2o\\b", "\\bwater\\b", "\\bpbs\\b", "\\bsaline\\b",
                    "\\bnone\\b", "no treatment", "\\bna\\b")

select_baseline_samples <- function(meta, gse) {
  # Search the parsed characteristic VALUES only. The raw_characteristics_*
  # columns hold "key: value" strings, so including them let a FIELD NAME
  # trigger exclusion -- GSE164868's vehicle controls read "treatment: h2o"
  # and every one of its 12 samples was dropped on the word "treatment".
  txt_cols <- grep("^(title|source|char_)", colnames(meta), value = TRUE)
  txt <- apply(meta[, txt_cols, drop = FALSE], 1,
               function(r) tolower(paste(r, collapse = " | ")))
  # strip any residual "key:" prefixes for the same reason
  txt <- gsub("[a-z0-9_. ]+:\\s*", "", txt)

  # Raw text WITH keys intact, for dose-style encodings where the magnitude is
  # the perturbation. GSE294621 records shear stress as "Dynes:0" (static,
  # baseline) vs "Dynes:10" (10 dyn/cm2 applied) -- stripping the key hides it,
  # so a non-zero dose must be detected explicitly or shear-stressed samples
  # would silently enter a "healthy baseline" atlas.
  txt_raw <- apply(meta[, txt_cols, drop = FALSE], 1,
                   function(r) tolower(paste(r, collapse = " | ")))
  DOSE_KEYS <- c("dyne", "shear", "dose", "concentration", "\\bconc\\b",
                 "ng/ml", "ug/ml", "\\bnm\\b", "\\bum\\b", "\\bmm\\b")
  dose_applied <- rep(FALSE, length(txt_raw))
  for (k in DOSE_KEYS) {
    m <- regmatches(txt_raw,
                    regexpr(paste0(k, "[^0-9-]{0,12}(-?[0-9]+\\.?[0-9]*)"), txt_raw))
    vals <- suppressWarnings(as.numeric(
      sub(paste0(".*?(-?[0-9]+\\.?[0-9]*)$"), "\\1", m)))
    idx <- which(nzchar(m))
    if (length(idx)) {
      hit <- rep(FALSE, length(txt_raw))
      hit[idx] <- !is.na(vals) & vals > 0
      dose_applied <- dose_applied | hit
    }
  }
  if (any(dose_applied)) {
    log_info("[%s] %d sample(s) carry a non-zero dose/shear value -> treated as stimulated",
             gse, sum(dose_applied))
  }

  is_stim <- Reduce(`|`, lapply(STIMULATION_TERMS, function(p) grepl(p, txt))) |
    dose_applied
  is_base <- Reduce(`|`, lapply(BASELINE_TERMS, function(p) grepl(p, txt))) &
    !dose_applied   # an applied dose overrides any "control" wording

  # a sample explicitly marked baseline overrides an incidental keyword hit
  # (e.g. "untreated control for TNF experiment")
  keep <- is_base | !is_stim

  audit <- data.frame(
    gse = gse,
    gsm_id = meta$gsm_id,
    title = meta$title,
    matched_stimulation = is_stim,
    matched_baseline = is_base,
    dose_applied = dose_applied,
    kept = keep,
    metadata_text = txt,
    stringsAsFactors = FALSE)

  log_info("[%s] baseline selection: %d/%d kept (%d excluded as stimulated)",
           gse, sum(keep), nrow(meta), sum(!keep))
  list(keep = keep, audit = audit)
}

# ---- cell type assignment -------------------------------------------------
#
# GEO names the same cell type a dozen different ways across submitters, so the
# free-text fields are normalised to one label per sample. This matters twice
# over: ComBat needs the cell-type variable in `mod` to protect real biological
# differences, and the correlations are reported per cell type.
#
# EDIT ME for your own cell types. Each entry is
#   Label = c(<lowercase substrings that identify it>)
# and the FIRST matching entry wins, so put the more specific patterns first --
# "umbilical artery" must be tested before "artery", or every HUAEC sample is
# silently relabelled. Anything unmatched becomes "Other".
#
# The default panel below covers human endothelial vascular beds.
CELL_TYPE_PATTERNS <- list(
  HUVEC        = c("huvec", "umbilical vein"),
  HUAEC        = c("huaec", "umbilical artery"),   # distinct bed from HUVEC
  HPAEC        = c("hpaec", "pulmonary artery", "pulmonary arterial"),
  HAEC         = c("haec", "aortic"),
  HCAEC        = c("coronary"),
  CarotidEC    = c("carotid"),
  SaphenousEC  = c("saphenous"),
  EndocardialEC = c("endocardial", "endocardium"),
  HMVEC        = c("hmvec", "microvascular", "micro-vascular", "dermal microvascular"),
  HLEC         = c("hlec", "lymphatic"),
  BrainEC      = c("brain", "bmec", "cerebral", "hbmec"),
  RenalEC      = c("renal", "kidney", "glomerular"),
  HepaticEC    = c("hepatic", "liver", "sinusoidal"),
  DermalEC     = c("dermal", "skin")
)

assign_cell_type <- function(meta) {
  txt_cols <- grep("^(title|source|char_)", colnames(meta), value = TRUE)
  txt <- apply(meta[, txt_cols, drop = FALSE], 1,
               function(r) tolower(paste(r, collapse = " | ")))
  out <- rep(NA_character_, nrow(meta))
  for (nm in names(CELL_TYPE_PATTERNS)) {
    hit <- Reduce(`|`, lapply(CELL_TYPE_PATTERNS[[nm]], function(p) grepl(p, txt)))
    out[is.na(out) & hit] <- nm
  }
  out[is.na(out)] <- "Other"
  log_info("cell types: %s",
           paste(sprintf("%s=%d", names(table(out)), table(out)), collapse = ", "))
  out
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

ec_gses <- DATASETS$gse[DATASETS$arm == "atlas"]
log_step("healthy EC atlas: %d datasets", length(ec_gses))

studies <- list()
audits  <- list()

for (gse in ec_gses) {
  res <- tryCatch(fetch_gse(gse), error = function(e) {
    log_warn("[%s] fetch failed: %s", gse, conditionMessage(e)); NULL })
  if (is.null(res) || is.null(res$counts)) {
    append_progress(gse, "Fetch (healthy EC)", "BLOCKED",
                    "No count matrix obtainable from GEO.", "-")
    next
  }
  meta <- res$metadata
  if (is.null(meta$sample_id)) meta$sample_id <- meta$gsm_id

  sel <- select_baseline_samples(meta, gse)
  audits[[gse]] <- sel$audit
  if (sum(sel$keep) < 4) {
    log_warn("[%s] only %d baseline samples -- too few for correlation",
             gse, sum(sel$keep))
    append_progress(gse, "Baseline selection", "SKIPPED",
                    sprintf("Only %d untreated samples after excluding stimulated conditions (need >=4).",
                            sum(sel$keep)), "-")
    next
  }

  meta <- meta[sel$keep, , drop = FALSE]

  # align by shared IDs -- count columns are not guaranteed to carry GSM IDs
  common <- intersect(meta$sample_id, colnames(res$counts))
  if (length(common) < 4) {
    log_warn("[%s] only %d samples shared between counts and metadata -- skipping",
             gse, length(common))
    append_progress(gse, "Baseline selection (healthy EC)", "BLOCKED",
                    sprintf("Count columns could not be matched to GEO samples (%d shared IDs).",
                            length(common)), "-")
    next
  }
  meta <- meta[match(common, meta$sample_id), , drop = FALSE]
  counts <- res$counts[, common, drop = FALSE]
  meta$cell_type <- assign_cell_type(meta)
  meta$study <- gse

  if (!res$scale$is_raw_counts) {
    log_warn("[%s] %s -- using values as-is on a log scale for correlation",
             gse, res$scale$verdict)
    vsd <- log2(as.matrix(counts) + 1)
  } else {
    meta$group <- factor("Control")   # placeholder; no contrast in the atlas
    dds <- build_dds(counts, meta, design = ~ 1)
    vsd <- vst_matrix(dds, blind = TRUE)
  }

  studies[[gse]] <- list(vst = vsd, meta = meta)

  append_progress(gse, "Baseline selection + VST (healthy EC)", "DONE",
                  sprintf("%d/%d samples kept as untreated baseline; EC types: %s.",
                          nrow(meta), nrow(res$metadata),
                          paste(names(table(meta$cell_type)), collapse = "/")),
                  file.path("04_processed", paste0(gse, "_counts.csv")))
}

if (length(studies) == 0) stop("No usable healthy EC datasets.")

# audit trail for sample selection -- METHODS S8.1 requires this to be explicit
save_result(do.call(rbind, audits), PATHS$coexpr, "atlas_sample_selection_audit")

# ---- within-study correlations (METHODS S8.3) ------------------------------------

log_step("within-study %s-%s correlations", GENES_OF_INTEREST[1], GENES_OF_INTEREST[2])
cors <- do.call(rbind, lapply(names(studies), function(g) {
  study_correlation(studies[[g]]$vst, g,
                    GENES_OF_INTEREST[1], GENES_OF_INTEREST[2],
                    group_by = "cell_type", meta = studies[[g]]$meta)
}))
save_result(cors, PATHS$coexpr, "atlas_within_study_correlations")

for (g in names(studies)) {
  p <- fig_correlation(studies[[g]]$vst, studies[[g]]$meta, g,
                       GENES_OF_INTEREST[1], GENES_OF_INTEREST[2], "cell_type")
  if (!is.null(p)) save_fig(p, sprintf("coexpr_%s_%s_vs_%s", g, GENES_OF_INTEREST[1],
                                      GENES_OF_INTEREST[2]),
                            width = 6.5, height = 5)
}

# ---- pooling (METHODS S4) -------------------------------------------------------

pooled <- pool_correlations(cors, subgroup = "all")
if (!is.null(pooled)) {
  save_result(pooled$summary, PATHS$coexpr, "atlas_pooled_correlation")
  p <- fig_forest(pooled, sprintf("%s-%s correlation", GENES_OF_INTEREST[1],
                                  GENES_OF_INTEREST[2]),
                  effect_label = "Pearson correlation coefficient (r)")
  save_fig(p, sprintf("coexpr_forest_%s_%s", GENES_OF_INTEREST[1],
                      GENES_OF_INTEREST[2]), width = 7, height = 4.5)

  append_progress("Co-expression atlas",
                  "Random-effects pooling of the gene-pair correlation",
                  "DONE",
                  sprintf("k=%d studies; pooled r=%.3f [%.3f, %.3f], p=%.3g; tau^2=%.4f, I^2=%.1f%%, Q=%.2f (p=%.3g). Within-study first, then pooled -- never pooled raw samples.",
                          pooled$summary$k, pooled$summary$pooled_r,
                          pooled$summary$ci_low_r, pooled$summary$ci_high_r,
                          pooled$summary$p_value, pooled$summary$tau2,
                          pooled$summary$I2, pooled$summary$Q, pooled$summary$Q_p),
                  "05_results/coexpression/atlas_pooled_correlation.csv")
}

# ---- ComBat with cell type protected (METHODS S8.2) ------------------------------

log_step("ComBat batch correction (cell type protected)")

common <- Reduce(intersect, lapply(studies, function(s) rownames(s$vst)))
log_info("%s genes common to all %d studies", fmt_n(length(common)), length(studies))

if (length(common) > 1000) {
  combined <- do.call(cbind, lapply(studies, function(s) s$vst[common, , drop = FALSE]))
  cmeta <- do.call(rbind, lapply(studies, function(s)
    s$meta[, c("sample_id", "study", "cell_type")]))
  cmeta$sample_id <- colnames(combined)

  # Merging studies quantified on different references leaves genes that are
  # absent or non-finite in some of them. prcomp() and ComBat both fail on
  # those, and a zero-variance gene contributes nothing but breaks the SVD.
  finite_rows <- apply(combined, 1, function(r) all(is.finite(r)))
  var_rows <- apply(combined, 1, function(r) stats::sd(r) > 1e-8)
  keep_rows <- finite_rows & var_rows
  log_info("gene filter for integration: %s -> %s (dropped %s non-finite, %s zero-variance)",
           fmt_n(nrow(combined)), fmt_n(sum(keep_rows)),
           fmt_n(sum(!finite_rows)), fmt_n(sum(finite_rows & !var_rows)))
  combined <- combined[keep_rows, , drop = FALSE]
  if (nrow(combined) < 500) {
    log_warn("only %s usable genes after filtering -- integration unreliable",
             fmt_n(nrow(combined)))
  }

  # variance BEFORE correction
  pcv_before <- pc_variance_explained(combined, cmeta, c("study", "cell_type"))
  vsum_before <- summarise_variance(pcv_before)
  log_info("before ComBat: %s",
           paste(sprintf("%s=%.1f%%", vsum_before$covariate,
                         vsum_before$total_var_explained_pct), collapse = ", "))

  # docs/METHODS.md S8.2 requires ComBat with cell type protected via `mod`.
  # That is only estimable when each protected cell type appears in more than
  # one study; here most EC types are nested inside a single study (e.g.
  # DermalEC only in GSE92724), which makes cell type confounded with batch and
  # ComBat refuses. Restrict to the EC types that genuinely span studies rather
  # than dropping the protection, which would let ComBat regress away real
  # biological differences between vascular beds.
  # Estimability requires BOTH: every cell type in >=2 studies, and every study
  # carrying >=2 cell types. Dropping one can break the other, so prune
  # iteratively until the design is stable (or nothing is left).
  keep_s <- rep(TRUE, nrow(cmeta))
  repeat {
    t <- table(cmeta$cell_type[keep_s], cmeta$study[keep_s])
    if (nrow(t) == 0 || ncol(t) == 0) break
    bad_types   <- rownames(t)[rowSums(t > 0) < 2]
    bad_studies <- colnames(t)[colSums(t > 0) < 2]
    if (length(bad_types) == 0 && length(bad_studies) == 0) break
    if (length(bad_types)) {
      log_info("  dropping single-study cell type(s): %s",
               paste(bad_types, collapse = ", "))
      keep_s <- keep_s & !(cmeta$cell_type %in% bad_types)
    }
    if (length(bad_studies)) {
      log_info("  dropping study/studies with only one cell type: %s",
               paste(bad_studies, collapse = ", "))
      keep_s <- keep_s & !(cmeta$study %in% bad_studies)
    }
  }

  spanning <- unique(cmeta$cell_type[keep_s])
  dropped <- setdiff(unique(cmeta$cell_type), spanning)
  log_info("EC types retained for integration: %s",
           if (length(spanning)) paste(sort(spanning), collapse = ", ") else "NONE")
  if (length(dropped)) {
    log_warn("excluded from ComBat integration: %s",
             paste(sort(dropped), collapse = ", "))
  }

  ok_for_combat <- length(spanning) >= 2 && sum(keep_s) >= 6 &&
    length(unique(cmeta$study[keep_s])) >= 2

  if (!ok_for_combat) {
    log_warn("SKIPPING ComBat: cell type is confounded with study in this dataset collection.")
    log_warn("The within-study-then-pool correlation (the primary analysis) is unaffected.")
    append_progress("Healthy EC atlas", "ComBat integration", "BLOCKED",
                    sprintf("Cell type is nested within study (EC types spanning >=2 studies: %s). ComBat with cell type protected is not estimable; integrating without protection would regress away real vascular-bed differences. Within-study-then-pool correlation is unaffected and remains the primary result.",
                            if (length(spanning)) paste(spanning, collapse = "/") else "none"),
                    "-")
    corrected <- NULL
  } else {
    combined <- combined[, keep_s, drop = FALSE]
    cmeta <- cmeta[keep_s, , drop = FALSE]
    log_info("ComBat on %d samples, %d studies, %d cell types",
             ncol(combined), length(unique(cmeta$study)),
             length(unique(cmeta$cell_type)))
    # recompute the "before" picture on the restricted set so before/after are
    # comparable (METHODS S3)
    pcv_before <- pc_variance_explained(combined, cmeta, c("study", "cell_type"))
    vsum_before <- summarise_variance(pcv_before)

    mod <- stats::model.matrix(~ cell_type, data = cmeta)
    set.seed(SEED)
    corrected <- sva::ComBat(dat = combined, batch = cmeta$study, mod = mod,
                             par.prior = TRUE, prior.plots = FALSE)
  }

  if (is.null(corrected)) {
    log_info("no integrated atlas produced; skipping before/after PCA and variance tables")
  } else {
  pcv_after <- pc_variance_explained(corrected, cmeta, c("study", "cell_type"))
  vsum_after <- summarise_variance(pcv_after)
  log_info("after ComBat:  %s",
           paste(sprintf("%s=%.1f%%", vsum_after$covariate,
                         vsum_after$total_var_explained_pct), collapse = ", "))

  cmp <- merge(vsum_before, vsum_after, by = "covariate",
               suffixes = c("_before", "_after"))
  save_result(cmp, PATHS$coexpr, "atlas_variance_before_after_combat")
  save_result(pcv_before, PATHS$coexpr, "atlas_pc_variance_before")
  save_result(pcv_after, PATHS$coexpr, "atlas_pc_variance_after")

  # METHODS S3 -- PCA before and after, side by side
  p <- fig_pca_before_after(combined, corrected, cmeta, "Healthy EC atlas",
                            colour_col = "cell_type", batch_col = "study")
  save_fig(p, "coexpr_pca_before_after_combat", width = 12, height = 5)

  p2 <- fig_variance_explained(vsum_after, "Healthy EC atlas (after ComBat)")
  save_fig(p2, "coexpr_variance_explained_after", width = 6.5, height = 4)

  saveRDS(list(combined = combined, corrected = corrected, meta = cmeta),
          file.path(PATHS$processed, "coexpression_atlas.rds"))

  append_progress("Healthy EC atlas", "ComBat (cell type protected) + variance decomposition",
                  "DONE",
                  sprintf("%s common genes, %d samples, %d studies. Study-batch variance %.1f%% -> %.1f%%; cell-type variance %.1f%% -> %.1f%%.",
                          fmt_n(length(common)), ncol(combined), length(studies),
                          cmp$total_var_explained_pct_before[cmp$covariate == "study"],
                          cmp$total_var_explained_pct_after[cmp$covariate == "study"],
                          cmp$total_var_explained_pct_before[cmp$covariate == "cell_type"],
                          cmp$total_var_explained_pct_after[cmp$covariate == "cell_type"]),
                  "05_results/coexpression/atlas_variance_before_after_combat.csv")
  }
} else {
  log_warn("only %d common genes -- skipping ComBat integration", length(common))
}

log_session_info(SCRIPT)
log_step("COMPLETE: healthy EC co-expression atlas")
