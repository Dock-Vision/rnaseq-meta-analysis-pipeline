# ---------------------------------------------------------------------------
# R/prepare.R -- from a raw GEO matrix to an analysis-ready DESeqDataSet.
#
# Handles the things that silently corrupt a meta-analysis if missed:
#   * lane-level GSMs (technical replicates) inflating n
#   * non-integer / normalized matrices reaching DESeq2
#   * low-count gene filtering tied to the smallest group size
#   * PCA + variance decomposition before and after batch correction
#     (METHODS S3)
# ---------------------------------------------------------------------------

# ---- technical replicate detection ---------------------------------------

# Detects whether several GSMs describe the same biological sample. Returns the
# metadata column to collapse on, or NULL if every GSM is already unique.

# Metadata fields that legitimately differ between technical replicates of the
# same library (per-run QC numbers), so variation here does NOT mean the
# samples are biologically distinct.
.TECHNICAL_FIELDS <- c("lib.size", "libsize", "norm.factors", "rin",
                       "rnaconcentration", "concentration", "volum", "volume",
                       "processingdate", "date", "batch", "lane", "run",
                       "flowcell", "barcode", "index", "readcount", "reads")

# Classifies what a repeated key actually represents. Two GSMs sharing a donor
# ID can be either:
#   * technical replicates  -- same library sequenced twice (GSE112087 lanes)
#   * biological sub-samples -- distinct specimens from one donor
#     (GSE213001: apex vs base of the same lung)
# Both must be aggregated to avoid pseudoreplication inflating n, but they are
# scientifically different and must be labelled correctly.
classify_replicates <- function(meta, key) {
  grp <- as.character(meta[[key]])
  other <- setdiff(grep("^char_", colnames(meta), value = TRUE), key)
  other <- other[!grepl(paste(.TECHNICAL_FIELDS, collapse = "|"),
                        tolower(other))]

  varying <- character(0)
  for (cl in other) {
    v <- as.character(meta[[cl]])
    n_within <- tapply(v, grp, function(x) length(unique(x)))
    if (any(n_within > 1, na.rm = TRUE)) varying <- c(varying, cl)
  }
  list(type = if (length(varying) == 0) "technical" else "biological_subsample",
       varying_fields = varying)
}

detect_replicate_key <- function(meta, candidates = NULL) {
  if (is.null(candidates)) {
    candidates <- grep("^char_(donor|subject|patient|individual|sample|id)",
                       colnames(meta), value = TRUE, ignore.case = TRUE)
    # a repeated title is also a reliable signal (e.g. GSE112087 lanes)
    candidates <- c(candidates, "title", "first_token_title")
  }
  meta$first_token_title <- sub("\\s.*$", "", meta$title)

  n <- nrow(meta)
  for (cand in candidates) {
    if (is.null(meta[[cand]])) next
    k <- length(unique(meta[[cand]]))
    if (k < n && k > 1) {
      cls <- classify_replicates(meta, cand)
      log_info("replicate key '%s': %d GSMs -> %d unique (%.1f per group) [%s]",
               cand, n, k, n / k, cls$type)
      if (cls$type == "biological_subsample") {
        log_warn("within-donor samples differ on: %s",
                 paste(cls$varying_fields, collapse = ", "))
        log_warn("these are distinct specimens, not technical replicates -- aggregating to donor level to avoid pseudoreplication, but the within-donor contrast is discarded")
      }
      return(structure(cand, type = cls$type,
                       varying_fields = cls$varying_fields))
    }
  }
  NULL
}

# Sums counts across technical replicates. Summing (not averaging) is correct
# for lane-level replicates: they are independent reads of the same library, so
# the pooled library is the sum.

# `first_token_title` is derived inside detect_replicate_key() on its own copy
# of the metadata, so any caller receiving that key back must derive it too.
.ensure_key_column <- function(meta, key) {
  key <- as.character(key)
  if (!is.null(meta[[key]])) return(meta)
  if (key == "first_token_title") {
    meta$first_token_title <- sub("\\s.*$", "", meta$title)
    return(meta)
  }
  stop("replicate key '", key, "' is not a column of the metadata")
}

collapse_technical_replicates <- function(counts, meta, key) {
  meta <- .ensure_key_column(meta, key)
  key <- as.character(key)
  # Align by name rather than assuming identical order -- remapping column
  # names to GSM IDs does not guarantee they come back in metadata order.
  if (!identical(colnames(counts), meta$gsm_id)) {
    common <- intersect(meta$gsm_id, colnames(counts))
    if (length(common) < 2) {
      stop("counts and metadata share only ", length(common),
           " sample IDs -- cannot collapse")
    }
    if (length(common) < nrow(meta)) {
      log_warn("%d metadata samples have no count column and are dropped",
               nrow(meta) - length(common))
    }
    if (length(common) < ncol(counts)) {
      log_warn("%d count columns have no metadata row and are dropped",
               ncol(counts) - length(common))
    }
    meta <- meta[match(common, meta$gsm_id), , drop = FALSE]
    counts <- counts[, common, drop = FALSE]
  }
  grp <- as.character(meta[[key]])

  collapsed <- t(rowsum(t(as.matrix(counts)), group = grp, reorder = FALSE))
  # one metadata row per group; take the first GSM's annotation, and record
  # which GSMs were merged
  first_idx <- match(colnames(collapsed), grp)
  new_meta <- meta[first_idx, , drop = FALSE]
  new_meta$n_technical_reps <- as.integer(table(grp)[colnames(collapsed)])
  new_meta$merged_gsms <- vapply(colnames(collapsed), function(g)
    paste(meta$gsm_id[grp == g], collapse = ";"), character(1))
  new_meta$sample_id <- colnames(collapsed)
  rownames(new_meta) <- NULL

  log_info("collapsed %d -> %d samples on '%s' (counts summed)",
           ncol(counts), ncol(collapsed), key)
  list(counts = collapsed, metadata = new_meta)
}

# ---- group assignment -----------------------------------------------------

# Maps free-text phenotype labels onto Case / Control. Explicit and inspectable
# rather than guessed inside the DE call.

# Preference order for the phenotype column. A first-match grep is not enough:
# GSE213001 has a `char_group` field holding "IPF.Apex"/"NDC.Base", i.e. disease
# crossed with anatomical site, which is not a case/control variable at all.
.GROUP_COL_PREF <- c("diseasenormal", "disease.?status", "diseasestate",
                     "phenotype", "diagnosis", "diseasegroup", "disease",
                     "condition", "status", "subject.?group", "case",
                     "treatment", "group")

# Picks the phenotype column that actually yields a usable two-level contrast.
# Columns that vary within a donor cannot be sample-level phenotype and are
# excluded.
pick_group_column <- function(meta, exclude = character(0), min_per_group = 2) {
  cands <- grep("^char_", colnames(meta), value = TRUE)
  cands <- setdiff(cands, exclude)
  if (length(cands) == 0) return(NULL)

  ranked <- character(0)
  for (pat in .GROUP_COL_PREF) {
    hit <- grep(pat, tolower(cands), value = TRUE)
    ranked <- c(ranked, setdiff(hit, ranked))
  }
  ranked <- c(ranked, setdiff(cands, ranked))

  for (cl in ranked) {
    g <- suppressWarnings(assign_groups(meta, cl, verbose = FALSE))
    tab <- table(g)
    if (length(tab) == 2 && all(tab >= min_per_group)) {
      log_info("phenotype column '%s' -> %s", cl,
               paste(sprintf("%s=%d", names(tab), tab), collapse = ", "))
      return(cl)
    }
  }
  log_warn("no column yields a two-level contrast with >=%d per group",
           min_per_group)
  NULL
}

assign_groups <- function(meta, group_col,
                          control_patterns = c("healthy", "control", "normal",
                                               "non-?diseased", "\\bndc\\b",
                                               "\\bhc\\b", "\\bctrl\\b",
                                               "untreated", "baseline",
                                               "unaffected", "\\bnc\\b",
                                               # negated-disease phrasings, e.g.
                                               # GSE202625's "no coronary artery
                                               # disease (CAD)"
                                               "^no[ _-]", "\\bwithout\\b",
                                               "disease[- ]free", "\\bnegative\\b"),
                          verbose = TRUE) {
  raw <- tolower(trimws(as.character(meta[[group_col]])))
  is_ctrl <- Reduce(`|`, lapply(control_patterns, function(p) grepl(p, raw)))
  grp <- ifelse(is_ctrl, "Control", "Case")

  tab <- table(grp)
  if (verbose) {
    log_info("group assignment from '%s': %s",
             group_col, paste(sprintf("%s=%d", names(tab), tab), collapse = ", "))
    if (length(unique(grp)) < 2) {
      log_warn("only one group present -- this dataset cannot support a DE contrast")
    }
  }
  factor(grp, levels = c("Control", "Case"))
}

# ---- batch covariate selection -------------------------------------------

# Picks a batch variable to include in the DE design, but ONLY when it is safe:
# a batch that is confounded with case/control cannot be adjusted for without
# removing the disease effect itself. GSE202625 is the motivating case --
# sequencing batch explains ~15% of variance versus ~1% for disease, yet is
# balanced across groups (chi-square p = 0.94), so modelling it recovers power
# that would otherwise be lost.
select_batch_covariate <- function(meta, group_col = "group",
                                   max_levels_frac = 0.5, min_per_level = 2,
                                   min_confound_p = 0.05) {
  cands <- grep("batch|run|lane|flowcell|site|centre|center|processingdate",
                grep("^char_", colnames(meta), value = TRUE),
                value = TRUE, ignore.case = TRUE)
  if (length(cands) == 0 || is.null(meta[[group_col]])) return(NULL)

  n <- nrow(meta)
  for (cl in cands) {
    v <- factor(as.character(meta[[cl]]))
    k <- nlevels(v)
    if (k < 2 || k > max(2, floor(n * max_levels_frac))) next
    if (any(table(v) < min_per_level)) next

    tb <- table(v, meta[[group_col]])
    # every batch level must contain both groups, or the design is not estimable
    if (any(rowSums(tb > 0) < 2)) {
      log_warn("batch '%s' has levels containing only one group -- not estimable, skipping",
               cl)
      next
    }
    p <- suppressWarnings(stats::chisq.test(tb)$p.value)
    if (is.na(p) || p < min_confound_p) {
      log_warn("batch '%s' is confounded with %s (chi-square p = %.3g) -- NOT adjusting; the disease effect cannot be separated from batch",
               cl, group_col, p)
      next
    }
    log_info("batch '%s' (%d levels) is balanced across %s (chi-square p = %.2f) -- including in the design",
             cl, k, group_col, p)
    return(cl)
  }
  NULL
}

# ---- DESeqDataSet construction -------------------------------------------

# Low-count filter: keep genes with >= min_count reads in at least as many
# samples as the smallest group. Standard, and avoids filtering that is biased
# toward the larger group.
filter_low_counts <- function(counts, group, min_count = MIN_COUNT) {
  min_samples <- min(table(group))
  keep <- rowSums(counts >= min_count) >= min_samples
  log_info("gene filter: %s -> %s genes (>=%d reads in >=%d samples)",
           fmt_n(nrow(counts)), fmt_n(sum(keep)), min_count, min_samples)
  counts[keep, , drop = FALSE]
}

build_dds <- function(counts, meta, design = ~ group, filter = TRUE) {
  require_pkg("DESeq2", "differential expression")

  mat <- as.matrix(counts)
  if (any(abs(mat - round(mat)) > 1e-8, na.rm = TRUE)) {
    log_warn("non-integer values present -- rounding for DESeq2 (estimated counts)")
  }
  storage.mode(mat) <- "integer"
  mat[is.na(mat)] <- 0L

  if (filter) mat <- filter_low_counts(mat, meta$group)

  # DESeq2 requires colData rownames to match countData colnames exactly.
  # Subsetting metadata during alignment leaves stale integer rownames, so
  # reset them here rather than relying on the caller.
  meta <- as.data.frame(meta)
  rownames(meta) <- NULL
  if (!is.null(meta$sample_id) && setequal(meta$sample_id, colnames(mat))) {
    meta <- meta[match(colnames(mat), meta$sample_id), , drop = FALSE]
    rownames(meta) <- meta$sample_id
  } else if (nrow(meta) == ncol(mat)) {
    rownames(meta) <- colnames(mat)
  } else {
    stop("metadata rows (", nrow(meta), ") != count columns (", ncol(mat), ")")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                        colData = meta,
                                        design = design)
  log_info("DESeqDataSet: %s genes x %d samples, design = %s",
           fmt_n(nrow(dds)), ncol(dds), paste(deparse(design), collapse = ""))
  dds
}

# ---- variance structure (METHODS S3) -------------------------------

# Percent of variance in each of the top PCs attributable to each covariate,
# via PC ~ covariate regression. This is the quantitative statement the Quality
# Bar requires alongside the before/after PCA plots.

# GEO returns every characteristic as a character string, so a genuinely
# continuous covariate like age arrives as text. Left alone, lm() treats it as
# a factor with one level per distinct age, burning a degree of freedom per
# level and inflating R^2 towards 1. Coerce anything that is fully numeric.
.coerce_covariate <- function(v, name) {
  if (is.numeric(v)) return(v)
  chr <- trimws(as.character(v))
  nonmiss <- chr[!is.na(chr) & nzchar(chr)]
  if (length(nonmiss) && !any(is.na(suppressWarnings(as.numeric(nonmiss))))) {
    log_info("covariate '%s': character but fully numeric -> treating as continuous",
             name)
    return(suppressWarnings(as.numeric(chr)))
  }
  factor(chr)
}

pc_variance_explained <- function(mat_vst, meta, covariates, n_pc = 5) {
  pca <- stats::prcomp(t(mat_vst), center = TRUE, scale. = FALSE)
  var_pct <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  n_pc <- min(n_pc, ncol(pca$x))
  n <- nrow(meta)

  out <- list()
  for (cv in covariates) {
    if (is.null(meta[[cv]])) next
    v <- .coerce_covariate(meta[[cv]], cv)
    if (length(unique(v[!is.na(v)])) < 2) next

    # A categorical covariate with many levels relative to n explains variance
    # by sheer degrees of freedom. Report adjusted R^2 in that case and warn.
    if (is.factor(v) && nlevels(droplevels(v)) > n / 3) {
      log_warn("covariate '%s' has %d levels on n=%d -- high-cardinality factor, using adjusted R^2",
               cv, nlevels(droplevels(v)), n)
    }
    for (i in seq_len(n_pc)) {
      fit <- stats::lm(pca$x[, i] ~ v)
      s <- summary(fit)
      f <- s$fstatistic
      p <- if (!is.null(f)) stats::pf(f[1], f[2], f[3], lower.tail = FALSE) else NA_real_
      # adjusted R^2 penalises the extra parameters of a multi-level factor;
      # floor at 0 since a negative adjusted R^2 means "no signal"
      r2_adj <- max(0, s$adj.r.squared)
      out[[length(out) + 1]] <- data.frame(
        covariate = cv, PC = paste0("PC", i),
        pc_variance_pct = round(var_pct[i], 2),
        df_used = if (is.factor(v)) nlevels(droplevels(v)) - 1L else 1L,
        r_squared = round(s$r.squared, 4),
        r_squared_adj = round(r2_adj, 4),
        # share of TOTAL dataset variance explained by this covariate via this
        # PC, using adjusted R^2 so covariates with different df are comparable
        total_variance_explained_pct = round(var_pct[i] * r2_adj, 3),
        p_value = signif(p, 3),
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) {
    res$p_adj <- signif(stats::p.adjust(res$p_value, method = "BH"), 3)
  }
  attr(res, "pca") <- pca
  attr(res, "var_pct") <- var_pct
  res
}

# Convenience: total variance attributable to a covariate across the top PCs.
summarise_variance <- function(pcv) {
  if (is.null(pcv)) return(NULL)
  agg <- stats::aggregate(total_variance_explained_pct ~ covariate,
                          data = pcv, FUN = sum)
  agg <- agg[order(-agg$total_variance_explained_pct), ]
  colnames(agg)[2] <- "total_var_explained_pct"
  agg$total_var_explained_pct <- round(agg$total_var_explained_pct, 2)
  agg
}

# ---- transformation -------------------------------------------------------

# VST for visualisation / PCA / correlation. blind = TRUE when inspecting for
# batch structure so the design does not shape the picture.
vst_matrix <- function(dds, blind = TRUE) {
  require_pkg("DESeq2", "variance stabilizing transform")
  n <- ncol(dds)
  vsd <- if (n >= 30) {
    DESeq2::vst(dds, blind = blind)
  } else {
    log_info("n=%d small -- using varianceStabilizingTransformation", n)
    DESeq2::varianceStabilizingTransformation(dds, blind = blind)
  }
  SummarizedExperiment::assay(vsd)
}
