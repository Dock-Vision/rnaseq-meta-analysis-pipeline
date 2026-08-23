#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 01_fetch_dataset.R
#
# Fetch ONE GEO series and report what it actually contains, before committing
# it to the pipeline. Run this on every candidate dataset first: it is much
# cheaper to find out here that a series ships FPKM rather than raw counts, or
# that its 120 "samples" are 60 donors sequenced twice, than to discover it
# halfway through a meta-analysis.
#
# Steps:
#   1. Fetch the GEO series record.
#   2. Extract sample metadata          -> 02_metadata/<GSE>_metadata.csv
#   3. List GSE-level supplementary files and download anything that looks
#      like a count matrix              -> 03_raw_data/case_control/<GSE>/
#   4. If GEO ships only per-sample files, download all of them and merge into
#      one gene x sample table          -> 04_processed/<GSE>_counts.csv
#   5. Report dimensions, a preview, a technical-replicate warning, and
#      whether the values are raw counts or a normalised measure
#      (TPM/FPKM/CPM) -- which decides whether DESeq2 is valid at all.
#
# Usage:
#   Rscript 01_scripts/01_fetch_dataset.R GSE112087
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(data.table)
})

# ---- configuration --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !grepl("^GSE[0-9]+$", args[1])) {
  stop("Usage: Rscript 01_scripts/01_fetch_dataset.R <GSE accession>\n",
       "  e.g. Rscript 01_scripts/01_fetch_dataset.R GSE112087")
}
GSE_ID <- args[1]

PROJ_ROOT <- normalizePath(
  file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                             value = TRUE)[1])), ".."),
  mustWork = FALSE
)
if (is.na(PROJ_ROOT) || !dir.exists(PROJ_ROOT)) PROJ_ROOT <- getwd()

META_DIR <- file.path(PROJ_ROOT, "02_metadata")
RAW_DIR  <- file.path(PROJ_ROOT, "03_raw_data", "case_control", GSE_ID)
PROC_DIR <- file.path(PROJ_ROOT, "04_processed")

for (d in c(META_DIR, RAW_DIR, PROC_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# GEO files are large and the NCBI FTP endpoint is often slow.
options(timeout = max(1800, getOption("timeout")))
options(download.file.method.GEOquery = "auto")

msg <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function(title = NULL) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  if (!is.null(title)) cat(title, "\n", strrep("=", 78), "\n", sep = "")
}

# ---- helpers --------------------------------------------------------------

# GEO characteristics fields arrive as "key: value" strings, one column per
# slot, and the key can differ between samples within the same slot. Parse each
# cell into key/value and re-assemble as one column per distinct key so the
# metadata table is actually usable downstream.
parse_characteristics <- function(pdata) {
  ch_cols <- grep("^characteristics_ch", colnames(pdata), value = TRUE)
  if (length(ch_cols) == 0) return(NULL)

  cells <- lapply(ch_cols, function(cl) as.character(pdata[[cl]]))
  n <- nrow(pdata)

  keys <- character(0)
  parsed <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    v <- cells[[i]]
    has_sep <- grepl(":", v, fixed = TRUE)
    k <- ifelse(has_sep, trimws(sub(":.*$", "", v)), ch_cols[i])
    val <- ifelse(has_sep, trimws(sub("^[^:]*:\\s*", "", v)), trimws(v))
    parsed[[i]] <- data.frame(key = k, value = val,
                              idx = seq_len(n), stringsAsFactors = FALSE)
    keys <- union(keys, unique(k[!is.na(k) & nzchar(k)]))
  }
  parsed <- do.call(rbind, parsed)
  parsed <- parsed[!is.na(parsed$key) & nzchar(parsed$key), , drop = FALSE]

  out <- lapply(keys, function(k) {
    sub <- parsed[parsed$key == k, , drop = FALSE]
    # collapse in the rare case a sample repeats the same key
    tapply(sub$value, factor(sub$idx, levels = seq_len(n)),
           function(x) paste(unique(x), collapse = "; "))
  })
  names(out) <- make.names(paste0("char_", keys), unique = TRUE)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# Read a count-like table from txt/tsv/csv, plain or gzipped. First column is
# assumed to be the gene identifier.
read_count_table <- function(path) {
  dt <- tryCatch(
    fread(path, header = TRUE, data.table = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || ncol(dt) < 2 || nrow(dt) == 0) return(NULL)
  gene_ids <- as.character(dt[[1]])
  mat <- dt[, -1, drop = FALSE]
  # keep only numeric columns (drops annotation columns like length/symbol)
  num <- vapply(mat, function(x) is.numeric(x) || is.integer(x), logical(1))
  if (!any(num)) return(NULL)
  mat <- mat[, num, drop = FALSE]
  rownames(mat) <- make.unique(gene_ids)
  mat
}

# Does a filename look like a combined count matrix (as opposed to a bigwig,
# a per-sample file, or a processed DE result)?
looks_like_matrix <- function(fnames) {
  grepl("count|matrix|expression|expr|rawdata|fpkm|tpm|rpkm|genes",
        fnames, ignore.case = TRUE) &
    grepl("\\.(txt|tsv|csv|tab)(\\.gz|\\.bz2|\\.zip)?$", fnames,
          ignore.case = TRUE)
}

# Classify a numeric matrix as raw counts vs a normalized measure.
classify_values <- function(mat) {
  v <- as.matrix(mat)
  storage.mode(v) <- "double"
  finite <- v[is.finite(v)]
  if (length(finite) == 0) {
    return(list(verdict = "UNKNOWN", reasons = "no finite values"))
  }

  frac_int  <- mean(abs(finite - round(finite)) < 1e-8)
  any_neg   <- any(finite < 0)
  max_val   <- max(finite)
  col_sums  <- colSums(v, na.rm = TRUE)
  med_sum   <- stats::median(col_sums)
  sum_cv    <- if (med_sum > 0) stats::sd(col_sums) / mean(col_sums) else NA_real_

  reasons <- c(
    sprintf("integer-valued entries: %.1f%%", 100 * frac_int),
    sprintf("value range: %.3f to %.3f", min(finite), max_val),
    sprintf("median column sum: %s", format(med_sum, big.mark = ",",
                                            scientific = FALSE, digits = 6)),
    sprintf("column-sum CV: %s",
            if (is.na(sum_cv)) "NA" else sprintf("%.3f", sum_cv)),
    sprintf("negative values present: %s", any_neg)
  )

  # TPM/CPM columns sum to a fixed constant (1e6), so the CV across samples is
  # essentially zero. Raw count libraries vary substantially in depth.
  near_1e6 <- !is.na(sum_cv) && sum_cv < 0.01 &&
    med_sum > 9e5 && med_sum < 1.1e6

  verdict <- if (any_neg) {
    "NORMALIZED / TRANSFORMED (negative values -> log-scale or batch-corrected)"
  } else if (near_1e6) {
    "NORMALIZED (columns sum to ~1e6 -> TPM or CPM)"
  } else if (frac_int > 0.999) {
    "RAW COUNTS (all-integer, variable library sizes)"
  } else if (frac_int > 0.90) {
    "LIKELY RAW COUNTS (mostly integer; possibly estimated counts, e.g. RSEM/salmon)"
  } else {
    "NORMALIZED (non-integer values -> TPM/FPKM/RPKM or similar)"
  }

  list(verdict = verdict, reasons = reasons,
       frac_int = frac_int, col_sums = col_sums)
}

# ---- 1/2. series record + sample metadata ---------------------------------

rule(sprintf("STEP 1-2: fetching %s and extracting sample metadata", GSE_ID))

gse <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) {
  msg("Series matrix objects returned: %d (%s)",
      length(gse), paste(names(gse), collapse = ", "))
  eset <- gse[[1]]
} else {
  eset <- gse
}

pdata <- pData(eset)
msg("Samples in series matrix: %d", nrow(pdata))

base_cols <- data.frame(
  gsm_id = as.character(pdata$geo_accession),
  title  = as.character(pdata$title),
  source = as.character(
    if (!is.null(pdata$source_name_ch1)) pdata$source_name_ch1 else NA
  ),
  organism = as.character(
    if (!is.null(pdata$organism_ch1)) pdata$organism_ch1 else NA
  ),
  platform_id = as.character(
    if (!is.null(pdata$platform_id)) pdata$platform_id else NA
  ),
  library_strategy = as.character(
    if (!is.null(pdata$`library strategy:ch1`)) pdata$`library strategy:ch1` else NA
  ),
  stringsAsFactors = FALSE
)

char_cols <- parse_characteristics(pdata)
metadata <- if (is.null(char_cols)) base_cols else cbind(base_cols, char_cols)

# also keep the raw characteristics_ch* strings verbatim for provenance
raw_ch <- grep("^characteristics_ch", colnames(pdata), value = TRUE)
if (length(raw_ch) > 0) {
  metadata <- cbind(
    metadata,
    setNames(lapply(raw_ch, function(cl) as.character(pdata[[cl]])),
             paste0("raw_", raw_ch))
  )
}

meta_path <- file.path(META_DIR, paste0(GSE_ID, "_metadata.csv"))
write.csv(metadata, meta_path, row.names = FALSE)
msg("Metadata written: %s  (%d rows x %d cols)",
    meta_path, nrow(metadata), ncol(metadata))
msg("Metadata columns: %s", paste(colnames(metadata), collapse = ", "))

cat("\n--- metadata preview (first 5 rows, first 6 cols) ---\n")
print(head(metadata[, seq_len(min(6, ncol(metadata))), drop = FALSE], 5))

# ---- 3. list + download GSE-level supplementary files ----------------------

rule("STEP 3: GSE-level supplementary files")

supp_list <- tryCatch(
  getGEOSuppFiles(GSE_ID, makeDirectory = FALSE, baseDir = RAW_DIR,
                  fetch_files = FALSE),
  error = function(e) { msg("Could not list supp files: %s", conditionMessage(e)); NULL }
)

supp_names <- character(0)
if (!is.null(supp_list) && nrow(supp_list) > 0) {
  supp_names <- if ("fname" %in% colnames(supp_list)) {
    as.character(supp_list$fname)
  } else {
    rownames(supp_list)
  }
  cat("\nSupplementary files listed for", GSE_ID, ":\n")
  for (i in seq_along(supp_names)) {
    sz <- if ("size" %in% colnames(supp_list)) supp_list$size[i] else NA
    msg("  [%d] %s%s", i, supp_names[i],
        if (is.na(sz)) "" else sprintf("  (%s bytes)", sz))
  }
} else {
  cat("\nNo GSE-level supplementary files reported.\n")
}

matrix_candidates <- supp_names[looks_like_matrix(supp_names)]
counts <- NULL
source_desc <- NA_character_

if (length(matrix_candidates) > 0) {
  msg("\n%d file(s) look like a combined count/expression matrix.",
      length(matrix_candidates))

  # Prefer a file whose name mentions counts over one mentioning TPM/FPKM.
  pref <- order(!grepl("count", matrix_candidates, ignore.case = TRUE),
                grepl("tpm|fpkm|rpkm|norm", matrix_candidates, ignore.case = TRUE))
  matrix_candidates <- matrix_candidates[pref]

  # escape regex metacharacters in the filenames (']' first, '\\' last inside
  # the class so TRE parses it as a literal set)
  esc <- gsub("([][{}().^$*+?|\\\\])", "\\\\\\1", matrix_candidates)
  invisible(getGEOSuppFiles(
    GSE_ID, makeDirectory = FALSE, baseDir = RAW_DIR,
    filter_regex = paste0("^(", paste(esc, collapse = "|"), ")$"),
    fetch_files = TRUE
  ))

  for (fn in matrix_candidates) {
    fp <- file.path(RAW_DIR, fn)
    if (!file.exists(fp)) next
    msg("Downloaded: %s (%s bytes)", fp,
        format(file.size(fp), big.mark = ","))
    cand <- read_count_table(fp)
    if (!is.null(cand) && ncol(cand) > 1) {
      counts <- cand
      source_desc <- sprintf("GSE-level combined matrix: %s", fn)
      msg("Parsed as matrix: %d rows x %d numeric columns",
          nrow(cand), ncol(cand))
      break
    } else {
      msg("Could not parse %s as a gene x sample matrix; trying next.", fn)
    }
  }
}

# ---- 4. fall back to per-sample files and merge ----------------------------

if (is.null(counts)) {
  rule("STEP 4: no combined matrix usable - downloading per-sample files")

  gsm_ids <- metadata$gsm_id
  gsm_dir <- file.path(RAW_DIR, "per_sample")
  if (!dir.exists(gsm_dir)) dir.create(gsm_dir, recursive = TRUE)

  per_sample <- list()
  for (i in seq_along(gsm_ids)) {
    g <- gsm_ids[i]
    msg("[%d/%d] %s", i, length(gsm_ids), g)
    files <- tryCatch(
      getGEOSuppFiles(g, makeDirectory = TRUE, baseDir = gsm_dir,
                      fetch_files = TRUE),
      error = function(e) { msg("   download failed: %s", conditionMessage(e)); NULL }
    )
    if (is.null(files) || nrow(files) == 0) next

    fps <- rownames(files)
    fps <- fps[grepl("\\.(txt|tsv|csv|tab|counts)(\\.gz)?$", fps, ignore.case = TRUE)]
    for (fp in fps) {
      tb <- read_count_table(fp)
      if (is.null(tb)) next
      # per-sample file: take the first numeric column as this sample's counts
      per_sample[[g]] <- setNames(
        data.frame(gene_id = rownames(tb), value = tb[[1]],
                   stringsAsFactors = FALSE),
        c("gene_id", g)
      )
      break
    }
  }

  if (length(per_sample) == 0) {
    stop("No per-sample count files could be downloaded or parsed for ", GSE_ID,
         ". Inspect ", RAW_DIR, " manually.")
  }

  msg("Merging %d per-sample tables...", length(per_sample))
  merged <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE),
                   per_sample)
  gene_ids <- merged$gene_id
  counts <- merged[, -1, drop = FALSE]
  rownames(counts) <- make.unique(as.character(gene_ids))
  source_desc <- sprintf("merged from %d per-sample GEO files", length(per_sample))
}

# ---- map column names back to GSM IDs where possible ----------------------

if (!is.null(counts)) {
  cn <- colnames(counts)
  if (!any(cn %in% metadata$gsm_id)) {
    # Count-matrix headers are usually the run/library name. GEO titles may be
    # the run name, or "<donor> <run>", so try progressively looser keys.
    keys <- list(
      title      = metadata$title,
      last_token = sub("^.*\\s", "", metadata$title),
      desc       = if (!is.null(pdata$description)) as.character(pdata$description) else NULL
    )
    for (k in names(keys)) {
      if (is.null(keys[[k]])) next
      hit <- match(cn, keys[[k]])
      if (sum(!is.na(hit)) > 0.5 * length(cn)) {
        msg("\nRemapped %d/%d column names to GSM IDs via '%s'.",
            sum(!is.na(hit)), length(cn), k)
        # keep the original run label alongside the GSM ID in the metadata
        metadata$count_column <- NA_character_
        metadata$count_column[hit[!is.na(hit)]] <- cn[!is.na(hit)]
        colnames(counts) <- ifelse(is.na(hit), cn, metadata$gsm_id[hit])
        write.csv(metadata, meta_path, row.names = FALSE)
        msg("Metadata updated with 'count_column' link: %s", meta_path)
        break
      }
    }
  }
}

counts_path <- file.path(PROC_DIR, paste0(GSE_ID, "_counts.csv"))
out <- cbind(gene_id = rownames(counts), counts)
fwrite(out, counts_path)
msg("\nCount table written: %s", counts_path)

# ---- 5. report -------------------------------------------------------------

rule(sprintf("STEP 5: SUMMARY - %s", GSE_ID))

msg("Source            : %s", source_desc)
msg("Number of genes   : %s", format(nrow(counts), big.mark = ","))
msg("Number of samples : %d", ncol(counts))
msg("Metadata samples  : %d", nrow(metadata))

overlap <- sum(colnames(counts) %in% metadata$gsm_id)
msg("Columns matching a GSM ID in metadata: %d / %d", overlap, ncol(counts))
if (overlap < ncol(counts)) {
  msg("NOTE: column names are not all GSM IDs. First few: %s",
      paste(head(colnames(counts), 5), collapse = ", "))
}

# GEO sometimes registers one GSM per sequencing lane. Flag it here so the
# downstream script knows whether lanes need collapsing per biological sample.
if (!is.null(metadata$char_donor)) {
  n_donor <- length(unique(metadata$char_donor))
  if (n_donor < nrow(metadata)) {
    msg("\nNOTE: %d GSMs map to %d unique donors (~%.1f runs/donor).",
        nrow(metadata), n_donor, nrow(metadata) / n_donor)
    msg("      Technical replicates likely need collapsing before DE analysis.")
  }
}

cat("\n--- first 5 rows (first 8 sample columns) ---\n")
print(head(counts[, seq_len(min(8, ncol(counts))), drop = FALSE], 5))

cls <- classify_values(counts)
cat("\n--- value scale diagnostics ---\n")
for (r in cls$reasons) msg("  %s", r)
cat("\n")
msg("VERDICT: values look like -> %s", cls$verdict)

if (grepl("^RAW|^LIKELY RAW", cls$verdict)) {
  cat("\nLibrary sizes (millions of reads), first 8 samples:\n")
  print(round(head(cls$col_sums, 8) / 1e6, 2))
  cat("\nThis matrix is suitable for DESeq2 input as-is.\n")
} else {
  cat("\nThis matrix is NOT suitable for DESeq2 (which requires raw counts).\n")
  cat("Use limma-voom / limma-trend on log2 values, or locate a raw-count file.\n")
}

rule("DONE")
