# ---------------------------------------------------------------------------
# R/download.R -- generalized GEO fetch.
#
# Generalises the inspection logic in 01_fetch_dataset.R so any series in the
# registry can be pulled with one call:
#
#     res <- fetch_gse("GSE213001")
#     res$metadata   # sample table
#     res$counts     # gene x sample matrix
#
# Re-runnable: cached files in 03_raw_data/ are reused rather than re-downloaded.
# ---------------------------------------------------------------------------

# ---- metadata -------------------------------------------------------------

# GEO characteristics arrive as "key: value" strings, one column per slot, and
# the key can differ between samples in the same slot. Parse into one tidy
# column per distinct key.
parse_characteristics <- function(pdata) {
  ch_cols <- grep("^characteristics_ch", colnames(pdata), value = TRUE)
  if (length(ch_cols) == 0) return(NULL)

  n <- nrow(pdata)
  keys <- character(0)
  parsed <- list()
  for (i in seq_along(ch_cols)) {
    v <- as.character(pdata[[ch_cols[i]]])
    has_sep <- grepl(":", v, fixed = TRUE)
    k   <- ifelse(has_sep, trimws(sub(":.*$", "", v)), ch_cols[i])
    val <- ifelse(has_sep, trimws(sub("^[^:]*:\\s*", "", v)), trimws(v))
    parsed[[i]] <- data.frame(key = k, value = val, idx = seq_len(n),
                              stringsAsFactors = FALSE)
    keys <- union(keys, unique(k[!is.na(k) & nzchar(k)]))
  }
  parsed <- do.call(rbind, parsed)
  parsed <- parsed[!is.na(parsed$key) & nzchar(parsed$key), , drop = FALSE]

  out <- lapply(keys, function(k) {
    sub <- parsed[parsed$key == k, , drop = FALSE]
    tapply(sub$value, factor(sub$idx, levels = seq_len(n)),
           function(x) paste(unique(x), collapse = "; "))
  })
  names(out) <- make.names(paste0("char_", keys), unique = TRUE)
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

get_gse_metadata <- function(gse, save = TRUE) {
  log_step("fetching series record: %s", gse)

  obj <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
  eset <- if (is.list(obj)) obj[[1]] else obj
  pdata <- Biobase::pData(eset)
  log_info("%d samples in series matrix", nrow(pdata))

  pick <- function(nm) {
    if (!is.null(pdata[[nm]])) as.character(pdata[[nm]]) else NA_character_
  }
  base <- data.frame(
    gsm_id      = as.character(pdata$geo_accession),
    title       = pick("title"),
    source      = pick("source_name_ch1"),
    organism    = pick("organism_ch1"),
    platform_id = pick("platform_id"),
    stringsAsFactors = FALSE
  )

  ch <- parse_characteristics(pdata)
  meta <- if (is.null(ch)) base else cbind(base, ch)

  # keep raw characteristics strings for provenance / manual auditing
  raw_ch <- grep("^characteristics_ch", colnames(pdata), value = TRUE)
  if (length(raw_ch)) {
    meta <- cbind(meta, setNames(
      lapply(raw_ch, function(cl) as.character(pdata[[cl]])),
      paste0("raw_", raw_ch)))
  }

  if (save) {
    path <- file.path(PATHS$metadata, paste0(gse, "_metadata.csv"))
    utils::write.csv(meta, path, row.names = FALSE)
    log_info("metadata -> %s (%d x %d)", path, nrow(meta), ncol(meta))
  }
  attr(meta, "pdata") <- pdata
  meta
}

# ---- supplementary files --------------------------------------------------

.looks_like_matrix <- function(f) {
  # include quantifier tool names -- several series name the combined matrix
  # after the tool that produced it (GSE131681_kallisto.GRCh37.xlsx) rather
  # than using the word "counts", which would otherwise force the much slower
  # and lossier per-sample download path.
  grepl(paste0("count|matrix|expression|expr|rawdata|fpkm|tpm|rpkm|genes|",
               "kallisto|salmon|rsem|htseq|featurecount|star|abundance|",
               "quant|est_count|gene_level|readspergene"), f,
        ignore.case = TRUE) &
    grepl("\\.(txt|tsv|csv|tab|xlsx|xls)(\\.gz|\\.bz2|\\.zip)?$", f, ignore.case = TRUE)
}

# Column names that are gene annotation, not samples. Some GEO matrices carry a
# full annotation block (GSE72509 ships CHR/START/END/SIZE/SYMBOL/DESC/...), and
# the numeric ones among them (START, END, SIZE) would otherwise be mistaken for
# extra samples and silently enter the analysis.
.ANNOT_COLS <- c("chr", "chrom", "chromosome", "start", "end", "str", "strand",
                 "coords", "size", "length", "width", "type", "biotype",
                 "gene_biotype", "symbol", "gene_symbol", "genesymbol", "desc",
                 "description", "genename", "gene_name", "feature_id",
                 "featureid", "gene", "gene_id", "geneid", "entrez",
                 "entrezid", "ensembl", "ensembl_id", "aliases", "alias",
                 "transcript", "transcript_id", "tx_id", "exonlength",
                 "effective_length", "gc", "gc_content")

# Preference order for which column supplies the gene identifier.
.ID_COL_PREF <- c("symbol", "gene_symbol", "genesymbol", "gene_name", "genename",
                  "gene_id", "geneid", "ensembl", "ensembl_id", "gene",
                  "feature_id", "featureid")

.read_count_table <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) {
    log_warn("empty or missing file: %s", basename(path))
    return(NULL)
  }

  # Salmon/kallisto per-sample output (quant.sf / abundance.tsv) is
  # TRANSCRIPT-level: Name, Length, EffectiveLength, TPM, NumReads. Taking it
  # as-is would put transcript IDs where gene IDs belong, so pick the read-count
  # column explicitly and mark the table for transcript->gene collapsing.
  # match any .sf file, not just files literally named quant.sf -- GSE294621
  # ships GSM8912599_D4M60C2.sf.gz
  if (grepl("\\.sf(\\.gz)?$|quant\\.sf|abundance\\.(tsv|h5)", basename(path),
            ignore.case = TRUE)) {
    tb <- tryCatch(data.table::fread(path, header = TRUE, data.table = FALSE,
                                     check.names = FALSE),
                   error = function(e) NULL)
    if (is.null(tb) || ncol(tb) < 2) return(NULL)
    val_col <- intersect(c("NumReads", "est_counts", "TPM", "tpm"), colnames(tb))
    if (length(val_col) == 0) return(NULL)
    mat <- data.frame(count = tb[[val_col[1]]])
    rownames(mat) <- make.unique(as.character(tb[[1]]))
    attr(mat, "transcript_level") <- TRUE
    attr(mat, "value_col") <- val_col[1]
    return(mat)
  }

  # STAR's ReadsPerGene.out.tab: no header, 4 columns
  # (gene, unstranded, forward, reverse) with 4 leading N_* summary rows.
  if (grepl("ReadsPerGene", basename(path), ignore.case = TRUE)) {
    tb <- tryCatch(data.table::fread(path, header = FALSE, data.table = FALSE),
                   error = function(e) NULL)
    if (is.null(tb) || ncol(tb) < 2) return(NULL)
    tb <- tb[!grepl("^N_", tb[[1]]), , drop = FALSE]   # drop N_unmapped etc.
    mat <- data.frame(count = tb[[2]])                  # unstranded column
    rownames(mat) <- make.unique(as.character(tb[[1]]))
    return(mat)
  }

  dt <- if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      log_warn("readxl not installed -- cannot read %s", basename(path))
      return(NULL)
    }
    tryCatch(as.data.frame(readxl::read_excel(path, sheet = 1),
                           check.names = FALSE),
             error = function(e) NULL)
  } else {
    # Some matrices carry comment/annotation lines above the real header
    # (GSE236730 starts with a '#group' row), which would otherwise be parsed
    # as the header and leave no numeric columns at all.
    n_skip <- 0L
    con <- tryCatch(gzfile(path, "rt"), error = function(e) NULL)
    if (!is.null(con)) {
      probe <- tryCatch(readLines(con, n = 10), error = function(e) character(0))
      close(con)
      while (n_skip < length(probe) && startsWith(probe[n_skip + 1], "#")) {
        n_skip <- n_skip + 1L
      }
    }
    if (n_skip > 0) log_info("skipping %d leading comment line(s)", n_skip)
    tryCatch(data.table::fread(path, header = TRUE, data.table = FALSE,
                               check.names = FALSE, skip = n_skip),
             error = function(e) NULL)
  }
  if (is.null(dt) || ncol(dt) < 2 || nrow(dt) == 0) return(NULL)

  cn_lower <- tolower(trimws(colnames(dt)))
  is_annot <- cn_lower %in% .ANNOT_COLS
  is_annot[1] <- TRUE   # first column is always an identifier, whatever its name

  # gene IDs: prefer a real symbol/ID column over a messy first column
  id_col <- 1L
  for (pref in .ID_COL_PREF) {
    hit <- which(cn_lower == pref)
    if (length(hit)) { id_col <- hit[1]; break }
  }
  ids <- as.character(dt[[id_col]])
  if (id_col != 1L) {
    log_info("using column '%s' as gene identifier (not column 1, '%s')",
             colnames(dt)[id_col], colnames(dt)[1])
  }

  # samples: numeric columns that are not part of the annotation block
  num <- vapply(dt, is.numeric, logical(1))
  keep <- num & !is_annot
  if (!any(keep)) return(NULL)
  dropped <- colnames(dt)[num & is_annot]
  if (length(dropped)) {
    log_info("dropped %d numeric annotation column(s): %s",
             length(dropped), paste(dropped, collapse = ", "))
  }

  mat <- dt[, keep, drop = FALSE]
  # blank/unmapped identifiers cannot be used downstream
  bad <- is.na(ids) | !nzchar(trimws(ids))
  if (any(bad)) {
    log_info("dropping %s rows with no gene identifier", fmt_n(sum(bad)))
    mat <- mat[!bad, , drop = FALSE]
    ids <- ids[!bad]
  }
  rownames(mat) <- make.unique(trimws(ids))
  mat
}

# Escape regex metacharacters so a literal filename can be used as a filter.
.esc <- function(x) gsub("([][{}().^$*+?|\\\\])", "\\\\\\1", x)

get_gse_counts <- function(gse, raw_dir = raw_dir_for(gse)) {
  log_step("supplementary files: %s", gse)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  supp <- tryCatch(
    GEOquery::getGEOSuppFiles(gse, makeDirectory = FALSE, baseDir = raw_dir,
                              fetch_files = FALSE),
    error = function(e) { log_warn("listing failed: %s", conditionMessage(e)); NULL })

  fnames <- character(0)
  if (!is.null(supp) && nrow(supp) > 0) {
    fnames <- if ("fname" %in% colnames(supp)) as.character(supp$fname) else rownames(supp)
    for (i in seq_along(fnames)) log_info("  [%d] %s", i, fnames[i])
  } else {
    log_warn("no GSE-level supplementary files listed")
  }

  cands <- fnames[.looks_like_matrix(fnames)]
  if (length(cands) == 0) return(list(counts = NULL, source = NA_character_,
                                      supp_files = fnames))

  # prefer a file named "counts" over one named TPM/FPKM
  cands <- cands[order(!grepl("count", cands, ignore.case = TRUE),
                       grepl("tpm|fpkm|rpkm|norm", cands, ignore.case = TRUE))]

  invisible(GEOquery::getGEOSuppFiles(
    gse, makeDirectory = FALSE, baseDir = raw_dir,
    filter_regex = paste0("^(", paste(.esc(cands), collapse = "|"), ")$"),
    fetch_files = TRUE))

  # Some series split the matrix across several files (GSE243193 ships
  # ..._GEO1.csv.gz and ..._GEO2.csv.gz, each holding a subset of samples).
  # Parse them all and merge on gene ID when their columns are disjoint.
  parsed <- list()
  for (fn in cands) {
    fp <- file.path(raw_dir, fn)
    if (!file.exists(fp)) next
    log_info("downloaded %s (%s bytes)", fn, fmt_n(file.size(fp)))
    m <- .read_count_table(fp)
    if (!is.null(m) && ncol(m) >= 1) {
      log_info("parsed %s: %s genes x %d columns", fn, fmt_n(nrow(m)), ncol(m))
      parsed[[fn]] <- m
    } else {
      log_warn("could not parse %s as a matrix", fn)
    }
  }
  if (length(parsed) == 0) {
    return(list(counts = NULL, source = NA_character_, supp_files = fnames))
  }

  if (length(parsed) == 1) {
    m <- parsed[[1]]
    if (ncol(m) < 2) {
      return(list(counts = NULL, source = NA_character_, supp_files = fnames))
    }
    return(list(counts = m,
                source = sprintf("GSE-level combined matrix: %s", names(parsed)[1]),
                supp_files = fnames))
  }

  overlap <- length(Reduce(intersect, lapply(parsed, colnames)))
  if (overlap == 0) {
    genes <- Reduce(intersect, lapply(parsed, rownames))
    log_info("merging %d matrix files on %s shared genes (disjoint samples)",
             length(parsed), fmt_n(length(genes)))
    # unname the list before cbind, or every column gets the filename prefixed
    mats <- unname(lapply(parsed, function(m) as.matrix(m[genes, , drop = FALSE])))
    merged <- as.data.frame(do.call(cbind, mats), check.names = FALSE)
    rownames(merged) <- genes
    return(list(counts = merged,
                source = sprintf("merged %d GSE-level files: %s",
                                 length(parsed), paste(names(parsed), collapse = ", ")),
                supp_files = fnames))
  }

  # overlapping columns: they are alternative renderings, take the largest
  best <- names(parsed)[which.max(vapply(parsed, ncol, integer(1)))]
  log_info("matrix files share columns; using the widest: %s", best)
  list(counts = parsed[[best]],
       source = sprintf("GSE-level combined matrix: %s", best),
       supp_files = fnames)
}

# Fallback: per-sample supplementary files, merged into one matrix.
get_per_sample_counts <- function(gse, gsm_ids, raw_dir = raw_dir_for(gse)) {
  log_step("per-sample download fallback: %s (%d GSMs)", gse, length(gsm_ids))
  gdir <- file.path(raw_dir, "per_sample")
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

  per <- list()
  for (i in seq_along(gsm_ids)) {
    g <- gsm_ids[i]
    if (i %% 10 == 1) log_info("  [%d/%d] %s", i, length(gsm_ids), g)
    files <- tryCatch(
      GEOquery::getGEOSuppFiles(g, makeDirectory = TRUE, baseDir = gdir,
                                fetch_files = TRUE),
      error = function(e) NULL)
    if (is.null(files) || nrow(files) == 0) next
    fps <- rownames(files)
    # include salmon (.sf) and kallisto (.tsv/.h5) per-sample quantification
    fps <- fps[grepl("\\.(txt|tsv|csv|tab|counts|sf|out\\.tab)(\\.gz)?$", fps,
                     ignore.case = TRUE)]
    fps <- fps[file.exists(fps) & file.size(fps) > 0]
    for (fp in fps) {
      tb <- .read_count_table(fp)
      if (is.null(tb)) next
      per[[g]] <- setNames(data.frame(gene_id = rownames(tb), v = tb[[1]],
                                      stringsAsFactors = FALSE),
                           c("gene_id", g))
      break
    }
  }
  if (length(per) == 0) return(list(counts = NULL, source = NA_character_))

  log_info("merging %d per-sample tables", length(per))
  merged <- Reduce(function(a, b) merge(a, b, by = "gene_id", all = TRUE), per)
  cnt <- merged[, -1, drop = FALSE]
  rownames(cnt) <- make.unique(as.character(merged$gene_id))
  # salmon/kallisto per-sample files are transcript-level -> sum to gene level
  cnt <- collapse_transcripts_to_genes(cnt)
  list(counts = cnt,
       source = sprintf("merged from %d per-sample GEO files", length(per)))
}

# ---- column -> GSM remapping ---------------------------------------------

# Count-matrix headers are usually run/library names, not GSM IDs. Try
# progressively looser keys and record the original label in the metadata.
remap_columns_to_gsm <- function(counts, meta, gse = NULL) {
  cn <- colnames(counts)
  if (any(cn %in% meta$gsm_id)) {
    meta$count_column <- meta$gsm_id
    return(list(counts = counts, metadata = meta, method = "already GSM"))
  }

  # explicit map from 00_config.R, for series no string rule can resolve
  if (!is.null(gse) && !is.null(MANUAL_COLUMN_MAP[[gse]])) {
    mp <- MANUAL_COLUMN_MAP[[gse]]
    hit <- mp[cn]
    if (sum(!is.na(hit)) > 0.5 * length(cn)) {
      log_info("remapped %d/%d columns to GSM IDs via MANUAL_COLUMN_MAP",
               sum(!is.na(hit)), length(cn))
      idx <- match(hit, meta$gsm_id)
      meta$count_column <- NA_character_
      meta$count_column[idx[!is.na(idx)]] <- cn[!is.na(idx)]
      colnames(counts) <- ifelse(is.na(hit), cn, hit)
      return(list(counts = counts, metadata = meta, method = "manual map"))
    }
    log_warn("MANUAL_COLUMN_MAP for %s did not match the columns present", gse)
  }

  pdata <- attr(meta, "pdata")
  # strip punctuation/case so "CAD_01 [EP00004]" and "EP00004" can meet
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(as.character(x)))
  # additionally drop zero-padding so "EC_013" matches "EC13" -- only zeros
  # that follow a non-digit, so "EC103" is left intact
  norm0 <- function(x) gsub("(?<=[^0-9])0+(?=[0-9])", "", norm(x), perl = TRUE)
  # text inside brackets, e.g. "CAD_01 [EP00004]" -> "EP00004"
  bracketed <- function(x) {
    out <- sub(".*\\[([^]]+)\\].*", "\\1", as.character(x))
    ifelse(grepl("\\[", as.character(x)), out, NA_character_)
  }
  keys <- list(
    title            = meta$title,
    last_token       = sub("^.*\\s", "", meta$title),
    first_token      = sub("\\s.*$", "", meta$title),
    bracketed_title  = bracketed(meta$title),
    description      = if (!is.null(pdata) && !is.null(pdata$description))
      as.character(pdata$description) else NULL,
    supp_file        = if (!is.null(pdata) && !is.null(pdata$supplementary_file_1))
      basename(as.character(pdata$supplementary_file_1)) else NULL
  )
  # any characteristics field can also hold the sample identifier
  # (GSE131681's count columns are its char_donor_id values)
  for (cl in grep("^char_", colnames(meta), value = TRUE)) {
    keys[[cl]] <- as.character(meta[[cl]])
  }
  keys <- keys[!vapply(keys, is.null, logical(1))]

  # normalized variants of every key, tried after the exact ones
  keys <- c(keys,
            setNames(lapply(keys, norm),  paste0("norm_",  names(keys))),
            setNames(lapply(keys, norm0), paste0("norm0_", names(keys))))

  for (k in names(keys)) {
    if (is.null(keys[[k]])) next
    probe <- if (startsWith(k, "norm0_")) norm0(cn) else
             if (startsWith(k, "norm_"))  norm(cn)  else cn
    hit <- match(probe, keys[[k]])
    if (sum(!is.na(hit)) > 0.5 * length(cn)) {
      log_info("remapped %d/%d columns to GSM IDs via '%s'",
               sum(!is.na(hit)), length(cn), k)
      meta$count_column <- NA_character_
      meta$count_column[hit[!is.na(hit)]] <- cn[!is.na(hit)]
      colnames(counts) <- ifelse(is.na(hit), cn, meta$gsm_id[hit])
      return(list(counts = counts, metadata = meta, method = k))
    }
  }
  log_warn("could not map count columns to GSM IDs -- inspect manually")
  meta$count_column <- NA_character_
  list(counts = counts, metadata = meta, method = "FAILED")
}

# ---- value-scale classification ------------------------------------------

classify_values <- function(mat) {
  v <- as.matrix(mat); storage.mode(v) <- "double"
  fin <- v[is.finite(v)]
  if (!length(fin)) return(list(verdict = "UNKNOWN", is_raw_counts = FALSE))

  frac_int <- mean(abs(fin - round(fin)) < 1e-8)
  any_neg  <- any(fin < 0)
  csum     <- colSums(v, na.rm = TRUE)
  med      <- stats::median(csum)
  cv       <- if (mean(csum) > 0) stats::sd(csum) / mean(csum) else NA_real_
  near1e6  <- !is.na(cv) && cv < 0.01 && med > 9e5 && med < 1.1e6

  verdict <- if (any_neg) {
    "NORMALIZED/TRANSFORMED (negatives -> log-scale or batch-corrected)"
  } else if (near1e6) {
    "NORMALIZED (columns sum to ~1e6 -> TPM or CPM)"
  } else if (frac_int > 0.999) {
    "RAW COUNTS (all-integer, variable library sizes)"
  } else if (frac_int > 0.90) {
    "LIKELY RAW COUNTS (mostly integer; estimated counts e.g. RSEM/salmon)"
  } else {
    "NORMALIZED (non-integer -> TPM/FPKM/RPKM)"
  }

  list(verdict = verdict,
       is_raw_counts = grepl("^RAW|^LIKELY RAW", verdict),
       frac_int = frac_int, col_sums = csum, median_sum = med, cv = cv,
       any_negative = any_neg)
}

# ---- top-level entry point ------------------------------------------------

fetch_gse <- function(gse, save_counts = TRUE) {
  info <- ds_info(gse)
  log_step("FETCH %s  [%s / %s / %s]", gse, info$disease, info$use, info$arm)

  meta <- get_gse_metadata(gse)
  res  <- get_gse_counts(gse)

  if (is.null(res$counts)) {
    res <- get_per_sample_counts(gse, meta$gsm_id)
  }
  if (is.null(res$counts)) {
    log_warn("NO COUNT MATRIX obtained for %s -- manual inspection needed", gse)
    return(list(gse = gse, metadata = meta, counts = NULL,
                source = NA_character_, scale = NULL))
  }

  rm <- remap_columns_to_gsm(res$counts, meta, gse = gse)
  counts <- rm$counts
  meta   <- rm$metadata
  utils::write.csv(meta, file.path(PATHS$metadata, paste0(gse, "_metadata.csv")),
                   row.names = FALSE)

  scale <- classify_values(counts)
  log_info("value scale: %s", scale$verdict)

  if (save_counts) {
    path <- file.path(PATHS$processed, paste0(gse, "_counts.csv"))
    data.table::fwrite(cbind(gene_id = rownames(counts), counts), path)
    log_info("counts -> %s (%s genes x %d samples)",
             path, fmt_n(nrow(counts)), ncol(counts))
  }

  list(gse = gse, metadata = meta, counts = counts,
       source = res$source, scale = scale,
       id_type = detect_id_type(rownames(counts)))
}
