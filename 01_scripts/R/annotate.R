# ---------------------------------------------------------------------------
# R/annotate.R -- gene identifier handling.
#
# Datasets arrive with either Ensembl gene IDs (e.g. GSE112087) or HGNC
# symbols. Everything downstream (marker panels, xCell, figures) is written in
# terms of SYMBOLS, so this module provides one place to detect the incoming ID
# type and convert.
#
# Mapping uses org.Hs.eg.db (offline, versioned) rather than biomaRt, which is
# a live web service and would make results non-reproducible and run-dependent.
# The built map is cached to 02_metadata/gene_id_map.rds.
# ---------------------------------------------------------------------------

# ---- ID type detection ----------------------------------------------------

detect_id_type <- function(ids) {
  ids <- as.character(ids)
  frac_ens <- mean(grepl("^ENSG[0-9]{11}", ids))
  frac_enst <- mean(grepl("^ENST[0-9]{11}", ids))
  frac_num <- mean(grepl("^[0-9]+$", ids))
  if (frac_enst > 0.5) return("ensembl_transcript")
  if (frac_ens > 0.5) return("ensembl")
  if (frac_num > 0.5) return("entrez")
  "symbol"
}

# Salmon/kallisto output is transcript-level. Summing transcript counts per gene
# is the standard way to get gene-level counts when tximport is not in play.
# Requires a transcript->gene map; org.Hs.eg.db does not carry one, so we use
# the ENSEMBLTRANS key where available and drop unmapped transcripts.
collapse_transcripts_to_genes <- function(mat) {
  if (detect_id_type(rownames(mat)) != "ensembl_transcript") return(mat)
  require_pkg("AnnotationDbi", "transcript -> gene collapsing")

  tx <- strip_ensembl_version(rownames(mat))

  # EnsDb is an Ensembl-native transcript database and maps essentially all
  # ENST IDs. org.Hs.eg.db's ENSEMBLTRANS key is NCBI-derived and covers only
  # ~8% of a salmon index -- it silently dropped the genes of interest from all three
  # salmon-quantified EC datasets, so it is only a fallback here.
  map <- NULL
  if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
    map <- tryCatch({
      db <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
      m <- AnnotationDbi::select(db, keys = unique(tx), keytype = "TXID",
                                 columns = c("TXID", "SYMBOL"))
      data.frame(ENSEMBLTRANS = m$TXID, SYMBOL = m$SYMBOL,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      log_warn("EnsDb lookup failed (%s); falling back to org.Hs.eg.db",
               conditionMessage(e)); NULL })
  }
  if (is.null(map) || !nrow(map)) {
    require_pkg("org.Hs.eg.db", "transcript -> gene collapsing")
    map <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(tx),
                            keytype = "ENSEMBLTRANS", columns = "SYMBOL"),
      error = function(e) NULL)
  }
  if (is.null(map) || !nrow(map)) {
    log_warn("no transcript->gene map available; leaving matrix transcript-level")
    return(mat)
  }
  map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), ]
  map <- map[!duplicated(map$ENSEMBLTRANS), ]
  sym <- map$SYMBOL[match(tx, map$ENSEMBLTRANS)]
  keep <- !is.na(sym)
  log_info("transcript->gene: %s / %s transcripts mapped, collapsing to %s genes",
           fmt_n(sum(keep)), fmt_n(length(tx)), fmt_n(length(unique(sym[keep]))))
  if (sum(keep) < 1000) {
    log_warn("only %s transcripts mapped -- gene-level matrix may be unreliable",
             fmt_n(sum(keep)))
  }
  rowsum(as.matrix(mat)[keep, , drop = FALSE], group = sym[keep], reorder = FALSE)
}

# Strip Ensembl version suffixes ("ENSG00000151702.12" -> "ENSG00000151702").
strip_ensembl_version <- function(ids) sub("\\..*$", "", as.character(ids))

# ---- mapping table --------------------------------------------------------

.gene_map_cache <- NULL

# Builds (and caches) an Ensembl <-> symbol table for the whole genome.
get_gene_map <- function(refresh = FALSE) {
  if (!is.null(.gene_map_cache) && !refresh) return(.gene_map_cache)

  cache_path <- file.path(PATHS$metadata, "gene_id_map.rds")
  if (file.exists(cache_path) && !refresh) {
    map <- readRDS(cache_path)
    .gene_map_cache <<- map
    log_info("gene map loaded from cache (%s rows)", fmt_n(nrow(map)))
    return(map)
  }

  require_pkg("org.Hs.eg.db", "gene ID mapping")
  require_pkg("AnnotationDbi", "gene ID mapping")

  log_info("building Ensembl<->symbol map from org.Hs.eg.db %s ...",
           packageVersion("org.Hs.eg.db"))

  db <- org.Hs.eg.db::org.Hs.eg.db
  keys <- AnnotationDbi::keys(db, keytype = "ENSEMBL")
  map <- AnnotationDbi::select(db, keys = keys, keytype = "ENSEMBL",
                               columns = c("SYMBOL", "GENENAME", "ENTREZID"))
  map <- map[!is.na(map$SYMBOL) & !is.na(map$ENSEMBL), ]
  map <- map[!duplicated(map$ENSEMBL), ]
  colnames(map) <- c("ensembl", "symbol", "gene_name", "entrez")

  saveRDS(map, cache_path)
  .gene_map_cache <<- map
  log_info("gene map built and cached: %s genes -> %s",
           fmt_n(nrow(map)), cache_path)
  map
}

# ---- conversion -----------------------------------------------------------

# Convert a vector of IDs to symbols. Unmapped IDs come back as NA.
ids_to_symbols <- function(ids) {
  type <- detect_id_type(ids)
  if (type == "symbol") return(as.character(ids))
  map <- get_gene_map()
  if (type == "ensembl") {
    return(map$symbol[match(strip_ensembl_version(ids), map$ensembl)])
  }
  map$symbol[match(as.character(ids), map$entrez)]
}

# Resolve the analysis gene panels (genes of interest + markers) to whatever ID space a
# given count matrix uses. Returns a data.frame with one row per requested
# symbol and the matching rowname, so a missing gene is visible rather than
# silently dropped.
resolve_panel <- function(symbols, matrix_rownames) {
  type <- detect_id_type(matrix_rownames)
  out <- data.frame(symbol = symbols, id = NA_character_,
                    found = FALSE, stringsAsFactors = FALSE)

  if (type == "symbol") {
    hit <- match(symbols, matrix_rownames)
  } else {
    map <- get_gene_map()
    if (type == "ensembl") {
      target <- map$ensembl[match(symbols, map$symbol)]
      hit <- match(target, strip_ensembl_version(matrix_rownames))
    } else {
      target <- map$entrez[match(symbols, map$symbol)]
      hit <- match(target, matrix_rownames)
    }
  }
  out$id    <- matrix_rownames[hit]
  out$found <- !is.na(hit)

  missing <- out$symbol[!out$found]
  if (length(missing)) {
    log_warn("panel genes not present in matrix: %s",
             paste(missing, collapse = ", "))
  }
  out
}

# Collapse a matrix from Ensembl IDs to unique symbols, summing counts for
# multiple Ensembl IDs mapping to one symbol. Needed for xCell, which requires
# a symbol-indexed matrix.
matrix_to_symbols <- function(mat) {
  if (detect_id_type(rownames(mat)) == "symbol") return(mat)
  sym <- ids_to_symbols(rownames(mat))
  keep <- !is.na(sym) & nzchar(sym)
  log_info("symbol conversion: %s / %s rows mapped",
           fmt_n(sum(keep)), fmt_n(nrow(mat)))
  mat <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  if (anyDuplicated(sym)) {
    log_info("collapsing %s duplicate symbols by sum", fmt_n(sum(duplicated(sym))))
    mat <- rowsum(as.matrix(mat), group = sym, reorder = FALSE)
  } else {
    rownames(mat) <- sym
  }
  mat
}
