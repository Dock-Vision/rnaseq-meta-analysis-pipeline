# Bulk RNA-seq Meta-Analysis Pipeline

**A reproducible R pipeline for asking whether a gene is genuinely dysregulated
across many published bulk RNA-seq studies — or whether the signal is a shift in
cell composition.**

Built and maintained by **DockVision**.

Point it at a list of GEO accessions and a gene of interest. It downloads each
series, works out what the data actually is, runs differential expression to a
publication standard, tests every result against cell-identity markers and
deconvolution, pools the studies with a formal random-effects meta-analysis, and
writes a consistent, self-documenting figure set.

---

## The problem this pipeline is built around

A single differential expression result is easy. A cross-study meta-analysis of
public bulk RNA-seq is where things go quietly wrong, and this pipeline is
organised around the failure modes rather than around the happy path:

| Failure mode | What it does to your result | How this pipeline handles it |
|---|---|---|
| **Cell composition** | A gene "falls in disease" because the cells expressing it are scarcer — not because it is regulated | Every result reported beside a marker panel, a covariate-adjusted refit, a marker-correlation test, and xCell deconvolution (S2) |
| **Technical replicates** | Lane-level GSMs inflate *n* and shrink p-values | Replicate keys detected and collapsed before any test (`R/prepare.R`) |
| **Normalised matrices** | DESeq2 silently accepts FPKM and returns nonsense | Value scale classified on download; normalised series routed to limma-trend + ashr |
| **Multi-arm series** | A "control vs everything else" rule pools unrelated diseases into the case group | Explicit `CONTRAST_SPEC`; unmatched samples are dropped, never absorbed |
| **Unshrunken effect sizes** | Noisy low-count LFCs distort inverse-variance pooling | `lfcShrink()` (apeglm/ashr) everywhere, estimator recorded |
| **Naive pooling** | Two opposite effects average into a null | `metafor` random-effects with τ², I², Q, plus stratified estimates and a moderator test |
| **Batch structure** | "Corrected for batch" with no evidence either way | PCA before/after **and** % variance by batch vs biology, both quantified |
| **Detection floor** | A barely-expressed gene produces a confident, meaningless effect | Dedicated sensitivity stage (`40_`) that re-tests on detected samples only |

The full standard, with the reasoning behind each rule, is in
**[`docs/METHODS.md`](docs/METHODS.md)**. The code refers to it by section
number, so any comment in a script can be traced to the rule it implements.

---

## Workflow

```
  config/datasets.csv          your studies, one row each
  01_scripts/00_config.R       your genes, your marker panels
            │
            ▼
  ┌─────────────────────────────────────────────────────────────┐
  │ 01_fetch_dataset.R      inspect ONE series before committing │
  │                         → is it counts or TPM? how many real │
  │                           samples? how many donors?          │
  └─────────────────────────────────────────────────────────────┘
            │
            ▼
  ┌─────────────────────────────────────────────────────────────┐
  │ 10_run_dataset.R        one dataset, end to end   (per GSE)  │
  │                                                              │
  │   fetch ─ collapse technical replicates ─ assign groups      │
  │         ─ filter ─ DESeq2 (or limma-trend) + LFC shrinkage   │
  │         ─ cell-composition: adjusted model, marker           │
  │           correlation, xCell deconvolution                   │
  │         ─ PCA + variance decomposition                       │
  │         ─ figures ─ sessionInfo ─ run-log row                │
  └─────────────────────────────────────────────────────────────┘
            │
            │   all four read 10_run_dataset.R's outputs; none depends
            │   on the others, so run whichever you need
            │
            ├──────────────────┬──────────────────┬──────────────────┐
            ▼                  ▼                  ▼                  ▼
  ┌──────────────────┐ ┌────────────────┐ ┌────────────────┐ ┌────────────────┐
  │ 20_meta_         │ │ 30_coexpres-   │ │ 40_sensitivity │ │ 60_celltype_   │
  │    analysis.R    │ │    sion_atlas.R│ │    _tests.R    │ │    atlas.R     │
  │                  │ │                │ │                │ │                │
  │ random-effects   │ │ baseline       │ │ is the effect  │ │ the gene vs    │
  │ pooling;         │ │ samples only,  │ │ real, or a     │ │ EVERY xCell    │
  │ τ², I², Q;       │ │ ComBat pro-    │ │ detection-     │ │ cell type,     │
  │ tissue-strati-   │ │ tecting cell   │ │ floor artifact?│ │ pooled across  │
  │ fied subgroups   │ │ type; within-  │ │ three inde-    │ │ cohorts — the  │
  │ + moderator      │ │ study corr.,   │ │ pendent tests  │ │ whole cell-type│
  │ test             │ │ then pooled    │ │                │ │ space, not one │
  │                  │ │                │ │                │ │ marker panel   │
  └──────────────────┘ └────────────────┘ └────────────────┘ └────────────────┘
            │                  │                                     │
            └────────┬─────────┘                                     │
                     ▼                                               ▼
  ┌──────────────────────────────┐              ┌──────────────────────────────┐
  │ 50_publication_figures.R     │              │ 61_celltype_atlas_figure.R   │
  │ Figures 1–7                  │              │ Figures 8–9 (one per gene)   │
  └──────────────────────────────┘              └──────────────────────────────┘
                     │                                          │
                     └────────────────────┬─────────────────────┘
                                          ▼
                        the nine-figure set, 600 dpi PNG
                                + vector PDF
```

### The nine figures

Figures 1–7 come from `50_publication_figures.R`; Figures 8–9 from
`61_celltype_atlas_figure.R`, one per gene in `GENES_OF_INTEREST`.

| # | Figure | Answers |
|---|---|---|
| 1 | Dot plot | What is the effect size in every study, side by side? |
| 2 | Volcano | Where do the genes of interest sit against the whole transcriptome? |
| 3 | Forest | What does the pooled random-effects estimate say, overall and by tissue? |
| 4 | Boxplots | What do the underlying expression distributions actually look like? |
| 5 | Marker heatmap | Regulation, or composition shift? |
| 6 | PCA | What is the sample structure, before and after batch correction? |
| 7 | Correlation | Is the gene's normal co-expression relationship preserved? |
| 8–9 | Cell-type association atlas | Across the *whole* cell-type space, which populations does each gene track with, and how consistently across cohorts? |

Figures 8–9 are the marker panel of Figure 5 generalised: instead of four
hand-picked genes, every cell type xCell can score, pooled across cohorts with
heterogeneity statistics. Read the caveat in the caption — a positive *r* means
*samples richer in that cell type express more of the gene*, which is an
association across samples, **not** expression measured inside that cell type.

---

## Quick start

```bash
# 1. clone
git clone https://github.com/Dock-Vision/rnaseq-meta-analysis-pipeline.git
cd rnaseq-meta-analysis-pipeline

# 2. install packages (system R, OUTSIDE any conda environment)
Rscript 01_scripts/00_install_packages.R
Rscript 01_scripts/00b_install_extras.R      # org.Hs.eg.db + xCell

# 3. tell the pipeline what to analyse
cp config/datasets.example.csv config/datasets.csv
$EDITOR config/datasets.csv                  # your GEO accessions
$EDITOR 01_scripts/00_config.R               # GENES_OF_INTEREST, MARKERS

# 4. look before you leap -- inspect one series
Rscript 01_scripts/01_fetch_dataset.R GSE112087

# 5. run it
Rscript 01_scripts/10_run_dataset.R GSE112087 2>&1 | tee 08_logs/GSE112087.log
Rscript 01_scripts/10_run_dataset.R GSE213001 2>&1 | tee 08_logs/GSE213001.log
#   ... one call per dataset

Rscript 01_scripts/20_meta_analysis.R        # pool the case/control studies
Rscript 01_scripts/30_coexpression_atlas.R   # baseline co-expression atlas
Rscript 01_scripts/40_sensitivity_tests.R    # optional: stress-test a result
Rscript 01_scripts/60_celltype_atlas.R       # gene vs cell-type association atlas

Rscript 01_scripts/50_publication_figures.R  # Figures 1-7
Rscript 01_scripts/61_celltype_atlas_figure.R # Figures 8-9
```

**Always run from the project root**, not from inside `01_scripts/` — the
scripts resolve paths relative to the working directory and will stop with a
clear error if you get this wrong.

Everything is re-runnable off the download cache in `03_raw_data/`, so a second
run costs minutes rather than a fresh download of every supplementary file.

### Configuring it for your own question

Three things, and nothing else:

1. **`GENES_OF_INTEREST`** in `01_scripts/00_config.R` — one or two HGNC
   symbols. With two, the first is the primary gene and the second its
   comparator, which is what the co-expression atlas correlates.
2. **`MARKERS`** in the same file — the cell-identity panel for the cell type
   that carries your gene. This is what licenses the regulation-versus-
   composition call, so it is worth choosing carefully.
3. **`config/datasets.csv`** — one row per GEO series. Columns are documented
   in `config/datasets.example.csv`.

Two optional escape hatches, both in `00_config.R` with worked examples:
`CONTRAST_SPEC` for series with more than two phenotype levels, and
`MANUAL_COLUMN_MAP` for series whose count-matrix headers cannot be matched to
GEO sample IDs by any string rule.

---

## Repository layout

```
01_scripts/
  00_config.R                 paths, gene panels, dataset registry  ← EDIT THIS
  00_install_packages.R       package installer
  00b_install_extras.R        org.Hs.eg.db + xCell (the awkward two)
  00_theme.R                  working figure style      → 06_figures/
  00_theme_publication.R      publication figure style  → 07_final_figures/
  01_fetch_dataset.R          inspect one GEO series before committing to it
  10_run_dataset.R            the per-dataset pipeline (the workhorse)
  20_meta_analysis.R          random-effects pooling across studies
  30_coexpression_atlas.R     baseline co-expression atlas
  40_sensitivity_tests.R      detection-floor and compartment stress tests
  50_publication_figures.R    Figures 1-7 of the final set
  60_celltype_atlas.R         gene vs every xCell cell type, pooled
  61_celltype_atlas_figure.R  Figures 8-9 of the final set
  R/
    utils.R                   logging, sessionInfo capture, run log
    annotate.R                Ensembl ↔ symbol mapping (offline, cached)
    download.R                GEO fetch with caching and format detection
    prepare.R                 replicate collapsing, filtering, PCA, variance
    de.R                      DESeq2 / limma-trend, shrinkage, composition
    deconv.R                  xCell deconvolution
    meta.R                    metafor random-effects pooling
    figures.R                 the seven core figure types

config/
  datasets.example.csv        the dataset registry format, documented

docs/
  METHODS.md                  the analysis standard the code enforces
```

Output folders (`02_metadata/` … `08_logs/`) are created on first run and are
git-ignored — data and results are never committed.

---

## Tools used

Developed and tested on **R 4.5.2** with **Bioconductor 3.22**.

| Package | Role |
|---|---|
| [GEOquery](https://bioconductor.org/packages/GEOquery/) | GEO series and sample metadata retrieval |
| [DESeq2](https://bioconductor.org/packages/DESeq2/) | Differential expression on raw counts |
| [limma](https://bioconductor.org/packages/limma/) | limma-trend path for series shipped as TPM/FPKM |
| [apeglm](https://bioconductor.org/packages/apeglm/) · [ashr](https://cran.r-project.org/package=ashr) | Log2 fold-change shrinkage |
| [metafor](https://cran.r-project.org/package=metafor) | Random-effects meta-analysis (REML, Fisher *z*) |
| [xCell](https://github.com/dviraran/xCell) | Cell-type enrichment / deconvolution |
| [sva](https://bioconductor.org/packages/sva/) | ComBat batch correction with biology protected |
| [variancePartition](https://bioconductor.org/packages/variancePartition/) | Variance attributable to batch vs biology |
| [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/) | Offline, versioned gene ID mapping |
| [ggplot2](https://ggplot2.tidyverse.org/) · patchwork · ggrepel · ragg | Figures |
| [data.table](https://rdatatable.gitlab.io/data.table/) | Fast matrix I/O |

`org.Hs.eg.db` is used for identifier mapping in preference to a live biomaRt
query, because a web service that changes between runs makes results
irreproducible.

**System dependency:** xCell pulls in GSVA, which needs ImageMagick headers. On
Debian/Ubuntu, install them before running the installer:

```bash
sudo apt-get install -y libmagick++-dev
```

Install packages with **system R, outside any conda environment**. Installing
the Bioconductor stack through conda shadows the system library and breaks
DESeq2/sva/GSVA in ways that are painful to diagnose.

---

## Reproducibility

- One seed (`SEED` in `00_config.R`), set once and inherited by every stochastic
  step.
- Every script writes `sessionInfo()` — with the seed, R version and library
  paths — to `08_logs/`.
- Every completed step appends a row to `08_logs/RUN_LOG.md` automatically,
  recording dimensions, group sizes, the value scale and anything unexpected.
  Re-running a dataset updates its row rather than appending a near-duplicate,
  so the log stays a true description of the current outputs.
- Raw downloads are kept exactly as downloaded in `03_raw_data/`; derived tables
  go to `04_processed/`. Every step is reproducible from source.

---

## Citation

If this pipeline contributed to work you are publishing, please cite it:

> DockVision. *Bulk RNA-seq Meta-Analysis Pipeline*: a reproducible R workflow
> for cross-study differential expression, cell-composition control and
> random-effects meta-analysis. Version 1.0.0, 2026.
> Available at: `https://github.com/Dock-Vision/rnaseq-meta-analysis-pipeline`

BibTeX:

```bibtex
@software{dockvision_rnaseq_meta_2026,
  author  = {{DockVision}},
  title   = {Bulk RNA-seq Meta-Analysis Pipeline: a reproducible R workflow
             for cross-study differential expression, cell-composition control
             and random-effects meta-analysis},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/Dock-Vision/rnaseq-meta-analysis-pipeline},
  note    = {Accessed: <date>}
}
```

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff), which GitHub
renders as a "Cite this repository" button on the repository page.

**Please also cite the underlying methods** — DESeq2, apeglm, metafor, xCell,
sva and limma each have their own publication, and the meta-analysis rests on
them.

---

## Contributing and support

Issues and pull requests are welcome. For collaboration, custom analysis or
consulting enquiries, please open an issue on this repository.

## License

[MIT](LICENSE) © 2026 DockVision.

The pipeline is licensed permissively, but the **data** it downloads is not
covered by this licence: GEO datasets carry their own terms, and every study you
analyse should be cited in your own work.
