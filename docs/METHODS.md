# Analysis standards

The rules this pipeline enforces, and why each one is there. The code refers to
these sections by number (`METHODS S4`, `METHODS S2.3(b)`, and so on), so a
comment in a script can always be traced to the standard it implements.

These are not style preferences. Each one closes a specific way that a bulk
RNA-seq meta-analysis produces a confident, publishable, wrong answer.

---

## S1. Statistics

**S1.1 — Report the raw p-value *and* the FDR-adjusted p-value, and name the
method.**
Every results table carries both `pvalue` and `padj` (Benjamini–Hochberg), and
every figure caption states the correction. A raw p-value presented alone from
a transcriptome-wide test of ~15,000 genes is not interpretable.

**S1.2 — Report shrunken effect sizes, never raw DESeq2 log2 fold changes.**
`lfcShrink()` with `type = "apeglm"`, falling back to `"ashr"` where the
contrast cannot be expressed as a model coefficient. Raw log2FC estimates for
low-count genes are dominated by noise, and a meta-analysis weights by the
inverse of their standard errors — so unshrunken estimates do not merely add
scatter, they systematically distort the pooled result. The script records
which estimator was used for each dataset.

---

## S2. Cell-composition control

**The central scientific control of this pipeline.**

In bulk RNA-seq a gene can fall simply because the cells expressing it are
scarcer in the sample. A drop in an endothelial transcription factor in
diseased lung may mean endothelial dropout, not downregulation per cell. These
two explanations have opposite biological meanings and identical differential
expression results.

**S2.1 — Every result is reported next to a cell-identity marker panel.**
Configured in `MARKERS` (`00_config.R`); `tissue_class` in the registry decides
which panel a dataset gets.

**S2.2 — The decision rule is stated, not implied.**
Gene down *and* markers all down together → likely composition shift. Gene down
*while* markers stay flat → likely genuine regulation. The pipeline writes the
verdict into the results and the run log.

**S2.3 — The relationship is quantified, never eyeballed.**
Two independent tests, both reported with statistics:
  - **(a)** refit the differential expression model with a marker score as a
    covariate, and report whether the effect survives adjustment;
  - **(b)** correlate the gene's change against the marker change across
    samples, with a test statistic, confidence interval and p-value.

**S2.4 — Run real deconvolution.**
xCell enrichment scores for every bulk dataset, reported as supplementary
results with a between-group test. A four-gene marker panel says "endothelial
genes moved"; deconvolution says by how much, with a p-value.

---

## S3. Batch correction and structure

**Show PCA before *and* after correction, and quantify what changed.** The
percentage of variance attributable to batch versus biological group is
reported for both states (PC regression against each covariate). "Batch was
corrected" without that quantification is an assertion, not a result — and
over-correction that removes real biological signal looks identical to success
unless you measure it.

Where cell type or another biological variable must be preserved, it is passed
through ComBat's `mod` argument so the correction cannot regress it away.

---

## S4. Meta-analysis

**Pool with formal random-effects models (`metafor`), never by averaging.**

Correlations are pooled on the Fisher *z* scale and back-transformed; effect
sizes are pooled on the log2 fold-change scale, inverse-variance weighted, REML.

**τ², I², Q and the pooled estimate with its confidence interval are all part of
the result**, not optional diagnostics. Where heterogeneity is high — as it
usually is across tissues and platforms — a single pooled estimate is reported
*with* the stratified estimates and a formal moderator test, because averaging
two opposite effects into a null is the most common way a meta-analysis
misleads.

Studies where a gene fell below the expression filter are reported as
*unmeasured*, distinct from *measured and null*.

---

## S5. Reproducibility

- Every script writes `sessionInfo()` to `08_logs/`, with the seed, the library
  paths and the R version.
- One project seed (`SEED` in `00_config.R`) set once and inherited by every
  stochastic step.
- Raw downloads are kept as-downloaded in `03_raw_data/`; derived tables are
  written to `04_processed/`. Every step is reproducible from source rather
  than from an intermediate someone has since edited.
- Every completed step appends a row to `08_logs/RUN_LOG.md` automatically,
  carrying the numbers a reader would otherwise have to recompute.

---

## S6. Figures

One shared visual system across the entire figure set — theme, palette, font
size and mark scale defined once in `00_theme.R` (working figures) and
`00_theme_publication.R` (final figures), and sourced everywhere.

Every figure carries:
- axis labels **with units**,
- annotated **n per group**,
- a caption stating the test, the multiple-testing correction and the shrinkage
  estimator.

A figure that leaves its own statistics to the surrounding prose is not
finished, because figures are what get reused, screenshotted and presented
separately from the text that qualified them.

Colour-blind-safe palettes throughout, and colour is assigned by the job it
does: categorical pairs for case/control, a diverging scale with a neutral
midpoint for direction of effect, a single-hue ramp for magnitude.

---

## S7. Interpretation

Every written interpretation must address biological plausibility, compare
against existing literature, and **state limitations explicitly** — cohort
imbalance, composition confounding, in vitro versus in vivo, absent controls,
platform and batch heterogeneity, multiple-testing burden, and any gene that
was near its detection floor in the tissue tested.

---

## S8. Baseline co-expression atlas

Applies to `30_coexpression_atlas.R`, which asks whether two genes' normal
relationship is preserved.

**S8.1 — Baseline samples only.** Untreated, healthy, static, vehicle-control.
Every stimulated condition is excluded — cytokines, hypoxia, shear stress,
knockdown, drug treatment. The selection is recorded per sample in an audit
table rather than assumed, because "healthy control" in a GEO series regularly
means "the untreated arm of a stimulation experiment".

**S8.2 — Batch correction protecting cell type.** ComBat with the cell-type
variable in `mod`, so genuine differences between cell types are not removed
along with the study batch.

**S8.3 — Correlate within each study first, then pool.** Pooling raw samples
across studies and correlating once manufactures correlation out of
between-study batch structure. Per-study correlations are computed first, then
combined by random-effects meta-analysis (S4).
