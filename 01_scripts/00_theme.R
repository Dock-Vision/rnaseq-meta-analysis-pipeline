# ---------------------------------------------------------------------------
# 00_theme.R -- working figure style for 06_figures/.
#
# docs/METHODS.md S6 requires ONE consistent visual system across the whole
# figure set, so every plotting script sources this file and nothing defines
# its own colours or theme inline.
#
#     source("01_scripts/00_theme.R")
#
# Also provides save_fig(), which writes matched PNG and PDF versions into
# 06_figures/ -- one to look at, one that stays sharp at any size.
#
# For the final, publication-grade figure set see 00_theme_publication.R.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- palettes -------------------------------------------------------------
# Colour-blind-safe (Okabe-Ito derived). Disease colours are fixed per disease
# so the same condition is the same colour in every figure of the set.

# Assigned from the dataset registry so the same condition keeps the same
# colour in every figure of the set, whatever conditions you happen to be
# analysing. Override individual entries after sourcing if you want specific
# colours: PAL_CONDITION["Sepsis"] <- "#A63D40"
OKABE_ITO <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7",
               "#E69F00", "#56B4E9", "#A63D40", "#7F7F7F")

PAL_CONDITION <- local({
  lv <- if (exists("CONDITION_LABELS")) names(CONDITION_LABELS) else character(0)
  if (!length(lv)) return(setNames(character(0), character(0)))
  setNames(rep(OKABE_ITO, length.out = length(lv)), lv)
})

# Case/control. Grey = control everywhere, so the eye reads it consistently.
PAL_GROUP <- c(
  "Control" = "#9AA0A6",
  "Case"    = "#C0392B"
)

# Direction of effect, used by volcano and dot plots.
PAL_DIRECTION <- c(
  "Up"   = "#C0392B",
  "Down" = "#2E6DB4",
  "NS"   = "#BFC4C9"
)

# Sequential scale for heatmaps / correlation (diverging, zero-centred).
PAL_DIVERGING <- c("#2166AC", "#F7F7F7", "#B2182B")

# ---- base theme -----------------------------------------------------------

BASE_SIZE <- 11
BASE_FAMILY <- ""   # leave default; avoids missing-font warnings on this box

theme_dv <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.25, colour = "grey90"),
      panel.border       = element_rect(colour = "grey30", linewidth = 0.5),
      axis.text          = element_text(colour = "grey20", size = rel(0.9)),
      axis.title         = element_text(colour = "grey10", face = "plain"),
      plot.title         = element_text(face = "bold", size = rel(1.15),
                                        hjust = 0, margin = margin(b = 4)),
      plot.subtitle      = element_text(colour = "grey35", size = rel(0.92),
                                        hjust = 0, margin = margin(b = 8)),
      plot.caption       = element_text(colour = "grey40", size = rel(0.78),
                                        hjust = 0, margin = margin(t = 8)),
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      legend.title       = element_text(size = rel(0.9)),
      strip.background   = element_rect(fill = "grey94", colour = "grey30",
                                        linewidth = 0.4),
      strip.text         = element_text(colour = "grey10", face = "bold",
                                        size = rel(0.9)),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

theme_set(theme_dv())

# ---- annotation helpers ---------------------------------------------------

# METHODS S6: every figure annotates n per group. Returns labels like
# "SLE\n(n=31)" for use as axis tick labels.
label_with_n <- function(group_vector) {
  tab <- table(group_vector)
  setNames(sprintf("%s\n(n=%d)", names(tab), as.integer(tab)), names(tab))
}

# Standard caption stating the statistical method, so no figure travels without
# its multiple-testing statement.
caption_stats <- function(extra = NULL,
                          method = "Wald test (DESeq2), Benjamini-Hochberg FDR",
                          shrink = "apeglm") {
  txt <- sprintf("Statistics: %s. log2FC shrunk with %s.", method, shrink)
  if (!is.null(extra)) txt <- paste(txt, extra)
  txt
}

# ---- saving ---------------------------------------------------------------

# Writes both PNG (300 dpi, for viewing) and PDF (vector, for publication).
# METHODS S6 -- consistent output geometry across the figure set.
save_fig <- function(plot, name, width = 7, height = 5, dpi = 300) {
  if (!exists("PATHS")) {
    stop("save_fig() needs PATHS from 00_config.R -- source it first.")
  }
  png_path <- file.path(PATHS$fig_png, paste0(name, ".png"))
  pdf_path <- file.path(PATHS$fig_pdf, paste0(name, ".pdf"))

  ggsave(png_path, plot, width = width, height = height, dpi = dpi,
         units = "in", bg = "white")
  ggsave(pdf_path, plot, width = width, height = height,
         units = "in", bg = "white", device = grDevices::cairo_pdf)

  message(sprintf("[fig] %s  (%.1f x %.1f in)", name, width, height))
  invisible(c(png = png_path, pdf = pdf_path))
}

message("[theme] theme_dv loaded; palettes: condition/group/direction/diverging")
