# ---------------------------------------------------------------------------
# 00_theme_publication.R -- publication-grade figure system, writing to
# 07_final_figures/.
#
# Deliberately separate from 00_theme.R, which stays as the working /
# exploratory style for 06_figures/. Nothing here writes to 06_figures/, so
# exploratory output can never be mistaken for a final figure.
#
# Design decisions, and why:
#
#   COLOUR is assigned by the JOB it does, not by taste:
#     * Control vs Case      -> categorical pair, blue #2a78d6 / red #e34948.
#                               Validated all-pairs: CVD dE 21.6, normal dE 32.3,
#                               both >= 3:1 on the light surface.
#     * Down / NS / Up       -> DIVERGING (polarity), so two poles + a neutral
#                               grey midpoint. NS is deliberately desaturated so
#                               it recedes; it is a midpoint, not a series.
#     * z-scores in heatmaps -> the same diverging pair, grey at zero.
#     * magnitude ramps      -> ONE hue, light -> dark. Never a rainbow.
#
#   TYPE: Lato throughout (humanist sans, multiple weights available system-
#   wide). Text always wears ink tokens, never a series colour.
#
#   MARKS: thin strokes, recessive axes, no panel gridlines competing with the
#   data, direct labels in preference to legends where only one series exists.
#
#   OUTPUT: 600 dpi PNG via ragg (correct font rasterisation) + vector PDF via
#   cairo_pdf, identical geometry.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

stopifnot(requireNamespace("ragg", quietly = TRUE))
stopifnot(requireNamespace("systemfonts", quietly = TRUE))

# ---- type ------------------------------------------------------------------

PUB_FONT <- if ("Lato" %in% systemfonts::system_fonts()$family) "Lato" else "sans"
PUB_BASE_SIZE <- 9    # journal single-column base; scaled per figure

# ---- ink tokens (text never wears a series colour) -------------------------

INK <- list(
  primary   = "#16161a",
  secondary = "#4a4a52",
  muted     = "#77777f",
  faint     = "#b9b9c0",
  rule      = "#d8d8de",
  surface   = "#ffffff",
  panel     = "#fbfbfc"
)

# ---- palettes (validated; see header) --------------------------------------

PUB_GROUP <- c(Control = "#2a78d6", Case = "#e34948")

# diverging: polarity, neutral grey midpoint
PUB_DIRECTION <- c(Down = "#2a78d6", NS = "#c2c2c8", Up = "#e34948")

PUB_DIVERGING <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")

# single-hue sequential ramp (blue 100 -> 700), for magnitude only
PUB_SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# tissue context: two levels only, so a validated categorical pair
PUB_TISSUE <- c(tissue = "#1baf7a", blood = "#4a3aa7")

# Boxplot fills: a neutral two-step grey. Red/blue are reserved across the
# figure set for DIRECTION (up/down in disease), so Case/Control must not reuse
# them. Lightness carries the distinction; group names appear on the axis and in
# the legend, so identity never rests on colour alone.
# Boxplot fills. The CASE box is coloured by DIRECTION, using the same
# diverging scheme as Figures 1/2/3/5, so colour carries the finding rather than
# merely separating two groups. Control stays a constant near-neutral so the eye
# reads it as the reference. NS grey is deliberately a step darker than Control
# so the two never blur together.
# Validated all-pairs (light surface): worst CVD dE 14.6, worst normal-vision
# dE 21.9 -- both clear. An earlier Control/NS pairing measured 13.7 normal
# vision, below the 15 floor, and was re-stepped rather than shipped.
PUB_BOX <- c(Control = "#f4f4f5",
             Up      = "#e34948",
             Down    = "#2a78d6",
             NS      = "#a8a8b3")
# darker companions for the jittered sample points
PUB_BOX_PT <- c(Control = "#8b8b95",
                Up      = "#9e2a29",
                Down    = "#1a4f92",
                NS      = "#63636d")

PUB_ACCENT   <- "#e34948"   # pooled estimates, emphasis
PUB_NEUTRAL  <- "#6f6f78"   # study-level marks

# ---- theme -----------------------------------------------------------------

theme_pub <- function(base_size = PUB_BASE_SIZE, base_family = PUB_FONT,
                      grid = c("y", "x", "both", "none")) {
  grid <- match.arg(grid)
  th <- theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text               = element_text(colour = INK$primary, family = base_family),
      plot.background    = element_rect(fill = INK$surface, colour = NA),
      panel.background   = element_rect(fill = INK$surface, colour = NA),
      panel.border       = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = INK$rule, linewidth = 0.25),
      axis.line.x        = element_line(colour = INK$secondary, linewidth = 0.35),
      axis.line.y        = element_line(colour = INK$secondary, linewidth = 0.35),
      axis.ticks         = element_line(colour = INK$secondary, linewidth = 0.3),
      axis.ticks.length  = unit(2.2, "pt"),
      axis.text          = element_text(colour = INK$secondary, size = rel(0.92)),
      axis.title         = element_text(colour = INK$primary, size = rel(1.0)),
      axis.title.x       = element_text(margin = margin(t = 6)),
      axis.title.y       = element_text(margin = margin(r = 6)),
      plot.title         = element_text(size = rel(1.30), face = "bold",
                                        colour = INK$primary, hjust = 0,
                                        margin = margin(b = 3)),
      plot.subtitle      = element_text(size = rel(0.97), colour = INK$secondary,
                                        hjust = 0, margin = margin(b = 9),
                                        lineheight = 1.18),
      plot.caption       = element_text(size = rel(0.80), colour = INK$muted,
                                        hjust = 0, margin = margin(t = 9),
                                        lineheight = 1.20),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      legend.title       = element_text(size = rel(0.90), colour = INK$secondary),
      legend.text        = element_text(size = rel(0.90), colour = INK$secondary),
      legend.key.size    = unit(9, "pt"),
      legend.margin      = margin(0, 0, 0, 0),
      strip.background   = element_blank(),
      strip.text         = element_text(size = rel(1.0), face = "bold",
                                        colour = INK$primary,
                                        margin = margin(b = 4)),
      panel.spacing      = unit(11, "pt"),
      plot.margin        = margin(11, 13, 9, 11)
    )
  if (grid == "y") th <- th + theme(panel.grid.major.x = element_blank())
  if (grid == "x") th <- th + theme(panel.grid.major.y = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major = element_blank())
  th
}

# ---- helpers ---------------------------------------------------------------

# Axis labels carry units; this keeps the wording identical across figures.
LAB_LFC   <- expression(log[2]~"fold change (shrunken, disease vs control)")
LAB_PADJ  <- expression(-log[10]~"BH-adjusted "*italic(p))
LAB_VST   <- "Normalised expression (VST, log₂ scale)"
LAB_R     <- "Pearson correlation coefficient (r)"

# "n = 31" style annotation, used on every panel per the brief.
n_lab <- function(x, n) sprintf("%s\n(n = %d)", x, n)

# Formats a p-value for display without ever printing "p = 0".
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-16, "< 1e-16",
                ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))))
}

# ---- output ----------------------------------------------------------------

PUB_DIR <- if (exists("PATHS")) PATHS$fig_pub else file.path(getwd(), "07_final_figures")
dir.create(PUB_DIR, showWarnings = FALSE)
dir.create(file.path(PUB_DIR, "png"), showWarnings = FALSE)
dir.create(file.path(PUB_DIR, "pdf"), showWarnings = FALSE)

# 600 dpi PNG through ragg (proper font hinting) + vector PDF, same geometry.
save_pub <- function(plot, name, width, height, dpi = 600) {
  png_path <- file.path(PUB_DIR, "png", paste0(name, ".png"))
  pdf_path <- file.path(PUB_DIR, "pdf", paste0(name, ".pdf"))

  ragg::agg_png(png_path, width = width, height = height, units = "in",
                res = dpi, background = INK$surface, scaling = 1)
  print(plot); invisible(grDevices::dev.off())

  grDevices::cairo_pdf(pdf_path, width = width, height = height,
                       family = PUB_FONT, bg = INK$surface)
  print(plot); invisible(grDevices::dev.off())

  message(sprintf("[pub] %-28s %.1f x %.1f in @ %d dpi  (%.1f MB png)",
                  name, width, height, dpi, file.size(png_path) / 1e6))
  invisible(c(png = png_path, pdf = pdf_path))
}

message(sprintf("[pub-theme] font=%s  output=%s", PUB_FONT, PUB_DIR))
