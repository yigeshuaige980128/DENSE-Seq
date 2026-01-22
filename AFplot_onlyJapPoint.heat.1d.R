# Load required packages
library(ggplot2)  
library(cowplot)   # Combine "title | legend"
library(grid)      # unit() to set legend bar size

# --------------------------
# Parameter settings
# --------------------------
setwd("E:/1.now/pollenPool/20251213/t65_kasalath")
getwd()

input_file  <- "somatic_comb.T65_Kasalath_7_amplificationFree.pairwiseAFcount.txt"
output_file <- "somatic_comb.T65_Kasalath_7_amplificationFree.pairwiseAFcount.txt.pdf"

pdf_width  <- 16
pdf_height <- 12
units      <- "in"
dpi        <- 300
point_size <- 0.5
n_columns  <- 3

# Bandwidth adjustment for 1D vertical (y) density (the larger, the smoother)
y_bw_adjust <- 1.0

# Number of equal-width x-bins per chromosome (for local grouping only)
x_bins_per_chr <- 2000
# Minimal number of points required per bin to compute local 1D density;
# otherwise fall back to chromosome-wide 1D density
min_bin_n <- 40

# Header layout fine-tuning
legend_top_pad        <- 0.30  # Relative top blank height (0–1) above the legend row
title_x               <- 0.56  # Title position in its column (0–1), 0.5 = centered
header_rel_widths     <- c(1, 0.34)  # Width ratio: title column : legend column
legend_bar_width_mm   <- 55   # Horizontal color bar length (mm)
legend_bar_height_mm  <- 4.5  # Horizontal color bar thickness (mm)

# --------------------------
# Data preparation
# --------------------------
df <- read.table(input_file, header = TRUE, sep = "", stringsAsFactors = FALSE)

# Standardize column names by position:
#   column 1      -> Chr
#   column 2      -> POS
#   last 4 columns -> CTRL_Parent1, CTRL_Parent2, POOL_Parent1, POOL_Parent2
n_col <- ncol(df)

colnames(df)[1] <- "Chr"
colnames(df)[2] <- "POS"
colnames(df)[(n_col - 3):n_col] <- c(
  "CTRL_Parent1",
  "CTRL_Parent2",
  "POOL_Parent1",
  "POOL_Parent2"
)

# Keep only Chr1–Chr12
keep_chr <- paste0("Chr", 1:12)
df <- subset(df, Chr %in% keep_chr)
df$Chr <- factor(df$Chr, levels = keep_chr)

# Convert numeric columns
num_cols <- c("POS", "CTRL_Parent1", "CTRL_Parent2", "POOL_Parent1", "POOL_Parent2")
df[num_cols] <- lapply(num_cols, function(col) {
  suppressWarnings(as.numeric(df[[col]]))
})

# Derived fields
df$POOL_total   <- df$POOL_Parent1 + df$POOL_Parent2
df$POOL_P1_freq <- ifelse(df$POOL_total > 0, df$POOL_Parent1 / df$POOL_total, NA_real_)
df$x_mb         <- df$POS / 1e6

# ==== Per-chromosome mean allele frequency + annotation ====
# 1) Per-chromosome mean allele frequency (AF)
chr_mean_af <- aggregate(
  POOL_P1_freq ~ Chr,
  data = df,
  FUN  = function(x) mean(x, na.rm = TRUE)
)
colnames(chr_mean_af)[colnames(chr_mean_af) == "POOL_P1_freq"] <- "mean_AF"

# 2) Per-chromosome minimal x position (used as left anchor for labels)
x_min_chr <- aggregate(
  x_mb ~ Chr,
  data = df,
  FUN  = function(x) min(x, na.rm = TRUE)
)
colnames(x_min_chr)[colnames(x_min_chr) == "x_mb"] <- "x_min"

# 3) Merge and set factor levels consistent with df
chr_mean_af <- merge(chr_mean_af, x_min_chr, by = "Chr", all.x = TRUE)
chr_mean_af$Chr <- factor(chr_mean_af$Chr, levels = levels(df$Chr))

# 4) Label text (two decimal places) and y position (slightly above dashed line, capped at 0.99)
chr_mean_af$label <- sprintf("%.2f", chr_mean_af$mean_AF)
chr_mean_af$y_lab <- pmin(chr_mean_af$mean_AF + 0.03, 0.99)
# ===============================================================

# --------------------------
# 1D density function (y only)
# --------------------------
.eval_y_density <- function(y, adjust = 1, from = 0, to = 1, n = 512) {
  yy <- y[is.finite(y)]
  if (length(yy) < 2) return(rep(NA_real_, length(y)))
  d  <- stats::density(yy, adjust = adjust, from = from, to = to, n = n)
  f  <- stats::approxfun(d$x, d$y, yleft = 0, yright = 0)
  out <- ifelse(is.finite(y), f(y), NA_real_)
  out
}

# First compute chromosome-wide 1D density (used as fallback)
df$y_dens_chr <- NA_real_
for (chr in levels(df$Chr)) {
  idx <- which(df$Chr == chr)
  dens_vals <- .eval_y_density(df$POOL_P1_freq[idx], adjust = y_bw_adjust, from = 0, to = 1)
  rng <- range(dens_vals, na.rm = TRUE)
  if (is.finite(rng[1]) && diff(rng) > 0) {
    df$y_dens_chr[idx] <- (dens_vals - rng[1]) / diff(rng)
  } else {
    df$y_dens_chr[idx] <- 0
  }
}

# Then compute local 1D density per x-bin (to avoid full horizontal bands)
df$x_bin <- NA_integer_
df$y_dens_local <- NA_real_

for (chr in levels(df$Chr)) {
  idx_chr <- which(df$Chr == chr)
  x_chr   <- df$x_mb[idx_chr]
  
  # Split each chromosome along x into equal-width bins
  bin_id  <- cut(
    x_chr,
    breaks = x_bins_per_chr,
    include.lowest = TRUE,
    labels = FALSE
  )
  df$x_bin[idx_chr] <- bin_id
  
  # For each (Chr, x_bin), compute 1D density on y and normalize to 0–1 within that bin
  bins <- sort(unique(bin_id))
  for (b in bins) {
    idx_bin <- idx_chr[bin_id == b]
    y_bin   <- df$POOL_P1_freq[idx_bin]
    n_bin   <- sum(is.finite(y_bin))
    
    if (n_bin >= min_bin_n) {
      dens_vals <- .eval_y_density(y_bin, adjust = y_bw_adjust, from = 0, to = 1)
      rng <- range(dens_vals, na.rm = TRUE)
      if (is.finite(rng[1]) && diff(rng) > 0) {
        df$y_dens_local[idx_bin] <- (dens_vals - rng[1]) / diff(rng)
      } else {
        df$y_dens_local[idx_bin] <- 0
      }
    } else {
      # Too few points: fall back to chromosome-wide 1D density to avoid noisy patterns
      df$y_dens_local[idx_bin] <- df$y_dens_chr[idx_bin]
    }
  }
}

# --------------------------
# Optional dark red vertical lines + left annotations (pf loci)
# (kept as comments for future use)
# --------------------------
# locus_positions <- data.frame(
#   Chr = factor(c("Chr3", "Chr12"), levels = levels(df$Chr)),
#   POS = c(7796992, 910354),                     # 7.80 Mb, 0.91 Mb
#   Label = c("", ""),
#   x_offset = c(-1.5, -1.5),                     # Shift label 1.5 Mb to the left
#   y_position = c(0.96, 0.96)                    # y position for the label (0–1)
# )
#
# darkred_col <- "#8B0000"  
# annot_line_size <- 1.2

# --------------------------
# Plotting (color from local vertical density)
# --------------------------
plot_title <- "T65_Kasalath_7_amplificationFree (T65 AF, set reference to NP7)"

af_plot_core <- ggplot(df, aes(x = x_mb, y = POOL_P1_freq)) +
  # Hollow points: shape = 21, fill = NA; stroke controls border width
  geom_point(
    aes(color = y_dens_local),
    shape = 21,
    fill = NA,
    stroke = 0.4,
    size = point_size,
    na.rm = TRUE
  ) +
  # Per-chromosome mean AF as black dashed line
  geom_hline(
    data = chr_mean_af,
    aes(yintercept = mean_AF),
    linetype  = "11",
    color     = "black",
    linewidth = 0.4
  ) +
  # Label mean AF text slightly above the dashed line at left endpoint
  geom_text(
    data = chr_mean_af,
    aes(x = x_min, y = y_lab, label = label),
    hjust = 0,    # Left-aligned
    vjust = 0,    # Bottom of text at given y (above the line)
    size  = 6
  ) +
  # Optional dark red vertical lines and labels for pf loci
  # geom_vline(
  #   data = locus_positions,
  #   aes(xintercept = POS/1e6),
  #   color = darkred_col, linetype = "solid", size = annot_line_size
  # ) +
  # geom_text(
  #   data = locus_positions,
  #   aes(x = POS/1e6 + x_offset, y = y_position, label = Label),
  #   color = darkred_col, size = 5, angle = 90, hjust = 1, vjust = 0.5
  # ) +
  labs(
    x     = "Genomic Position (Mb)",
    y     = "Allele Frequency (jap)",
    color = "Local Vertical-Density Colored Scatter\n(per Chr)"
  ) +
  scale_color_gradientn(
    colours = c("grey95", "#2ca02c", "#f1c40f", "#F28E2B"),
    values  = c(0, 0.50, 0.85, 1.00),
    oob     = scales::squish,
    na.value = "grey85"
  ) +
  facet_wrap(~ Chr, ncol = n_columns, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none",   # Hide legend in the main plot (legend extracted separately)
    axis.text = element_text(size = 8),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

# Extract legend (horizontal color bar) to be placed on the right
legend_grob <- cowplot::get_legend(
  af_plot_core +
    guides(
      color = guide_colorbar(
        title.position = "top",
        direction = "horizontal",                 # Horizontal color bar
        label.position = "bottom",
        barwidth = unit(legend_bar_width_mm, "mm"),
        barheight = unit(legend_bar_height_mm, "mm")
      )
    ) +
    theme(
      legend.position = "bottom",                 # For legend extraction only
      legend.justification = "center",
      legend.box.margin = margin(t = 2, r = 2, b = 2, l = 2),
      legend.title = element_text(hjust = 0.5)
    )
)

# Top blank spacer: an empty ggplot as placeholder
spacer_plot <- ggplot() + theme_void()

# Right column: blank space + horizontal legend
legend_col <- plot_grid(
  spacer_plot,
  legend_grob,
  ncol = 1,
  rel_heights = c(legend_top_pad, 1)  # Relative height of top blank space
)

# Title (slightly shifted to appear visually centered)
title_grob <- ggdraw() + 
  draw_label(plot_title, x = title_x, hjust = 0.5, fontface = "bold", size = 20)

# Combine title (left) and legend (right) into the same row
header_row <- plot_grid(
  title_grob, legend_col,
  ncol = 2, align = "h",
  rel_widths = header_rel_widths
)

# Combine header and main plot
final_plot <- plot_grid(
  header_row,
  af_plot_core,
  ncol = 1,
  rel_heights = c(0.12, 1)  # Relative height of header row
)

# Export to PDF
ggsave(
  output_file,
  plot   = final_plot,
  width  = pdf_width,
  height = pdf_height,
  units  = units,
  dpi    = dpi,
  limitsize = FALSE
)
