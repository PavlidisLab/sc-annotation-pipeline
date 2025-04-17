library(ggplot2)
library(cowplot)
library(patchwork)
library(ggExtra)

base_theme <- theme(
    text = element_text(size = 25),            # Base text size
    axis.title = element_text(size = 25),      # Axis title size
    axis.text = element_text(size = 25),       # Axis tick label size
    plot.title = element_text(size = 25, face = "bold")  # Plot title size
)

plot_and_save <- function(y, hue) {

  dotplot <- ggplot(query_df, aes_string(x = "log1p_n_genes_by_counts", y = y, color = hue)) +
    geom_point(alpha = 0.6, size = 0.8) +
    theme_minimal() +
    labs(x = "log1p_n_genes_by_counts", y = y, color = hue)
  
  histo <- ggMarginal(dotplot, type = "histogram", fill="lightgrey", bins=100, size=2)

  return(histo)
}

plot_jointplots <- function(query_df, study_name, sample_name) {

  counts <- plot_and_save("log1p_total_counts", "counts_outlier")
  mito <- plot_and_save("pct_counts_mito", "outlier_mito")
  ribo <- plot_and_save("pct_counts_ribo", "outlier_ribo")
  hb <- plot_and_save("pct_counts_hb", "outlier_hb")
  total <- plot_and_save("log1p_total_counts", "total_outlier")

  combined <- patchwork::wrap_plots(counts, mito, ribo, hb, nrow = 2) & base_theme
  output_dir <- file.path(study_name, sample_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = file.path(output_dir, "outlier_jointplots_mqc.png"), plot = combined, width = 10, height = 5)
}

args <- commandArgs(trailingOnly = TRUE)

query_file <- args[1]
study_name <- args[2]
sample_name <- args[3]

# Read the TSV
query_df <- read.delim(query_file, header = TRUE, sep = "\t")

# Run the plotting
plot_jointplots(query_df, study_name, sample_name)
