#!/usr/bin/env Rscript

# -------------------------------
# Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(stringr)
})

# -------------------------------
# Parameters
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_taxa_heatmap.R <clustered_taxa_table.csv> <output_pdf>")
}

input_file <- args[1]
output_file <- args[2]

cat("------------------------------------------------------------\n")
cat("Generating taxa heatmap\n")
cat("Input file: ", input_file, "\n")
cat("Output file: ", output_file, "\n")
cat("------------------------------------------------------------\n")

# -------------------------------
# Load clustered taxa table
# -------------------------------
df <- read_csv(input_file, show_col_types = FALSE)
cat("Loaded", nrow(df), "lines.\n")

# Expecting columns: taxid, taxon, rank, n_sequences, sequences
# Split sequences into individual sequence IDs
df_long <- df %>%
  separate_rows(sequences, sep = ";") %>%
  distinct()

# -------------------------------
# Count abundance (per taxon)
# -------------------------------
abundance <- df_long %>%
  count(taxon = !!sym(names(df_long)[2]), sequence_id = sequences) %>%
  pivot_wider(names_from = taxon, values_from = n, values_fill = 0)

# If abundance matrix is too sparse, limit to most frequent taxa
if (ncol(abundance) > 50) {
  taxa_sums <- colSums(abundance[,-1])
  top_taxa <- names(sort(taxa_sums, decreasing = TRUE))[1:50]
  abundance <- abundance %>%
    select(sequence_id, all_of(top_taxa))
  cat("Reduced to top 50 taxa for visualization.\n")
}

# -------------------------------
# Convert to matrix and scale
# -------------------------------
mat <- as.data.frame(abundance)
rownames(mat) <- mat$sequence_id
mat <- mat[,-1]

# Scale across taxa for visualization
mat_scaled <- t(scale(t(as.matrix(mat))))

# -------------------------------
# Plot heatmap
# -------------------------------
pdf(output_file, width = 10, height = 8)
pheatmap(
  mat_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  border_color = NA,
  main = "Taxon Abundance Heatmap",
  fontsize = 9
)
dev.off()

cat("✅ Heatmap saved to:", output_file, "\n")
cat("------------------------------------------------------------\n")