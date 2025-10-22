#!/usr/bin/env Rscript

# -------------------------------
# Combined Multi-Project Taxa Heatmap
# Visualize taxa abundance across multiple projects
# -------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(RColorBrewer)

# Get parameters from Snakemake
input_file <- snakemake@input[[1]]
output_linear <- snakemake@output$linear
output_log <- snakemake@output$log

min_projects <- snakemake@params$min_projects
min_total_reads <- snakemake@params$min_total_reads
min_identity <- snakemake@params$min_identity
top_n_taxa <- snakemake@params$top_n_taxa
color_scheme <- snakemake@params$color_scheme
width <- snakemake@params$width
height <- snakemake@params$height
cluster_projects <- snakemake@params$cluster_projects == "TRUE"
cluster_taxa <- snakemake@params$cluster_taxa == "TRUE"

cat("------------------------------------------------------------\n")
cat("Multi-Project Combined Heatmap\n")
cat("Input file:", input_file, "\n")
cat("Min projects:", min_projects, "\n")
cat("Min total reads:", min_total_reads, "\n")
cat("Min identity:", min_identity, "\n")
cat("Top N taxa:", top_n_taxa, "\n")
cat("------------------------------------------------------------\n")

# Load combined data
data <- read_csv(input_file, show_col_types = FALSE)

cat("Loaded", nrow(data), "rows from", length(unique(data$project)), "projects\n")

# Filter by minimum identity
data_filtered <- data %>%
  filter(!is.na(taxon), taxon != "", total_reads > 0)

cat("After removing NA taxa:", nrow(data_filtered), "rows\n")

# Create sample identifier combining project, core, depth, etc.
data_filtered <- data_filtered %>%
  mutate(
    sample_id = paste(project, core, depth, sampling_batch, isolation_batch, library, sep = "_"),
    sample_id = gsub("_NA", "", sample_id),  # Remove NA values
    sample_id = gsub("_+", "_", sample_id),  # Clean up multiple underscores
    sample_id = gsub("_$", "", sample_id)    # Remove trailing underscores
  )

# Aggregate reads per taxon per sample
data_agg <- data_filtered %>%
  group_by(project, sample_id, taxon, taxid, obitag_rank) %>%
  summarise(reads = sum(total_reads), .groups = "drop")

cat("Aggregated to", nrow(data_agg), "taxon-sample combinations\n")

# Filter taxa by minimum appearance across projects
taxa_summary <- data_agg %>%
  group_by(taxon) %>%
  summarise(
    n_projects = n_distinct(project),
    total_reads = sum(reads),
    .groups = "drop"
  ) %>%
  filter(n_projects >= min_projects, total_reads >= min_total_reads)

cat("Kept", nrow(taxa_summary), "taxa meeting filter criteria\n")

# Filter to top N taxa by total reads
top_taxa <- taxa_summary %>%
  arrange(desc(total_reads)) %>%
  head(top_n_taxa) %>%
  pull(taxon)

cat("Selecting top", length(top_taxa), "taxa for visualization\n")

# Filter data to top taxa
data_top <- data_agg %>%
  filter(taxon %in% top_taxa)

# Create matrix for heatmap (samples x taxa)
heatmap_matrix <- data_top %>%
  select(sample_id, taxon, reads) %>%
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

# Transpose so taxa are rows, samples are columns
heatmap_matrix <- t(heatmap_matrix)

cat("Created heatmap matrix:", nrow(heatmap_matrix), "taxa x", ncol(heatmap_matrix), "samples\n")

# Create annotation for samples (project)
sample_annotation <- data_top %>%
  select(sample_id, project) %>%
  distinct() %>%
  column_to_rownames("sample_id")

# Set color palette
if (color_scheme %in% rownames(brewer.pal.info)) {
  n_colors <- brewer.pal.info[color_scheme, "maxcolors"]
  n_projects <- length(unique(sample_annotation$project))
  if (n_projects <= n_colors) {
    project_colors <- brewer.pal(max(3, n_projects), color_scheme)
  } else {
    project_colors <- colorRampPalette(brewer.pal(n_colors, color_scheme))(n_projects)
  }
  names(project_colors) <- unique(sample_annotation$project)
  annotation_colors <- list(project = project_colors)
} else {
  annotation_colors <- NA
}

# Generate linear scale heatmap
cat("Generating linear scale heatmap...\n")
pdf(output_linear, width = width, height = height)
pheatmap(
  heatmap_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = cluster_taxa,
  cluster_cols = cluster_projects,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  annotation_col = sample_annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Combined Taxa Abundance Across Projects (Linear Scale)"
)
dev.off()

# Generate log10 scale heatmap
cat("Generating log10 scale heatmap...\n")
heatmap_matrix_log <- log10(heatmap_matrix + 1)

pdf(output_log, width = width, height = height)
pheatmap(
  heatmap_matrix_log,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = cluster_taxa,
  cluster_cols = cluster_projects,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  annotation_col = sample_annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Combined Taxa Abundance Across Projects (Log10 Scale)"
)
dev.off()

cat("Heatmaps generated successfully!\n")
cat("Linear scale:", output_linear, "\n")
cat("Log10 scale:", output_log, "\n")
cat("------------------------------------------------------------\n")
