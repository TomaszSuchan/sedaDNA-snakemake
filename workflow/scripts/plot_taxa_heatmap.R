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

log_transform <- as.logical(args[1])    # Set to TRUE for log10 transformation
top_n_taxa <- as.integer(args[2])       # Maximum number of taxa to show in heatmap
plot_width <- as.numeric(args[3])          # Width of output PDF
plot_height <- as.numeric(args[4])         # Height of output PDF
input_file <- args[5]
output_file <- args[6]

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

# -------------------------------
# Validate data structure
# -------------------------------
required_cols <- c("core", "depth", "taxon", "total_reads")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}


# Check for duplicates (multiple entries per core-depth-taxon)
dup_check <- df %>%
  count(core, depth, taxon) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  stop(
    paste0(
      "ERROR: Input table is not yet summarized by core-depth-taxon.\n",
      "Found ", nrow(dup_check), " duplicate entries.\n",
      "Please ensure clustering/aggregation is completed before this step."
    )
  )
}


# -------------------------------
# Prepare matrix for heatmap
# -------------------------------
heatmap_data <- df %>%
  mutate(sample = paste(core, depth, sep = "_")) %>%
  select(sample, taxon, total_reads) %>%
  pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0)

mat <- as.data.frame(heatmap_data)
rownames(mat) <- mat$taxon
mat <- mat %>% select(-taxon)

# Optional log-transform
if (log_transform) {
  mat <- log10(mat + 1)
  cat("Applied log10 transformation to reads.\n")
}

# Limit to top N taxa
if (nrow(mat) > top_n_taxa) {
  mat <- mat %>%
    mutate(total = rowSums(.)) %>%
    arrange(desc(total)) %>%
    slice_head(n = top_n_taxa) %>%
    select(-total)
  cat("Reduced to top", top_n_taxa, "taxa for visualization.\n")
}

# Extract depth from column names: assuming column format is "core_depth"
sample_order <- colnames(mat) %>%
  str_split_fixed("_", 2) %>%
  as.data.frame() %>%
  mutate(depth_num = as.numeric(V2)) %>%
  arrange(depth_num) %>%
  pull(V1)  # just get the core names, adjust if needed

# Actually reorder columns
mat <- mat[, order(colnames(mat))] 

# -------------------------------
# Draw and save heatmap
# -------------------------------
cat("Drawing heatmap...\n")

pdf(output_file, width = plot_width, height = plot_height)
pheatmap(
  mat,
  cluster_rows = TRUE,       # cluster taxa
  cluster_cols = FALSE,      # skip clustering columns
  scale = "row",
  fontsize = 8,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = ""
)
dev.off()
cat("Heatmap saved to:", output_file, "\n")