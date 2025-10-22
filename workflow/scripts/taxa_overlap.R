#!/usr/bin/env Rscript

# -------------------------------
# Taxa Overlap Analysis
# Generate Venn diagram showing shared taxa across projects
# -------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(ggvenn)

# Get parameters from Snakemake
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
min_identity <- snakemake@params$min_identity
min_reads <- snakemake@params$min_reads

cat("------------------------------------------------------------\n")
cat("Taxa Overlap Analysis\n")
cat("Input file:", input_file, "\n")
cat("Min identity:", min_identity, "\n")
cat("Min reads:", min_reads, "\n")
cat("------------------------------------------------------------\n")

# Load combined data
data <- read_csv(input_file, show_col_types = FALSE)

cat("Loaded", nrow(data), "rows from", length(unique(data$project)), "projects\n")

# Filter and aggregate
taxa_by_project <- data %>%
  filter(!is.na(taxon), taxon != "") %>%
  group_by(project, taxon) %>%
  summarise(total_reads = sum(total_reads), .groups = "drop") %>%
  filter(total_reads >= min_reads)

cat("After filtering:", nrow(taxa_by_project), "taxon-project combinations\n")

# Create list of taxa sets per project
taxa_sets <- taxa_by_project %>%
  group_by(project) %>%
  summarise(taxa = list(unique(taxon)), .groups = "drop")

# Convert to named list for ggvenn
taxa_list <- setNames(taxa_sets$taxa, taxa_sets$project)

# Count unique taxa per project
for (proj in names(taxa_list)) {
  cat(proj, ":", length(taxa_list[[proj]]), "unique taxa\n")
}

# Generate Venn diagram
# Note: ggvenn supports up to 4 sets optimally
n_projects <- length(taxa_list)

if (n_projects > 5) {
  cat("Warning: More than 5 projects detected. Venn diagram may be complex.\n")
  cat("Showing only first 5 projects.\n")
  taxa_list <- taxa_list[1:5]
}

pdf(output_file, width = 10, height = 10)

if (n_projects <= 1) {
  # Just show a barplot
  plot.new()
  text(0.5, 0.5, paste("Only", n_projects, "project found.\nVenn diagram requires at least 2 projects."),
       cex = 1.5, col = "red")
} else {
  p <- ggvenn(
    taxa_list,
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF")[1:min(n_projects, 5)],
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 4
  ) +
    ggtitle("Shared Taxa Across Projects") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

  print(p)
}

dev.off()

# Calculate overlap statistics
overlap_stats <- data.frame(
  project = names(taxa_list),
  unique_taxa = sapply(taxa_list, length)
)

# Find shared taxa across all projects
if (n_projects > 1) {
  shared_all <- Reduce(intersect, taxa_list)
  cat("\nTaxa shared across all projects:", length(shared_all), "\n")

  if (length(shared_all) > 0 && length(shared_all) <= 20) {
    cat("Shared taxa:\n")
    print(shared_all)
  }
}

cat("\nVenn diagram saved to:", output_file, "\n")
cat("------------------------------------------------------------\n")
