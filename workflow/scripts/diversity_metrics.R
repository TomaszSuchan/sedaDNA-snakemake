#!/usr/bin/env Rscript

# -------------------------------
# Diversity Metrics Calculation
# Calculate alpha and beta diversity across projects
# -------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(vegan)

# Get parameters from Snakemake
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
min_identity <- snakemake@params$min_identity

cat("------------------------------------------------------------\n")
cat("Diversity Metrics Calculation\n")
cat("Input file:", input_file, "\n")
cat("Min identity:", min_identity, "\n")
cat("------------------------------------------------------------\n")

# Load combined data
data <- read_csv(input_file, show_col_types = FALSE)

cat("Loaded", nrow(data), "rows from", length(unique(data$project)), "projects\n")

# Filter and prepare data
data_filtered <- data %>%
  filter(!is.na(taxon), taxon != "", total_reads > 0)

# Create sample identifier
data_filtered <- data_filtered %>%
  mutate(
    sample_id = paste(project, core, depth, sampling_batch, isolation_batch, library, sep = "_"),
    sample_id = gsub("_NA", "", sample_id),
    sample_id = gsub("_+", "_", sample_id),
    sample_id = gsub("_$", "", sample_id)
  )

# Aggregate reads per taxon per sample
data_agg <- data_filtered %>%
  group_by(project, sample_id, taxon) %>%
  summarise(reads = sum(total_reads), .groups = "drop")

# Create community matrix (samples x taxa)
community_matrix <- data_agg %>%
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) %>%
  column_to_rownames("sample_id")

# Extract project information
project_info <- data_agg %>%
  select(sample_id, project) %>%
  distinct()

# Remove project column from matrix
community_matrix <- community_matrix %>%
  select(-project)

cat("Community matrix:", nrow(community_matrix), "samples x", ncol(community_matrix), "taxa\n")

# Calculate alpha diversity metrics per sample
alpha_diversity <- data.frame(
  sample_id = rownames(community_matrix),
  richness = specnumber(community_matrix),  # Species richness
  shannon = diversity(community_matrix, index = "shannon"),  # Shannon diversity
  simpson = diversity(community_matrix, index = "simpson"),  # Simpson diversity
  evenness = diversity(community_matrix, index = "shannon") / log(specnumber(community_matrix))  # Pielou's evenness
)

# Add project information
alpha_diversity <- alpha_diversity %>%
  left_join(project_info, by = "sample_id")

# Summarize by project
project_summary <- alpha_diversity %>%
  group_by(project) %>%
  summarise(
    n_samples = n(),
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness = sd(richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    sd_simpson = sd(simpson, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE),
    sd_evenness = sd(evenness, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate beta diversity (Bray-Curtis dissimilarity)
if (nrow(community_matrix) > 1) {
  beta_diversity <- vegdist(community_matrix, method = "bray")

  # Add project labels for beta diversity
  sample_projects <- project_info$project[match(rownames(community_matrix), project_info$sample_id)]

  # Calculate average within-project and between-project dissimilarity
  beta_matrix <- as.matrix(beta_diversity)

  beta_summary <- data.frame()
  projects <- unique(sample_projects)

  for (i in 1:length(projects)) {
    for (j in i:length(projects)) {
      proj_i <- projects[i]
      proj_j <- projects[j]

      # Get indices for each project
      idx_i <- which(sample_projects == proj_i)
      idx_j <- which(sample_projects == proj_j)

      if (proj_i == proj_j) {
        # Within-project dissimilarity
        if (length(idx_i) > 1) {
          within_diss <- beta_matrix[idx_i, idx_i]
          within_diss <- within_diss[upper.tri(within_diss)]
          beta_summary <- rbind(beta_summary, data.frame(
            project_1 = proj_i,
            project_2 = proj_j,
            comparison_type = "within",
            mean_dissimilarity = mean(within_diss, na.rm = TRUE),
            sd_dissimilarity = sd(within_diss, na.rm = TRUE),
            n_comparisons = length(within_diss)
          ))
        }
      } else {
        # Between-project dissimilarity
        between_diss <- as.vector(beta_matrix[idx_i, idx_j])
        beta_summary <- rbind(beta_summary, data.frame(
          project_1 = proj_i,
          project_2 = proj_j,
          comparison_type = "between",
          mean_dissimilarity = mean(between_diss, na.rm = TRUE),
          sd_dissimilarity = sd(between_diss, na.rm = TRUE),
          n_comparisons = length(between_diss)
        ))
      }
    }
  }
} else {
  cat("Warning: Not enough samples for beta diversity calculation\n")
  beta_summary <- data.frame()
}

# Combine results
cat("\n=== ALPHA DIVERSITY SUMMARY BY PROJECT ===\n")
print(project_summary)

if (nrow(beta_summary) > 0) {
  cat("\n=== BETA DIVERSITY SUMMARY ===\n")
  print(beta_summary)
}

# Save all results to single CSV
results <- list(
  alpha_per_sample = alpha_diversity,
  alpha_per_project = project_summary,
  beta_diversity = beta_summary
)

# Write with section headers
write_lines("# Alpha Diversity Per Sample", output_file)
write_csv(alpha_diversity, output_file, append = TRUE)

write_lines("\n# Alpha Diversity Per Project", output_file, append = TRUE)
write_csv(project_summary, output_file, append = TRUE)

if (nrow(beta_summary) > 0) {
  write_lines("\n# Beta Diversity (Bray-Curtis)", output_file, append = TRUE)
  write_csv(beta_summary, output_file, append = TRUE)
}

cat("\nDiversity metrics saved to:", output_file, "\n")
cat("------------------------------------------------------------\n")
