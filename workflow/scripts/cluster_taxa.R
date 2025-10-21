#!/usr/bin/env Rscript

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# -------------------------------
# PARAMETERS FROM ARGS
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript cluster_taxa.R <combined_classification_table.csv> <min_identity> <db_prefix> <output_file>")
}

input_file <- args[1]
min_identity <- as.numeric(args[2])
db_prefix <- args[3]
output_file <- args[4]

cat("------------------------------------------------------------\n")
cat("Running cluster_taxa.R\n")
cat("Input file: ", input_file, "\n")
cat("Min identity: ", min_identity, "\n")
cat("DB prefix: ", db_prefix, "\n")
cat("Output file: ", output_file, "\n")
cat("------------------------------------------------------------\n")

# -------------------------------
# Load combined classification table
# -------------------------------
data <- read_csv(input_file, show_col_types = FALSE)

cat("Loaded table with", nrow(data), "rows and", ncol(data), "columns.\n")

# -------------------------------
# Identify relevant classification columns
# -------------------------------
# Example: PhyloAlps_obitag_bestid, PhyloAlps_taxid, PhyloAlps_taxon, PhyloAlps_obitag_rank
prefix_cols <- names(data)[grepl(paste0("^", db_prefix, "_"), names(data))]

if (length(prefix_cols) == 0) {
  stop(paste0("No columns found for db_prefix '", db_prefix, "'. Check your config or file."))
}

cat("Using columns with prefix:", db_prefix, "\n")

# Select only relevant classification subset
id_col <- paste0(db_prefix, "_obitag_bestid")

db_data <- data %>%
  select(core, depth, sampling_batch, isolation_batch, library, blank_type, total_reads, remove,
         all_of(prefix_cols)
  )

# -------------------------------
# Extract and clean taxonomic info
# -------------------------------
# Typically, obitag_bestid contains the % identity and possibly other fields
id_col <- paste0(db_prefix, "_obitag_bestid")
rank_col <- paste0(db_prefix, "_obitag_rank")
taxid_col <- paste0(db_prefix, "_taxid")
taxon_col <- paste0(db_prefix, "_taxon")

if (!all(c(id_col, rank_col, taxid_col, taxon_col) %in% names(db_data))) {
  stop(paste0("Missing expected columns for prefix '", db_prefix, "'."))
}

db_data<- db_data %>%
  rename_with(~ str_remove(., paste0("^", db_prefix, "_")))

# -------------------------------
# Filter by min_identity and filter flag
# -------------------------------
filtered_data <- db_data %>%
  filter(remove==FALSE)
cat("Kept", nrow(filtered_data), "rows after removing flagged sequences.\n")


filtered_data <- filtered_data %>%
  filter(obitag_bestid >= min_identity)
cat("Kept", nrow(filtered_data), "rows after applying min_identity filter.\n")

# drop unnecesary columns
filtered_data <- filtered_data %>%
  select(core, depth, sampling_batch, isolation_batch, library, blank_type, total_reads,
         obitag_bestid, taxid, obitag_rank, taxon)

# -------------------------------
# Cluster taxa by identical taxid or taxon
# -------------------------------
# Example: group by taxid (or taxon if missing), summarise sequences
clustered <- filtered_data %>%
  group_by(core, depth, sampling_batch, isolation_batch, library, blank_type, taxid, taxon, obitag_rank) %>%
  summarise(
    total_reads = sum(total_reads),
    .groups = "drop"
  ) 

cat("Clustered into", nrow(clustered), "unique taxa.\n")

# -------------------------------
# Filter by taxonimic level, minimum to family
# ------------------------------- 

clustered <- clustered %>% 
  filter(obitag_rank %in% c("species", "subgenus", "section", "genus", "family", "subfamily", "tribe"))

cat("Filtered to", nrow(clustered), "rows at or above family level.\n")

# -------------------------------
# Save clustered output
# -------------------------------
write_csv(clustered, output_file)

cat("Saved clustered taxa table to:", output_file, "\n")
cat("Done.\n")
cat("------------------------------------------------------------\n")