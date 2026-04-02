#!/usr/bin/env Rscript

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Merge replicate-level read summaries across sequences belonging to one taxon.
merge_replicate_summaries <- function(replicate_summaries) {
  parsed <- lapply(replicate_summaries, function(x) {
    if (is.na(x) || x == "") {
      return(numeric(0))
    }
    as.numeric(str_split(x, ";", simplify = FALSE)[[1]])
  })
  
  max_len <- max(lengths(parsed), na.rm = TRUE)
  if (!is.finite(max_len) || max_len == 0) {
    return(list(
      n_replicates_present = 0L,
      replicate_summary = NA_character_
    ))
  }
  
  padded <- lapply(parsed, function(v) {
    length(v) <- max_len
    v[is.na(v)] <- 0
    v
  })
  
  merged_counts <- Reduce(`+`, padded)
  list(
    n_replicates_present = as.integer(sum(merged_counts > 0)),
    replicate_summary = paste(as.integer(merged_counts), collapse = ";")
  )
}

# -------------------------------
# PARAMETERS FROM ARGS
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript cluster_taxa.R <combined_classification_table.csv> <min_identity> <output_file>")
}

input_file <- args[1]
min_identity <- as.numeric(args[2])
output_file <- args[3]

cat("------------------------------------------------------------\n")
cat("Running cluster_taxa.R\n")
cat("Input file: ", input_file, "\n")
cat("Min identity: ", min_identity, "\n")
cat("Output file: ", output_file, "\n")
cat("------------------------------------------------------------\n")

# -------------------------------
# Load combined classification table
# -------------------------------
data <- read_csv(input_file, show_col_types = FALSE)

cat("Loaded table with", nrow(data), "rows and", ncol(data), "columns.\n")

# -------------------------------
# Check required columns exist
# -------------------------------
required_cols <- c("core", "depth", "sampling_batch", "isolation_batch", "library", 
                   "blank_type", "total_reads", "remove", 
                   "obitag_bestid", "taxid", "obitag_rank", "taxon", "replicate_summary", "replicate_summary_filtered")

missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
}

# -------------------------------
# Filter by min_identity and filter flag
# -------------------------------
filtered_data <- data %>%
  filter(remove == FALSE)
cat("Kept", nrow(filtered_data), "rows after removing flagged sequences.\n")

filtered_data <- filtered_data %>%
  filter(obitag_bestid >= min_identity)
cat("Kept", nrow(filtered_data), "rows after applying min_identity filter.\n")

# Select only necessary columns
filtered_data <- filtered_data %>%
  select(core, depth, sampling_batch, isolation_batch, library, blank_type, 
         total_reads, replicate_summary, replicate_summary_filtered, obitag_bestid, taxid, obitag_rank, taxon)

# -------------------------------
# Cluster taxa by identical taxid or taxon
# -------------------------------
clustered <- filtered_data %>%
  group_by(core, depth, sampling_batch, isolation_batch, library, blank_type, 
           taxid, taxon, obitag_rank) %>%
  group_modify(~{
    merged_rep <- merge_replicate_summaries(.x$replicate_summary)
    merged_rep_filtered <- merge_replicate_summaries(.x$replicate_summary_filtered)
    tibble(
      total_reads = sum(.x$total_reads, na.rm = TRUE),
      n_replicates_present = merged_rep$n_replicates_present,
      replicate_summary = merged_rep$replicate_summary,
      n_replicates_present_filtered = merged_rep_filtered$n_replicates_present,
      replicate_summary_filtered = merged_rep_filtered$replicate_summary
    )
  }) %>%
  ungroup()

cat("Clustered into", nrow(clustered), "unique taxa.\n")

# -------------------------------
# Filter by taxonomic level, minimum to family
# ------------------------------- 
clustered <- clustered %>% 
  filter(obitag_rank %in% c("species", "subgenus", "section", "genus", 
                            "family", "subfamily", "tribe"))

cat("Filtered to", nrow(clustered), "rows at or above family level.\n")

# -------------------------------
# Save clustered output
# -------------------------------
write_csv(clustered, output_file)

cat("Saved clustered taxa table to:", output_file, "\n")
cat("Done.\n")
cat("------------------------------------------------------------\n")