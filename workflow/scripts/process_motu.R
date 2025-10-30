#!/usr/bin/env Rscript

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(jsonlite)

# -------------------------------
# PARAMETERS FROM ARGS
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 8){
  stop("Usage: Rscript process_motu.R <motu_table> <classification_file> <demux_stats_json> <db_name> <reads_within> <reads_across> <reads_replicates> <output_file>")
}

motu_file <- args[1]
classification_file <- args[2]
demux_stats_json <- args[3]
db_name <- args[4]
reads_within <- as.numeric(args[5])
reads_across <- as.numeric(args[6])
reads_replicates <- as.numeric(args[7])
output_file <- args[8]

cat("------------------------------------------------------------\n")
cat("Running process_motu.R\n")
cat("MOTU table: ", motu_file, "\n")
cat("Classification file: ", classification_file, "\n")
cat("Demux stats: ", demux_stats_json, "\n")
cat("Database: ", db_name, "\n")
cat("Reads within: ", reads_within, "\n")
cat("Reads across: ", reads_across, "\n")
cat("Reads replicates: ", reads_replicates, "\n")
cat("Output file: ", output_file, "\n")
cat("------------------------------------------------------------\n")

# -------------------------------
# LOAD DEMUX STATS
# -------------------------------
demux_stats <- fromJSON(demux_stats_json)
sample_totals <- demux_stats$samples$sample_stats

# Convert to dataframe with sample name and total reads
sample_read_counts <- tibble(
  sample_id = names(sample_totals),
  total_reads_in_sample = map_dbl(sample_totals, ~.$reads)
)

cat("Loaded demux stats for", nrow(sample_read_counts), "samples.\n")

# -------------------------------
# LOAD MOTU TABLE
# -------------------------------
motu <- read_csv(motu_file, show_col_types = FALSE) %>%
  rename(id = 1)

# -------------------------------
# STEP 0: Parse sample names
# -------------------------------
motu_parsed <- motu %>%
  separate(
    id,
    into = c("core", "depth", "sampling_batch", "isolation_batch", "library", "replicate"),
    sep = "_",
    remove = FALSE,
    fill = "left"
  ) %>%
  mutate(
    blank_type = case_when(
      str_detect(isolation_batch, "^LB") ~ "LB",
      str_detect(isolation_batch, "^PB") ~ "PB",
      str_detect(sampling_batch, "^IB") ~ "IB",
      str_detect(depth, "^SB") ~ "SB",
      TRUE ~ "SAMPLE"
    ),
    sampling_batch = if_else(blank_type == "SAMPLE",
                             paste(core, sampling_batch, sep = "_"),
                             sampling_batch),
    core = if_else(blank_type == "SB", NA_character_, core),
    depth = if_else(blank_type == "SB", NA_character_, depth),
    sampling_batch = if_else(blank_type %in% c("IB","PB","LB"), NA_character_, sampling_batch),
    isolation_batch = if_else(blank_type %in% c("PB","LB"), NA_character_, isolation_batch)
  )

# -------------------------------
# Identify sequence columns
# -------------------------------
seq_cols <- motu_parsed %>%
  select(-id, -core, -depth, -sampling_batch, -isolation_batch, -library, -replicate, -blank_type) %>%
  select(where(is.numeric)) %>%
  colnames()

# -------------------------------
# Long format
# -------------------------------
motu_long <- motu_parsed %>%
  pivot_longer(
    cols = all_of(seq_cols),
    names_to = "sequence_id",
    values_to = "reads"
  )

# Join with sample totals to get proportions
motu_long <- motu_long %>%
  left_join(sample_read_counts, by = c("id" = "sample_id")) %>%
  mutate(
    proportion = reads / total_reads_in_sample,
    proportion = if_else(is.na(proportion) | is.infinite(proportion), 0, proportion)
  )

# -------------------------------
# Summarize reads by sample grouping with weighted proportions
# -------------------------------
motu_summary <- motu_long %>%
  group_by(core, depth, sampling_batch, isolation_batch, library, blank_type, sequence_id) %>%
  summarise(
    replicate_reads = list(reads),
    replicate_proportions = list(proportion),
    replicate_total_reads = list(total_reads_in_sample),
    total_reads = sum(reads[reads >= reads_within], na.rm = TRUE),
    n_replicates_present = sum(reads >= reads_within),
    # CORRECTED: Weighted average proportion (MergeAndFilter style)
    weighted_avg_proportion = sum(proportion * total_reads_in_sample, na.rm = TRUE) / 
                              sum(total_reads_in_sample, na.rm = TRUE),
    # Average proportion across replicates (unweighted)
    mean_proportion = mean(proportion[reads >= reads_within], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    replicate_summary = map_chr(replicate_reads, ~ paste(.x, collapse = ";")),
    proportion_summary = map_chr(replicate_proportions, ~ paste(round(.x, 4), collapse = ";")),
    not_replicated = !(n_replicates_present >= reads_replicates & total_reads >= reads_across)
  ) %>%
  filter(total_reads > 0)

cat("Calculated weighted proportions for", nrow(motu_summary), "sequence-sample combinations.\n")

# -------------------------------
# Compute blank flags
# -------------------------------
lb_flags <- motu_summary %>%
  filter(blank_type == "LB") %>%
  group_by(library, sequence_id) %>%
  summarise(in_LB = any(total_reads > 0), .groups = "drop")

ib_flags <- motu_summary %>%
  filter(blank_type == "IB") %>%
  group_by(isolation_batch, sequence_id) %>%
  summarise(in_IB = any(total_reads > 0), .groups = "drop")

sb_flags <- motu_summary %>%
  filter(blank_type == "SB") %>%
  group_by(sampling_batch, sequence_id) %>%
  summarise(in_SB = any(total_reads > 0), .groups = "drop")

motu_flagged <- motu_summary %>%
  left_join(lb_flags, by = c("library","sequence_id")) %>%
  left_join(ib_flags, by = c("isolation_batch","sequence_id")) %>%
  left_join(sb_flags, by = c("sampling_batch","sequence_id")) %>%
  mutate(across(starts_with("in_"), ~replace_na(., FALSE))) %>%
  mutate(remove = not_replicated | in_LB | in_IB | in_SB)

# -------------------------------
# Load and process classification table
# -------------------------------
classification <- read_csv(classification_file, show_col_types = FALSE) %>%
  select(id, obitag_bestid, taxid, obitag_rank)

# Extract numeric taxid and taxon name
numeric_taxid <- as.integer(str_extract(classification$taxid, "\\d+"))
taxon_name <- str_extract(classification$taxid, "(?<=\\[)[^\\]]+(?=\\])")

# Add taxon column (NO PREFIX)
classification <- classification %>%
  mutate(
    taxon = taxon_name,
    taxid = numeric_taxid
  )

cat("Loaded classification table with", nrow(classification), "sequences.\n")

# Join with MOTU flagged
motu_flagged_classified <- motu_flagged %>%
  left_join(classification, by = c("sequence_id" = "id"))

cat("Final table has", nrow(motu_flagged_classified), "rows.\n")

# -------------------------------
# Save final table
# -------------------------------
write_csv(motu_flagged_classified, output_file)

cat("Saved combined classification table to:", output_file, "\n")
cat("Done.\n")
cat("------------------------------------------------------------\n")