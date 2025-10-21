#!/usr/bin/env Rscript

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# -------------------------------
# PARAMETERS FROM ARGS
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 7){
  stop("Usage: Rscript process_motu.R <motu_table> <classification_file1> ... <classification_fileN> <prefixes_comma> <reads_within> <reads_across> <reads_replicates> <output_file>")
}

motu_file <- args[1]                                   # first argument
output_file <- args[length(args)]                      # last argument

# numeric thresholds
reads_replicates <- as.numeric(args[length(args)-1])
reads_across     <- as.numeric(args[length(args)-2])
reads_within     <- as.numeric(args[length(args)-3])

# classification prefixes (comma-separated)
classification_prefixes <- strsplit(args[length(args)-4], ",")[[1]]

# classification files (everything between motu_file and prefixes)
classification_files <- args[2:(length(args)-5)]

if(length(classification_files) != length(classification_prefixes)){
  stop("Number of classification files does not match number of prefixes")
}

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

# -------------------------------
# Summarize reads by sample grouping
# -------------------------------
motu_summary <- motu_long %>%
  group_by(core, depth, sampling_batch, isolation_batch, library, blank_type, sequence_id) %>%
  summarise(
    replicate_reads = list(reads),
    total_reads = sum(reads[reads >= reads_within], na.rm = TRUE),
    n_replicates_present = sum(reads >= reads_within),
    .groups = "drop"
  ) %>%
  mutate(
    replicate_summary = map_chr(replicate_reads, ~ paste(.x, collapse = ";")),
    not_replicated = !(n_replicates_present >= reads_replicates & total_reads >= reads_across)
  ) %>%
  filter(total_reads > 0)

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
# Merge classifications with prefixes
# -------------------------------
classification_tables <- map2(classification_files, classification_prefixes, ~{
  df <- read_csv(.x, show_col_types = FALSE) %>%
    select(id, obitag_bestid, taxid, obitag_rank)
  
  # Extract numeric taxid and taxon name
  numeric_taxid <- as.integer(str_extract(df$taxid, "\\d+"))
  taxon_name <- str_extract(df$taxid, "(?<=\\[)[^\\]]+(?=\\])")
  
  # Add taxon column
  df <- df %>%
    mutate(
      taxon = taxon_name,
      taxid = numeric_taxid
    )
  
  # Explicitly rename columns with the prefix
  df <- df %>% rename(
    !!paste0(.y,"_obitag_bestid") := obitag_bestid,
    !!paste0(.y,"_taxid") := taxid,
    !!paste0(.y,"_taxon") := taxon,
    !!paste0(.y,"_obitag_rank") := obitag_rank
  )
  
  df
})

classification_combined <- reduce(classification_tables, left_join, by = "id")

# Join with MOTU flagged
motu_flagged_classified <- motu_flagged %>%
  left_join(classification_combined, by = c("sequence_id"="id"))

# -------------------------------
# Save final table
# -------------------------------
write_csv(motu_flagged_classified, output_file)