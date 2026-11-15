#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# cluster_taxa.R
#
# Taxon table for one reference database, from the sediment-sample summary
# written by process_motu.R (<project>-bio_sample_summary.tsv: one row per
# bio_sample = core_depth, all isolations and libraries pooled; reads removed by
# a blank are already 0 in the replicates of that isolation/library):
#   1. keep sequences whose assignment in THIS database has
#      bestid >= min_identity (process_motu.R kept a sequence if ANY database
#      passed; here the database's own assignment must pass)
#   2. keep assignments at family level or below (log what is dropped per rank)
#   3. merge the sequences of one taxon within each sample; replicate reads are
#      added replicate by replicate (they are stored in a fixed order with
#      zeros), then totrep / proprep / avgpread / weightrep are recomputed
#
# Usage:
#   Rscript cluster_taxa.R <bio_sample_summary.tsv> <db> <min_identity> <output.csv>
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript cluster_taxa.R <bio_sample_summary.tsv> <db> <min_identity> <output_file>")
}
input_file   <- args[1]
db           <- args[2]
min_identity <- as.numeric(args[3])
output_file  <- args[4]

cat("------------------------------------------------------------\n")
cat("Running cluster_taxa.R\n")
cat("Input file:   ", input_file, "\n")
cat("Database:     ", db, "\n")
cat("Min identity: ", min_identity, "\n")
cat("Output file:  ", output_file, "\n")
cat("------------------------------------------------------------\n")

data <- read_tsv(input_file, show_col_types = FALSE, progress = FALSE,
                 col_types = cols(replicates = col_character(),
                                  replicate_reads = col_character(),
                                  replicate_totals = col_character(),
                                  depth = col_character()))

db_cols <- paste0(c("bestid_", "taxid_", "rank_", "taxon_"), db)
missing_cols <- setdiff(c("bio_sample", "sequence_id", "replicate_reads", "replicate_totals",
                          "use_weighting", db_cols), names(data))
if (length(missing_cols) > 0) {
  stop("Missing columns in ", input_file, ": ", paste(missing_cols, collapse = ", "))
}
data <- data %>%
  rename(obitag_bestid = !!db_cols[1], taxid = !!db_cols[2],
         obitag_rank = !!db_cols[3], taxon = !!db_cols[4])
cat("Loaded", nrow(data), "bio_sample x sequence rows.\n")

# 1. identity in this database
data <- data %>% filter(!is.na(obitag_bestid), obitag_bestid >= min_identity)
cat("Kept", nrow(data), "rows with", db, "identity >=", min_identity, "\n")

# 2. ranks: family level and below. "clade" and "no rank" cannot be placed
# without the lineage and are dropped (counts are logged).
keep_ranks <- c("family", "subfamily", "tribe", "subtribe", "genus", "subgenus",
                "section", "subsection", "series", "species group", "species subgroup",
                "species", "subspecies", "varietas", "forma")
cat("Rows per rank before the rank filter:\n")
print(count(data, obitag_rank, sort = TRUE), n = Inf)
data <- data %>% filter(obitag_rank %in% keep_ranks)
cat("Kept", nrow(data), "rows at family level or below.\n")

# 3. merge sequences per taxon within a bio_sample (replicates are aligned)
to_num <- function(x) as.numeric(str_split(x, ";", simplify = FALSE)[[1]])
group_cols <- c("bio_sample", "lab_sample_ids", "core", "depth", "n_units", "units",
                "replicates", "replicate_totals", "use_weighting", "taxid", "taxon", "obitag_rank")

clustered <- data %>%
  group_by(across(all_of(group_cols))) %>%
  summarise(
    n_sequences    = n(),
    sequence_ids   = paste(sequence_id, collapse = ";"),
    max_bestid     = max(obitag_bestid),
    removed_by     = paste(sort(unique(unlist(str_split(na.omit(removed_by[removed_by != ""]), ";")))), collapse = ";"),
    reads_vec      = list(Reduce(`+`, map(replicate_reads, to_num))),
    .groups = "drop"
  ) %>%
  mutate(
    totals_vec      = map(replicate_totals, to_num),
    replicate_reads = map_chr(reads_vec, ~ paste(.x, collapse = ";")),
    totread         = map_dbl(reads_vec, sum),
    totrep          = map_int(reads_vec, ~ sum(.x > 0)),
    nrep            = lengths(reads_vec),
    proprep         = totrep / nrep,
    props           = map2(reads_vec, totals_vec, ~ ifelse(.y > 0, .x / .y, 0)),
    avgpread        = map_dbl(props, mean),
    sdpread         = map_dbl(props, ~ if (length(.x) > 1) sd(.x) else 0),
    weightrep       = pmap_dbl(list(reads_vec, totals_vec, proprep, use_weighting),
                               function(r, t, p, w) if (w && sum(t) > 0) sum(t[r > 0]) / sum(t) else p)
  ) %>%
  select(-reads_vec, -totals_vec, -props) %>%
  relocate(replicate_reads, .after = replicates) %>%
  arrange(bio_sample, desc(totread))

cat("Clustered into", nrow(clustered), "bio_sample x taxon rows (",
    n_distinct(clustered$taxid), "taxa ).\n")

write_csv(clustered, output_file)
cat("Saved clustered taxa table to:", output_file, "\n")
cat("Done.\n")
cat("------------------------------------------------------------\n")
