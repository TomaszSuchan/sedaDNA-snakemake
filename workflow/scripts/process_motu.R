#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# process_motu.R
#
# Merge the OBITools MOTU table with the classifications from all reference
# databases and filter it, following the design of MergeAndFilter
# (Lammers, https://github.com/Y-Lammers/MergeAndFilter):
#
#   1. per observation (one sequence in one PCR replicate):
#        reads < reads_within  -> set to 0
#   2. per sequence, over the whole dataset (all SAMPLE PCR replicates):
#        total filtered reads   >= reads_across
#        number of observations >= reads_replicates   (replicates with reads > 0)
#   3. per sequence: kept if the identity threshold is met in ANY database
#   4. per sample: sequences found in a matching blank are flagged and removed
#        LB, PB -> same library;  IB -> same isolation batch;
#        SB     -> same core + sampling batch (e.g. CSG-8_D1)
#   5. per sample summary (MergeAndFilter "_summary"): totread, avgpread,
#      sdpread, totrep, proprep, weightrep, with the replicate reads kept in a
#      fixed replicate order (zeros included) so that later merging is exact.
#
# One run per project: the filtering does not depend on the database, so all
# databases get identical filtering; classification columns are kept per
# database (bestid_<db>, taxid_<db>, rank_<db>, taxon_<db>).
#
# Usage:
#   Rscript process_motu.R <motu_table.csv> <demux_stats_combined.json> \
#       <reads_within> <reads_across> <reads_replicates> <min_identity> \
#       <output_prefix> <db1>=<classification_table1.csv> [<db2>=<...> ...]
#
# Outputs (<output_prefix> = results/{project}/tables/{project}):
#   <prefix>-replicates.tsv.gz    every observation: pcr_id, sequence_id, reads
#                                 (raw obiclean weights, blanks included)
#   <prefix>-pcrs.tsv             one row per PCR replicate: parsed name, raw
#                                 reads and the share of reads removed per step
#   <prefix>-samples.tsv          one row per sample (all its replicates): read
#                                 totals and replication of the top-10 sequences
#   <prefix>-sequences.tsv        one row per sequence: dataset totals, filter
#                                 results and the assignment from each database
#   <prefix>-sample_summary.tsv   one row per sample (isolation x library) x retained sequence
#   <prefix>-bio_sample_summary.tsv  one row per sediment sample (core_depth, all its
#                                 isolations and libraries pooled) x sequence
#   <prefix>-blank_report.tsv     sequences found in blanks and the samples they
#                                 removed
#   <prefix>-filter_params.json   the parameters of this run
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste("Usage: Rscript process_motu.R <motu_table> <demux_stats_json>",
             "<reads_within> <reads_across> <reads_replicates> <min_identity>",
             "<output_prefix> <db>=<classification_table> [...]"))
}

motu_file        <- args[1]
demux_stats_json <- args[2]
reads_within     <- as.numeric(args[3])
reads_across     <- as.numeric(args[4])
reads_replicates <- as.numeric(args[5])
min_identity     <- as.numeric(args[6])
out_prefix       <- args[7]
db_args          <- args[-(1:7)]

if (any(!str_detect(db_args, "^[^=]+=.+$"))) {
  stop("Classification tables must be given as <db>=<path>, got: ",
       paste(db_args, collapse = " "))
}
db_files <- setNames(str_replace(db_args, "^[^=]+=", ""),
                     str_extract(db_args, "^[^=]+"))

out <- function(suffix) paste0(out_prefix, "-", suffix)

cat("------------------------------------------------------------\n")
cat("Running process_motu.R (MergeAndFilter design)\n")
cat("MOTU table:       ", motu_file, "\n")
cat("Demux stats:      ", demux_stats_json, "\n")
cat("Databases:        ", paste(names(db_files), collapse = ", "), "\n")
cat("reads_within:     ", reads_within, " (min reads per observation)\n")
cat("reads_across:     ", reads_across, " (min filtered reads per sequence, whole dataset)\n")
cat("reads_replicates: ", reads_replicates, " (min observations per sequence, whole dataset)\n")
cat("min_identity:     ", min_identity, " (in at least one database)\n")
cat("Output prefix:    ", out_prefix, "\n")
cat("------------------------------------------------------------\n")

# ------------------------------------------------------------------------------
# PCR names. Formats (same as rule validate_samples in demultiplex.smk):
#   LB_<library>_<rep>                              library blank
#   PB_<library>_<rep>                              PCR blank
#   IB_<isolation>_<library>_<rep>                  isolation blank
#   <core>_SB_<Dn>_<isolation>_<library>_<rep>      sampling blank
#   <core>_<depth>_<Dn>_<isolation>_<library>_<rep> sample
# ------------------------------------------------------------------------------
parse_pcr_ids <- function(id) {
  m_lp <- str_match(id, "^(LB|PB)_([^_]+)_(\\d+)$")
  m_ib <- str_match(id, "^IB_([^_]+)_([^_]+)_(\\d+)$")
  m_sb <- str_match(id, "^([^_]+)_SB_([^_]+)_([^_]+)_([^_]+)_(\\d+)$")
  m_s  <- str_match(id, "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(\\d+)$")
  type <- case_when(!is.na(m_lp[, 1]) ~ m_lp[, 2],
                    !is.na(m_ib[, 1]) ~ "IB",
                    !is.na(m_sb[, 1]) ~ "SB",
                    !is.na(m_s[, 1])  ~ "SAMPLE",
                    TRUE ~ NA_character_)
  pick <- function(lp, ib, sb, s) {
    case_when(type %in% c("LB", "PB") ~ lp, type == "IB" ~ ib,
              type == "SB" ~ sb, type == "SAMPLE" ~ s, TRUE ~ NA_character_)
  }
  na <- rep(NA_character_, length(id))
  tibble(
    pcr_id          = id,
    blank_type      = type,
    core            = pick(na, na, m_sb[, 2], m_s[, 2]),
    depth           = pick(na, na, na, m_s[, 3]),
    sampling_batch  = pick(na, na, paste(m_sb[, 2], m_sb[, 3], sep = "_"),
                           paste(m_s[, 2], m_s[, 4], sep = "_")),
    isolation_batch = pick(na, m_ib[, 2], m_sb[, 4], m_s[, 5]),
    library         = pick(m_lp[, 3], m_ib[, 3], m_sb[, 5], m_s[, 6]),
    replicate       = as.integer(pick(m_lp[, 4], m_ib[, 4], m_sb[, 6], m_s[, 7]))
  ) %>%
    mutate(
      sample_id  = str_replace(pcr_id, "_\\d+$", ""),
      # the sediment sample: one depth of one core (core_depth), samples only
      bio_sample = if_else(blank_type == "SAMPLE", paste(core, depth, sep = "_"), NA_character_),
      # the sample as named in the lab database (core_depth_Dn)
      lab_sample_id = if_else(blank_type == "SAMPLE",
                              paste(core, depth, str_replace(sampling_batch, "^[^_]+_", ""), sep = "_"),
                              NA_character_)
    )
}

# ------------------------------------------------------------------------------
# 1. Demultiplexing stats = the full list of PCR replicates and their raw reads
# ------------------------------------------------------------------------------
demux_stats <- fromJSON(demux_stats_json, simplifyVector = FALSE)
raw_counts <- map_dfr(names(demux_stats), function(lib_name) {
  st <- demux_stats[[lib_name]]$samples$sample_stats
  tibble(pcr_id = names(st),
         raw_reads = map_dbl(st, ~ as.numeric(.x$reads)),
         library_file = lib_name)
})
if (anyDuplicated(raw_counts$pcr_id)) {
  dup <- unique(raw_counts$pcr_id[duplicated(raw_counts$pcr_id)])
  stop("PCR names occur in more than one library: ", paste(head(dup, 10), collapse = ", "))
}
cat("Demux stats:", nrow(raw_counts), "PCR replicates in", length(demux_stats), "libraries.\n")

# ------------------------------------------------------------------------------
# 2. MOTU table (rows = PCR replicates, columns = sequences) -> long
# ------------------------------------------------------------------------------
motu <- read_csv(motu_file, show_col_types = FALSE, progress = FALSE)
names(motu)[1] <- "pcr_id"
seq_cols <- names(motu)[-1][map_lgl(motu[-1], is.numeric)]
obs <- motu %>%
  select(pcr_id, all_of(seq_cols)) %>%
  pivot_longer(-pcr_id, names_to = "sequence_id", values_to = "reads") %>%
  filter(!is.na(reads), reads > 0)
rm(motu); invisible(gc())
cat("MOTU table:", n_distinct(obs$pcr_id), "PCRs,", length(seq_cols), "sequences,",
    nrow(obs), "observations.\n")

not_in_stats <- setdiff(unique(obs$pcr_id), raw_counts$pcr_id)
if (length(not_in_stats) > 0) {
  cat("WARNING:", length(not_in_stats), "PCRs of the MOTU table are not in the demux stats",
      "(raw reads unknown):", paste(head(not_in_stats, 10), collapse = ", "), "\n")
}

pcrs <- parse_pcr_ids(union(raw_counts$pcr_id, unique(obs$pcr_id))) %>%
  left_join(raw_counts, by = "pcr_id")
bad <- pcrs %>% filter(is.na(blank_type))
if (nrow(bad) > 0) {
  writeLines(bad$pcr_id, out("unparsed_pcr_names.txt"))
  stop(nrow(bad), " PCR names do not follow the naming scheme (listed in ",
       out("unparsed_pcr_names.txt"), "), e.g. ", paste(head(bad$pcr_id, 5), collapse = ", "))
}
cat("PCR replicates by type:\n"); print(count(pcrs, blank_type))
nrep_tab <- pcrs %>% count(sample_id, name = "n_reps")
cat("Replicates per sample:\n"); print(count(nrep_tab, n_reps, name = "n_samples"))

write_tsv(obs %>% arrange(pcr_id, sequence_id), out("replicates.tsv.gz"))

obs <- obs %>%
  left_join(pcrs %>% select(pcr_id, sample_id, replicate, blank_type, library,
                            isolation_batch, sampling_batch), by = "pcr_id") %>%
  mutate(reads_f = if_else(reads >= reads_within, reads, 0))   # step 1

# ------------------------------------------------------------------------------
# 3. Dataset-level filters per sequence (samples only; blanks do not give support)
# ------------------------------------------------------------------------------
seqs <- obs %>%
  group_by(sequence_id) %>%
  summarise(
    reads_samples      = sum(reads[blank_type == "SAMPLE"]),
    reads_samples_filt = sum(reads_f[blank_type == "SAMPLE"]),
    obs_samples        = sum(reads_f[blank_type == "SAMPLE"] > 0),
    reads_blanks       = sum(reads[blank_type != "SAMPLE"]),
    .groups = "drop"
  ) %>%
  mutate(pass_reads_across     = reads_samples_filt >= reads_across,
         pass_reads_replicates = obs_samples >= reads_replicates,
         pass_technical        = pass_reads_across & pass_reads_replicates)

# ------------------------------------------------------------------------------
# 4. Classifications of all databases; identity filter = any database
# ------------------------------------------------------------------------------
read_classification <- function(db, path) {
  cl <- read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    select(id, obitag_bestid, taxid, obitag_rank)
  tibble(
    sequence_id = cl$id,
    bestid = as.numeric(cl$obitag_bestid),
    taxid  = as.integer(str_extract(cl$taxid, "\\d+")),
    rank   = cl$obitag_rank,
    taxon  = str_extract(cl$taxid, "(?<=\\[)[^\\]]+(?=\\])")
  ) %>% rename_with(~ paste0(.x, "_", db), -sequence_id)
}
for (db in names(db_files)) {
  cl <- read_classification(db, db_files[[db]])
  cat("Classification", db, ":", nrow(cl), "sequences.\n")
  miss <- setdiff(seqs$sequence_id, cl$sequence_id)
  if (length(miss) > 0) {
    cat("WARNING:", length(miss), "sequences of the MOTU table have no", db, "classification.\n")
  }
  seqs <- left_join(seqs, cl, by = "sequence_id")
}
bestid_cols <- paste0("bestid_", names(db_files))
bestid_mat <- as.matrix(seqs[bestid_cols])
best_idx <- max.col(replace(bestid_mat, is.na(bestid_mat), -1), ties.method = "first")
taxon_mat <- as.matrix(seqs[paste0("taxon_", names(db_files))])
seqs <- seqs %>%
  mutate(
    best_bestid   = suppressWarnings(apply(bestid_mat, 1, max, na.rm = TRUE)),
    best_bestid   = if_else(is.finite(best_bestid), best_bestid, NA_real_),
    best_db       = names(db_files)[best_idx],
    best_taxon    = taxon_mat[cbind(seq_along(best_idx), best_idx)],
    pass_identity = !is.na(best_bestid) & best_bestid >= min_identity,
    retained      = pass_technical & pass_identity
  )
retained_ids <- seqs$sequence_id[seqs$retained]
cat("Sequences: ", nrow(seqs), " | pass reads_across: ", sum(seqs$pass_reads_across),
    " | pass reads_replicates: ", sum(seqs$pass_reads_replicates),
    " | pass both: ", sum(seqs$pass_technical),
    " | + identity in any db: ", length(retained_ids), "\n", sep = "")

# ------------------------------------------------------------------------------
# 5. Blank flags: presence (reads >= reads_within) in a matching blank
# ------------------------------------------------------------------------------
blank_keys <- c(LB = "library", PB = "library", IB = "isolation_batch", SB = "sampling_batch")
blank_presence <- function(type, key) {
  obs %>% filter(blank_type == type, reads_f > 0) %>%
    distinct(key = .data[[key]], sequence_id) %>%
    mutate(!!paste0("in_", type) := TRUE)
}
flags <- imap(blank_keys, ~ blank_presence(.y, .x))

sample_meta <- pcrs %>%
  group_by(sample_id) %>%
  summarise(across(c(blank_type, core, depth, bio_sample, lab_sample_id, sampling_batch,
                     isolation_batch, library), first),
            n_reps = n(),
            raw_reads_sample = sum(raw_reads, na.rm = TRUE),
            .groups = "drop")

present <- obs %>%
  filter(sequence_id %in% retained_ids, reads_f > 0) %>%
  distinct(sample_id, sequence_id) %>%
  left_join(sample_meta %>% select(sample_id, blank_type, library, isolation_batch, sampling_batch),
            by = "sample_id")
for (bt in names(blank_keys)) {
  present <- present %>%
    left_join(flags[[bt]], by = setNames(c("key", "sequence_id"), c(blank_keys[[bt]], "sequence_id")))
}
present <- present %>%
  mutate(across(c(in_LB, in_PB, in_IB, in_SB), ~ replace_na(.x, FALSE) & blank_type == "SAMPLE"),
         n_blank_flags = in_LB + in_PB + in_IB + in_SB,
         remove_blank  = n_blank_flags > 0) %>%
  select(sample_id, sequence_id, in_LB, in_PB, in_IB, in_SB, n_blank_flags, remove_blank)

sample_keys <- pcrs %>% filter(blank_type == "SAMPLE")
for (bt in names(blank_keys)) {
  key <- blank_keys[[bt]]
  bl_keys <- unique(pcrs[[key]][pcrs$blank_type == bt])
  orphan <- setdiff(bl_keys, sample_keys[[key]])
  n_fl <- sum(present[[paste0("in_", bt)]])
  cat(sprintf("Blank %s: %d PCRs, %d sample x sequence rows flagged\n",
              bt, sum(pcrs$blank_type == bt), n_fl))
  if (length(orphan) > 0) {
    # a blank whose key (library / isolation batch / core_Dn) matches no sample
    # cannot flag anything: almost always a naming problem
    cat("  WARNING:", bt, "blanks with no matching samples (", key, "):",
        paste(orphan, collapse = ", "), "\n")
  }
}

# ------------------------------------------------------------------------------
# 6. Per-sample summary with replicates in fixed order (zeros included)
# ------------------------------------------------------------------------------
reps <- pcrs %>% select(sample_id, pcr_id, replicate)
grid <- present %>%
  select(sample_id, sequence_id, remove_blank) %>%
  inner_join(reps, by = "sample_id", relationship = "many-to-many") %>%
  left_join(obs %>% select(pcr_id, sequence_id, reads_f), by = c("pcr_id", "sequence_id")) %>%
  mutate(reads_f = coalesce(reads_f, 0))

# retained reads per replicate (denominator of the proportions, as in MergeAndFilter)
rep_totals <- grid %>% filter(!remove_blank) %>%
  group_by(pcr_id) %>% summarise(rep_total = sum(reads_f), .groups = "drop")
rep_weights <- reps %>%
  left_join(rep_totals, by = "pcr_id") %>%
  mutate(rep_total = coalesce(rep_total, 0)) %>%
  group_by(sample_id) %>%
  mutate(rep_weight = if (sum(rep_total) > 0) rep_total / sum(rep_total) else 0) %>%
  ungroup() %>%
  select(pcr_id, rep_total, rep_weight)
grid <- grid %>%
  left_join(rep_weights, by = "pcr_id") %>%
  mutate(prop = if_else(rep_total > 0, reads_f / rep_total, 0)) %>%
  arrange(sample_id, sequence_id, replicate)

summary_tbl <- grid %>%
  group_by(sample_id, sequence_id) %>%
  summarise(
    replicates       = paste(replicate, collapse = ";"),
    replicate_reads  = paste(reads_f, collapse = ";"),
    replicate_totals = paste(rep_total, collapse = ";"),
    totread          = sum(reads_f),
    totrep           = sum(reads_f > 0),
    nrep             = n(),
    avgpread         = mean(prop),
    sdpread          = if (n() > 1) sd(prop) else 0,
    weightrep_raw    = sum(rep_weight[reads_f > 0]),
    .groups = "drop"
  ) %>%
  mutate(proprep = totrep / nrep)

# MergeAndFilter: use the weighted proportion only when the sample is well
# replicated (mean replication of its barcodes >= 0.33) and has >= 10 barcodes
weighting <- summary_tbl %>%
  inner_join(present %>% filter(!remove_blank) %>% select(sample_id, sequence_id),
             by = c("sample_id", "sequence_id")) %>%
  group_by(sample_id) %>%
  summarise(mean_proprep = mean(proprep), n_barcodes = n(), .groups = "drop") %>%
  mutate(use_weighting = mean_proprep >= 0.33 & n_barcodes >= 10)

summary_tbl <- summary_tbl %>%
  left_join(weighting %>% select(sample_id, use_weighting), by = "sample_id") %>%
  mutate(use_weighting = coalesce(use_weighting, FALSE),
         weightrep = if_else(use_weighting, weightrep_raw, proprep)) %>%
  select(-weightrep_raw) %>%
  left_join(present, by = c("sample_id", "sequence_id")) %>%
  left_join(sample_meta, by = "sample_id") %>%
  mutate(prop_of_raw = if_else(raw_reads_sample > 0, totread / raw_reads_sample, NA_real_)) %>%
  left_join(seqs %>% select(sequence_id, best_bestid, best_db, best_taxon,
                            starts_with("bestid_"), starts_with("taxid_"),
                            starts_with("rank_"), starts_with("taxon_")),
            by = "sequence_id") %>%
  select(sample_id, bio_sample, lab_sample_id, blank_type, core, depth, sampling_batch, isolation_batch,
         library, sequence_id, n_reps, replicates, replicate_reads, replicate_totals,
         totread, avgpread, sdpread, totrep, nrep, proprep, weightrep, use_weighting,
         prop_of_raw, in_LB, in_PB, in_IB, in_SB, n_blank_flags, remove_blank,
         everything(), -raw_reads_sample) %>%
  arrange(sample_id, desc(totread))

write_tsv(summary_tbl, out("sample_summary.tsv"))
cat("Sample summary:", nrow(summary_tbl), "sample x sequence rows,",
    sum(!summary_tbl$remove_blank), "not removed by blanks.\n")

# ------------------------------------------------------------------------------
# 6b. Pool to the sediment sample (bio_sample = core_depth)
#
# All PCR replicates of all isolations / libraries of one bio_sample become the
# replicates of one sample. Blank removal stays at the level of each blank:
# a sequence removed by a blank in one isolation/library counts as 0 reads in
# the replicates of that isolation/library and keeps its reads in the others.
# ------------------------------------------------------------------------------
pool_reps <- pcrs %>%
  filter(blank_type == "SAMPLE") %>%
  select(bio_sample, sample_id, pcr_id, replicate) %>%
  left_join(rep_weights %>% select(pcr_id, rep_total), by = "pcr_id") %>%
  mutate(rep_total = coalesce(rep_total, 0),
         # replicate label inside the bio_sample: <Dn>_<isolation>_<library>_<rep>
         rep_label = substring(pcr_id, nchar(bio_sample) + 2)) %>%
  group_by(bio_sample) %>%
  mutate(rep_weight = if (sum(rep_total) > 0) rep_total / sum(rep_total) else 0) %>%
  ungroup() %>%
  arrange(bio_sample, sample_id, replicate)

kept_present <- present %>%
  filter(!remove_blank) %>%
  inner_join(sample_meta %>% filter(blank_type == "SAMPLE") %>% select(sample_id, bio_sample),
             by = "sample_id")

pool_grid <- kept_present %>%
  distinct(bio_sample, sequence_id) %>%
  inner_join(pool_reps, by = "bio_sample", relationship = "many-to-many") %>%
  left_join(kept_present %>% transmute(sample_id, sequence_id, kept = TRUE),
            by = c("sample_id", "sequence_id")) %>%
  left_join(obs %>% select(pcr_id, sequence_id, reads_f), by = c("pcr_id", "sequence_id")) %>%
  mutate(reads_f = if_else(coalesce(kept, FALSE), coalesce(reads_f, 0), 0),
         prop = if_else(rep_total > 0, reads_f / rep_total, 0)) %>%
  arrange(bio_sample, sequence_id, sample_id, replicate)

pool_units <- sample_meta %>%
  filter(blank_type == "SAMPLE") %>%
  group_by(bio_sample) %>%
  summarise(lab_sample_ids = paste(sort(unique(lab_sample_id)), collapse = ";"),
            core = first(core), depth = first(depth),
            n_units = n(),                                 # isolation x library combinations
            units = paste(sample_id, collapse = ";"),
            raw_reads_bio = sum(raw_reads_sample),
            .groups = "drop")

bio_summary <- pool_grid %>%
  group_by(bio_sample, sequence_id) %>%
  summarise(
    replicates       = paste(rep_label, collapse = ";"),
    replicate_reads  = paste(reads_f, collapse = ";"),
    replicate_totals = paste(rep_total, collapse = ";"),
    totread          = sum(reads_f),
    totrep           = sum(reads_f > 0),
    nrep             = n(),
    avgpread         = mean(prop),
    sdpread          = if (n() > 1) sd(prop) else 0,
    weightrep_raw    = sum(rep_weight[reads_f > 0]),
    .groups = "drop"
  ) %>%
  mutate(proprep = totrep / nrep)

bio_weighting <- bio_summary %>%
  group_by(bio_sample) %>%
  summarise(mean_proprep = mean(proprep), n_barcodes = n(), .groups = "drop") %>%
  mutate(use_weighting = mean_proprep >= 0.33 & n_barcodes >= 10)

bio_flags <- present %>%
  inner_join(sample_meta %>% filter(blank_type == "SAMPLE") %>% select(sample_id, bio_sample),
             by = "sample_id") %>%
  group_by(bio_sample, sequence_id) %>%
  summarise(units_present = n(), units_removed_by_blank = sum(remove_blank),
            removed_by = paste(sort(unique(c(if (any(in_LB)) "LB", if (any(in_PB)) "PB",
                                             if (any(in_IB)) "IB", if (any(in_SB)) "SB"))),
                               collapse = ";"),
            .groups = "drop")

bio_summary <- bio_summary %>%
  left_join(bio_weighting %>% select(bio_sample, use_weighting), by = "bio_sample") %>%
  mutate(weightrep = if_else(use_weighting, weightrep_raw, proprep)) %>%
  select(-weightrep_raw) %>%
  left_join(pool_units, by = "bio_sample") %>%
  left_join(bio_flags, by = c("bio_sample", "sequence_id")) %>%
  mutate(prop_of_raw = if_else(raw_reads_bio > 0, totread / raw_reads_bio, NA_real_)) %>%
  left_join(seqs %>% select(sequence_id, best_bestid, best_db, best_taxon,
                            starts_with("bestid_"), starts_with("taxid_"),
                            starts_with("rank_"), starts_with("taxon_")),
            by = "sequence_id") %>%
  select(bio_sample, lab_sample_ids, core, depth, n_units, units, sequence_id,
         replicates, replicate_reads, replicate_totals, totread, avgpread, sdpread,
         totrep, nrep, proprep, weightrep, use_weighting, prop_of_raw,
         units_present, units_removed_by_blank, removed_by, everything(), -raw_reads_bio) %>%
  arrange(bio_sample, desc(totread))

write_tsv(bio_summary, out("bio_sample_summary.tsv"))
cat("Bio-sample summary:", n_distinct(bio_summary$bio_sample), "sediment samples,",
    nrow(bio_summary), "bio_sample x sequence rows (pooled over",
    sum(pool_units$n_units), "isolation/library units).\n")

# ------------------------------------------------------------------------------
# 7. Sequences table
# ------------------------------------------------------------------------------
write_tsv(seqs %>% arrange(desc(reads_samples)), out("sequences.tsv"))

# ------------------------------------------------------------------------------
# 8. Per-PCR read accounting (MergeAndFilter "_samplestat")
# ------------------------------------------------------------------------------
kept_rows <- present %>% filter(!remove_blank) %>% select(sample_id, sequence_id)
pcr_acc <- obs %>%
  left_join(seqs %>% select(sequence_id, pass_technical, retained), by = "sequence_id") %>%
  left_join(kept_rows %>% mutate(kept = TRUE), by = c("sample_id", "sequence_id")) %>%
  mutate(kept = coalesce(kept, FALSE)) %>%
  group_by(pcr_id) %>%
  summarise(reads_motu        = sum(reads),
            reads_technical   = sum(reads_f[pass_technical]),
            reads_identity    = sum(reads_f[retained]),
            reads_final       = sum(reads_f[retained & kept]),
            n_sequences_final = sum(retained & kept & reads_f > 0),
            .groups = "drop")
pcr_tbl <- pcrs %>%
  left_join(pcr_acc, by = "pcr_id") %>%
  mutate(across(c(reads_motu, reads_technical, reads_identity, reads_final, n_sequences_final),
                ~ coalesce(.x, 0)),
         prop_removed_technical = if_else(raw_reads > 0, (raw_reads - reads_technical) / raw_reads, NA_real_),
         prop_removed_identity  = if_else(raw_reads > 0, (reads_technical - reads_identity) / raw_reads, NA_real_),
         prop_removed_blank     = if_else(raw_reads > 0, (reads_identity - reads_final) / raw_reads, NA_real_),
         prop_retained          = if_else(raw_reads > 0, reads_final / raw_reads, NA_real_)) %>%
  arrange(sample_id, replicate)
write_tsv(pcr_tbl, out("pcrs.tsv"))

# per sample: replication of the 10 most abundant sequences before and after
# filtering (MergeAndFilter avg_rep / avg_filt_rep / overlap)
top10 <- function(d) {
  d %>% group_by(sample_id, sequence_id) %>%
    summarise(tot = sum(reads), nrep_present = sum(reads > 0), .groups = "drop") %>%
    filter(tot > 0) %>%
    group_by(sample_id) %>% slice_max(tot, n = 10, with_ties = FALSE) %>%
    summarise(ids = list(sequence_id), reps = list(nrep_present), .groups = "drop")
}
raw_top  <- top10(obs %>% select(sample_id, sequence_id, reads))
filt_top <- top10(obs %>% semi_join(kept_rows, by = c("sample_id", "sequence_id")) %>%
                    transmute(sample_id, sequence_id, reads = reads_f))
samples_tbl <- pcr_tbl %>%
  group_by(sample_id) %>%
  summarise(across(c(blank_type, bio_sample, lab_sample_id, core, depth, sampling_batch, isolation_batch, library), first),
            n_reps        = n(),
            raw_reads     = sum(raw_reads, na.rm = TRUE),
            reads_final   = sum(reads_final),
            prop_retained = if_else(raw_reads > 0, reads_final / raw_reads, NA_real_),
            .groups = "drop") %>%
  left_join(kept_rows %>% count(sample_id, name = "n_sequences_final"), by = "sample_id") %>%
  mutate(n_sequences_final = coalesce(n_sequences_final, 0L)) %>%
  left_join(raw_top %>% transmute(sample_id, raw_ids = ids, raw_reps = reps), by = "sample_id") %>%
  left_join(filt_top %>% transmute(sample_id, filt_ids = ids, filt_reps = reps), by = "sample_id") %>%
  mutate(avg_rep      = map2_dbl(raw_reps, n_reps, ~ if (is.null(.x)) NA_real_ else mean(.x) / .y),
         sd_rep       = map2_dbl(raw_reps, n_reps, ~ if (length(.x) < 2) NA_real_ else sd(.x / .y)),
         avg_filt_rep = map2_dbl(filt_reps, n_reps, ~ if (is.null(.x)) NA_real_ else mean(.x) / .y),
         sd_filt_rep  = map2_dbl(filt_reps, n_reps, ~ if (length(.x) < 2) NA_real_ else sd(.x / .y)),
         overlap      = map2_int(raw_ids, filt_ids, ~ length(intersect(.x, .y)))) %>%
  select(-raw_ids, -raw_reps, -filt_ids, -filt_reps)
write_tsv(samples_tbl, out("samples.tsv"))

# ------------------------------------------------------------------------------
# 9. Blank report: what each blank removed
# ------------------------------------------------------------------------------
blank_reads <- obs %>% filter(blank_type != "SAMPLE", reads_f > 0) %>%
  mutate(blank_id = case_when(blank_type %in% c("LB", "PB") ~ library,
                              blank_type == "IB" ~ isolation_batch,
                              blank_type == "SB" ~ sampling_batch)) %>%
  group_by(blank_type, blank_id, sequence_id) %>%
  summarise(blank_reads = sum(reads_f), blank_pcrs = n_distinct(pcr_id), .groups = "drop")
blank_report <- map_dfr(names(blank_keys), function(bt) {
  key <- blank_keys[[bt]]
  summary_tbl %>%
    filter(.data[[paste0("in_", bt)]]) %>%
    transmute(blank_type = bt, blank_id = .data[[key]], sequence_id, sample_id, totread,
              removed_solely = n_blank_flags == 1)
})
if (nrow(blank_report) > 0) {
  blank_report <- blank_report %>%
    group_by(blank_type, blank_id, sequence_id) %>%
    summarise(n_samples_flagged = n_distinct(sample_id),
              n_samples_removed_solely = n_distinct(sample_id[removed_solely]),
              sample_reads_removed = sum(totread),
              samples = paste(sort(unique(sample_id)), collapse = ";"),
              .groups = "drop") %>%
    left_join(blank_reads, by = c("blank_type", "blank_id", "sequence_id")) %>%
    left_join(seqs %>% select(sequence_id, best_taxon, best_bestid, best_db), by = "sequence_id") %>%
    arrange(blank_type, blank_id, desc(sample_reads_removed))
}
write_tsv(blank_report, out("blank_report.tsv"))
cat("Blank report:", nrow(blank_report), "rows.\n")

# ------------------------------------------------------------------------------
# 10. Parameters of this run
# ------------------------------------------------------------------------------
write_json(list(
  script = "process_motu.R (MergeAndFilter design)",
  date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  motu_table = motu_file, demux_stats = demux_stats_json,
  classification_tables = as.list(db_files),
  reads_within = reads_within, reads_across = reads_across,
  reads_replicates = reads_replicates, min_identity = min_identity,
  n_pcrs = nrow(pcrs), n_samples = nrow(sample_meta),
  n_sequences = nrow(seqs), n_sequences_retained = length(retained_ids)
), out("filter_params.json"), auto_unbox = TRUE, pretty = TRUE)

cat("Done.\n")
cat("------------------------------------------------------------\n")
