# sedaDNA Metabarcoding Pipeline

## Overview

This pipeline processes sedimentary ancient DNA (sedaDNA) metabarcoding data using OBITools, applying multi-stage quality filtering, blank contamination control, and taxonomic classification. The workflow handles multiple projects, libraries, and reference databases.

---

## 1. Configuration

### Edit the Configuration File

Edit `config/config.yaml` to define your project structure, input files, and parameters:
```yaml
## Global processing parameters
# These parameters will be applied to all projects and libraries,
# unless overridden in the specific project or library.

parameters: &global_params
  # General settings:
  max-cpu: 8            # <INTEGER> maximum number of CPU threads to use for processing (default 1)
  mem-per-cpu: 3850     # <INTEGER> amount of memory (in MB) to allocate per CPU thread (default 3850 - for plgrid ARES cluster)
  
  # Reference databases for taxonomic assignment:
  reference_dbs:
    PhyloAlps: "data/ref_db/Alps_GH_clean_edit.cl.final_20-07-21.fasta"
    ArctBorBryo: "data/ref_db/arctborbryo.gh.fixed.final_20-07-21.fasta"
    EMBL143: "data/ref_db/GH_143_no-env_clean.final_20-07-21.fasta"
    PhyloNorway: "data/ref_db/Norway_GH_clean_edit.cl.final_20-07-21.fasta"
  
  # Merging parameters for obipairing:
  obipairing:
    gap-penalty: 2.0       # <FLOAT64> gap penalty expressed as the multiply factor applied to the mismatch score between two nucleotides with a quality of 40 (default 2.0)
    min-identity: 0.9      # <FLOAT64> minimum identity between overlapped regions of the reads to consider the alignment (default 0.9)
    min-overlap: 20        # <INTEGER> minimum overlap between both the reads to consider the alignment (default 20)
    penalty-scale: 1.0     # <FLOAT64> scale factor applied to the mismatch score and the gap penalty (default 1.0)
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to merge a read pair (default 60)
  
  # Demultiplexing parameters:
  obimultiplex:
    matching: strict       # <STRING> matching mode: strict or fuzzy (default strict)
    primer_mismatches: 2   # <INTEGER> maximum number of mismatches allowed in primer sequences (default 2)
    indels: false          # <BOOLEAN> if true, allow insertions/deletions in primer matching (default false)
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to demultiplex a library (default 60)
  
  # Dereplication parameters:
  obiuniq:
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to dereplicate sequences (default 60)
  
  # Minimum count filtering:
  filtering:
    min-count: 2           # <INTEGER> minimum count to keep a sequence (default 2)
    min-length: 10         # <INTEGER> minimum length to keep a sequence (default 10)
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to filter sequences (default 60)
  
  # Filtering parameters for obiclean:
  obiclean:
    ratio: 0.05            # <FLOAT64> ratio to consider a sequence as an error of another sequence (default 0.05)
    distance: 1            # <INTEGER> maximum distance to consider two sequences as related (default 1)
    detect_chimera: false  # <BOOLEAN> if true, detect and flag chimeric sequences (default false)
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to clean sequences (default 60)
  
  # Taxonomic assignment parameters for obitag:
  obitag:
    time: 60               # <INTEGER> maximum time (in minutes) to spend trying to classify a sequence (default 60)
  
  # Final filters on OTU tables to flag low-abundance and low-replicability sequences:
  seq_filters:
    reads_within: 3        # <INTEGER> minimum reads of a sequence in one PCR replicate; lower counts are set to 0 (default 3)
    reads_across: 10       # <INTEGER> minimum total reads of a sequence over all sample PCR replicates of the project (default 10)
    reads_replicates: 3    # <INTEGER> minimum number of sample PCR replicates (whole project) in which a sequence occurs (default 3)
  
  # Final filters for taxonomic assignments to remove low confidence assignments:
  tax_filters:
    min_identity: 0.97     # <FLOAT64> minimum identity, reached in at least one database to keep a sequence; each database's taxon table uses only its own assignments that reach it (0-1)
  
  # Plotting parameters:
  plotting:
    top_n_taxa: 200        # <INTEGER> number of top taxa to show in heatmap (default 200)
    width: 10              # <INTEGER> width of heatmap plot in inches (default 10)
    height: 15             # <INTEGER> height of heatmap plot in inches (default 15)

## Projects configuration
# Each project can have multiple reference databases and multiple libraries.

projects:
  p6_loop:
    libraries:
      ZSSG3_MET1:
        forward: "/path/to/forward.fq.gz"
        reverse: "/path/to/reverse.fq.gz"
        barcode_file: "/path/to/barcodes.txt"
      ZSSG3_MET5:
        forward: "/path/to/forward.fq.gz"
        reverse: "/path/to/reverse.fq.gz"
        barcode_file: "/path/to/barcodes.txt"
    parameters:
        <<: *global_params # inherit global parameters
        reference_dbs:
          PhyloAlps: "data/ref_db/custom_db.fasta"  # override global database
```

Parameters can be configured globally (under `parameters`) but can be overridden for specific projects as shown above.

### Prepare Barcode Files

OBITools requires separate processing for barcodes of different lengths. The pipeline automatically detects barcode lengths and splits processing accordingly.

**Barcode file format** (`barcodes.txt`):
```csv
experiment,sample,sample_tag,forward_primer,reverse_primer
ZSG3,LB_ZSG3MET1_1,AACAAGCC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,LB_ZSG3MET1_2,TGAGAGCT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_ZSG3MET1_1,GGATAGCA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_ZSG3MET1_2,ACACACAG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,IB_250129_ZSG3MET1_1,TGAGTTCCT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG-3_SB_D1_250129_ZSG3MET1_1,GTGTAGTC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG-3_025_D1_250129_ZSG3MET1_1,AGGAATGAG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG-3_025_D1_250129_ZSG3MET1_2,ACTGACCTT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
```

**Sample naming convention:**
- **Library blanks**: `LB_<library>_<replicate>` (e.g., `LB_ZSG3MET1_1`)
- **PCR blanks**: `PB_<library>_<replicate>` (e.g., `PB_ZSG3MET1_1`)
- **Isolation blanks**: `IB_<batch>_<library>_<replicate>` (e.g., `IB_250129_ZSG3MET1_1`)
- **Sampling blanks**: `<core>_SB_<horizon>_<batch>_<library>_<replicate>` (e.g., `ZSG-3_SB_D1_250129_ZSG3MET1_1`)
- **Environmental samples**: `<core>_<depth>_<horizon>_<batch>_<library>_<replicate>` (e.g., `ZSG-3_025_D1_250129_ZSG3MET1_1`)

The pipeline automatically identifies blank types based on these naming patterns.

---

## 2. Run the Pipeline

### Dry Run (Test Configuration)
```bash
snakemake -n
```

### Local Execution
```bash
# Run with 4 cores
snakemake --cores 4 --use-conda
```

### HPC Execution (SLURM)
```bash
# Run on HPC cluster with 48 parallel jobs
snakemake --executor slurm \
  --default-resources slurm_account=your-account slurm_partition=your-partition \
  --jobs 48 --use-conda --retries 3 --rerun-incomplete
```

---

## 3. Pipeline Stages

### Stage 1: Barcode Validation and Splitting

**Rule:** `validate_barcodes`, `validate_samples`, `split_barcodes`

**Process:**
1. Validates barcode file format (required columns, duplicate checks)
2. Analyzes barcode length distribution
3. Splits barcode files by length (e.g., 8bp, 9bp)
4. Creates OBITools-compatible barcode files for each length
5. Writes validation reports and stores PASS/FAIL status for sample naming and duplicate identity checks

**Output:**
- `results/{project}/sequences/barcodes-{library}_{length}bp_only.txt`
- `results/{project}/validation/{library}.barcode_validation.txt`
- `results/{project}/validation/{project}.all_replicates.tsv`
- `results/{project}/validation/{project}.sample_replicates.tsv`
- `results/{project}/validation/{project}.LB_PB_replicates.tsv`
- `results/{project}/validation/{project}.SB_replicates.tsv`
- `results/{project}/validation/{project}.IB_replicates.tsv`
- `results/{project}/validation/{project}.duplicate_sample_identity.tsv`
- `results/{project}/validation/{project}.sample_name_errors.tsv`
- `results/{project}/validation/{project}.validation_passed.txt`
- `results/{project}/validation/{project}.validation_asserted.txt` (created in full runs after all downstream rules finish)

### Stage 2: Read Pairing

**Rule:** `pair_reads`

**Process:**
1. Pairs forward and reverse reads using `obipairing`
2. Applies quality-based alignment parameters
3. Keeps only successfully merged reads (filters by `mode="alignment"`)

**Parameters** (from `config["parameters"]["obipairing"]`):
- **gap-penalty** `<FLOAT64>`: Gap penalty expressed as multiply factor applied to mismatch score between two nucleotides with quality of 40 (default: 2.0)
- **min-identity** `<FLOAT64>`: Minimum identity between overlapped regions to consider alignment valid (default: 0.9)
- **min-overlap** `<INTEGER>`: Minimum overlap length between forward and reverse reads (default: 20 bp)
- **penalty-scale** `<FLOAT64>`: Scale factor applied to mismatch score and gap penalty (default: 1.0)
- **time** `<INTEGER>`: Maximum time in minutes to spend merging a read pair (default: 60)

**Output:**
- `results/{project}/sequences/{library}.paired.fastq.gz`

### Stage 3: Demultiplexing

**Rule:** `demultiplex`

**Process:**
1. Demultiplexes reads by sample barcode using `obimultiplex`
2. Processes each barcode length separately
3. Concatenates results across barcode lengths
4. Generates read count statistics per sample

**Parameters** (from `config["parameters"]["obimultiplex"]`):
- **matching** `<STRING>`: Matching mode for barcode/primer identification - "strict" or "fuzzy" (default: strict)
- **primer_mismatches** `<INTEGER>`: Maximum number of mismatches allowed in primer sequences (default: 2)
- **indels** `<BOOLEAN>`: If true, allow insertions/deletions in primer matching (default: false)
- **time** `<INTEGER>`: Maximum time in minutes to spend demultiplexing a library (default: 60)

**Output:**
- `results/{project}/sequences/{library}.demux_{length}bp.fastq.gz` (per barcode length)
- `results/{project}/sequences/{library}.demux.fastq.gz` (concatenated)
- `results/{project}/stats/{library}.demux_stats.json` (read counts per sample)

### Stage 4: Dereplication and Filtering

**Rules:** `dereplicate`, `filter_counts`, `denoise`

**Process:**

#### 4.1 Dereplication (`obiuniq`)
- Collapses identical sequences
- Counts abundance of each unique sequence
- Retains quality scores and metadata

**Parameters** (from `config["parameters"]["obiuniq"]`):
- **time** `<INTEGER>`: Maximum time in minutes to spend dereplicating sequences (default: 60)

#### 4.2 Count Filtering (`obigrep`)
- Removes low-abundance sequences (likely errors)
- Removes sequences below minimum length

**Parameters** (from `config["parameters"]["filtering"]`):
- **min-count** `<INTEGER>`: Minimum read count to keep a sequence (default: 2)
- **min-length** `<INTEGER>`: Minimum sequence length in bp to keep (default: 10)
- **time** `<INTEGER>`: Maximum time in minutes to spend filtering sequences (default: 60)

#### 4.3 Denoising (`obiclean`)
- Identifies PCR/sequencing errors using sequence similarity
- Flags low-abundance variants as potential errors of more abundant sequences
- Optional chimera detection

**Parameters** (from `config["parameters"]["obiclean"]`):
- **ratio** `<FLOAT64>`: Abundance ratio threshold to consider a sequence as potential error of another sequence (default: 0.05)
- **distance** `<INTEGER>`: Maximum edit distance to consider two sequences as potentially related (default: 1)
- **detect_chimera** `<BOOLEAN>`: If true, detect and flag chimeric sequences (default: false)
- **time** `<INTEGER>`: Maximum time in minutes to spend cleaning sequences (default: 60)

**Output:**
- `results/{project}/sequences/{library}.demux.uniq.filtered.denoised.fasta.gz`

### Stage 5: Library Merging and Classification

**Rules:** `merge_all_libraries`, `classify`

**Process:**
1. Merges all libraries within a project using `obiuniq`
2. Performs taxonomic classification against each reference database using `obitag`
3. Creates MOTU tables (read count matrices) and classification tables

**Parameters** (from `config["parameters"]["obitag"]`):
- **time** `<INTEGER>`: Maximum time in minutes to spend classifying a sequence (default: 60)

**Reference databases** (from `config["parameters"]["reference_dbs"]`):
- User-defined databases for taxonomic assignment
- Multiple databases can be specified per project

**Output:**
- `results/{project}/sequences/{project}-merged.fasta.gz`
- `results/{project}/classified/{project}-{db}.classified.fasta`
- `results/{project}/classified/{project}-{db}.motu_table.csv`
- `results/{project}/classified/{project}-{db}.classification_table.csv`

---

## 4. MOTU Processing and Filtering

### Overview

`process_motu.R` merges the MOTU table with the classifications from **all** reference databases and filters it. The design follows MergeAndFilter (Lammers, https://github.com/Y-Lammers/MergeAndFilter):

1. per observation (one sequence in one PCR replicate): counts below `reads_within` are set to 0;
2. per sequence, over the whole project: the sequence must have at least `reads_across` filtered reads and occur in at least `reads_replicates` sample PCR replicates;
3. per sequence: it is kept if its identity reaches `min_identity` in at least one database;
4. per sample: sequences found in a matching blank are flagged and removed;
5. per sample: replicate statistics (`totread`, `avgpread`, `sdpread`, `totrep`, `proprep`, `weightrep`);
6. pooling: all isolations and libraries of one sediment sample (`bio_sample` = `<core>_<depth>`) become one sample whose replicates are all their PCR replicates.

The script runs **once per project**, not once per database. The read and replicate filters do not depend on the database, so every database gets identical filtering, and the classification of each database is kept in its own columns (`bestid_<db>`, `taxid_<db>`, `rank_<db>`, `taxon_<db>`). The MOTU tables of all databases contain the same counts and sequence ids, so the first database's MOTU table is used.

### Sample Name Parsing

PCR names are parsed with one explicit pattern per type (the same formats as `validate_samples`). A name that matches none of them stops the run and is listed in `{project}-unparsed_pcr_names.txt`.

| Type | Format | Example | Blank key |
|------|--------|---------|-----------|
| Sample | `<core>_<depth>_<Dn>_<isolation>_<library>_<rep>` | `ZSG-3_025_D1_250129_ZSG3MET1_1` | – |
| SB | `<core>_SB_<Dn>_<isolation>_<library>_<rep>` | `ZSG-3_SB_D1_250129_ZSG3MET1_1` | `sampling_batch` = `<core>_<Dn>` |
| IB | `IB_<isolation>_<library>_<rep>` | `IB_250129_ZSG3MET1_1` | `isolation_batch` |
| LB | `LB_<library>_<rep>` | `LB_ZSG3MET1_1` | `library` |
| PB | `PB_<library>_<rep>` | `PB_ZSG3MET1_1` | `library` |

A *sample* is all PCR replicates sharing the name without `_<rep>`. `bio_sample` (`<core>_<depth>`, e.g. `ZSG-3_025`) is the sediment sample, one depth of one core; it can occur in several isolations or libraries, so it is the column to pool or plot by. `lab_sample_id` (`<core>_<depth>_<Dn>`, e.g. `ZSG-3_025_D1`) is the sample as named in the lab database.

The list of PCR replicates comes from the demultiplexing statistics, so a replicate without any reads still counts (as zeros) in `nrep`, `proprep` and `avgpread`.

### Filters

**Parameters** (from `config["parameters"]["seq_filters"]` and `["tax_filters"]`):

| Parameter | Default | Applies to | Meaning |
|-----------|---------|------------|---------|
| `reads_within` | 3 | one observation | counts below this in a PCR replicate are set to 0 |
| `reads_across` | 10 | one sequence, whole project | minimum sum of filtered reads over all sample PCR replicates |
| `reads_replicates` | 3 | one sequence, whole project | minimum number of sample PCR replicates with filtered reads > 0 |
| `min_identity` | 1.0 | one sequence | kept if reached in at least one database |

Blanks do not count towards `reads_across` and `reads_replicates`.

**Example** (`reads_within = 3`, `reads_across = 10`, `reads_replicates = 3`):
```
Sequence A: sample 1 reps 5, 2, 6, 4 reads; sample 2 reps 0, 3, 0, 0 reads
  -> after reads_within: 5, 0, 6, 4 and 0, 3, 0, 0
  -> total 18 reads (>= 10), 4 observations (>= 3): KEEP
Sequence B: sample 1 reps 4, 2, 3, 1 reads
  -> after reads_within: 4, 0, 3, 0
  -> total 7 reads (< 10), 2 observations (< 3): REMOVE
```

Replication within a single sample is not used as a filter; it is reported (`totrep`, `proprep`, `weightrep`) so that it can be applied later.

### Blank Contamination

After the filters above, a sequence is flagged in a sample when it has reads (≥ `reads_within`) in a blank with the same key:

| Flag | Blank | Matched by |
|------|-------|------------|
| `in_LB` | library blank | library |
| `in_PB` | PCR blank | library |
| `in_IB` | isolation blank | isolation batch |
| `in_SB` | sampling blank | core + sampling batch (`<core>_<Dn>`) |

`remove_blank` is TRUE when any flag is set. Blanks are not filtered against themselves. The log warns about blanks whose key matches no sample (usually a naming problem).

### Replicate Statistics (MergeAndFilter `_summary`)

For each sample and retained sequence, the replicate reads are stored in a fixed replicate order with zeros (`replicates` = `1;2;3;4`, `replicate_reads` = `50;0;40;0`), together with the retained reads of each replicate (`replicate_totals`, the denominator of the proportions). Because the order is fixed, sequences can be added replicate by replicate later (see section 5).

| Column | Meaning |
|--------|---------|
| `totread` | sum of the filtered reads over the replicates |
| `avgpread`, `sdpread` | mean and SD over the replicates of the proportion of the replicate's retained reads (replicates without the sequence count as 0) |
| `totrep`, `nrep` | replicates with the sequence, replicates of the sample |
| `proprep` | `totrep / nrep` |
| `weightrep` | sum of the weights of the replicates with the sequence; a replicate's weight is its share of the sample's retained reads. As in MergeAndFilter it is used only when the sample's mean `proprep` ≥ 0.33 and it has ≥ 10 sequences (`use_weighting`); otherwise `weightrep = proprep` |
| `prop_of_raw` | `totread` / raw demultiplexed reads of the sample |

### Pooling to the Sediment Sample

Re-isolations and re-amplifications of one depth end up as one sample. All PCR replicates of all isolation × library units of a `bio_sample` become its replicates (labelled `<Dn>_<isolation>_<library>_<rep>`), and the statistics above are recomputed over them.

Blank removal is done **before** pooling, at the level of each blank: a sequence found in the PCR blank of library A is set to 0 in the replicates of library A only, and keeps its reads from library B. In the same way SB acts only on the samples of its sampling session (`<core>_<Dn>`), IB only on its isolation batch, LB only on its library. `units_removed_by_blank` and `removed_by` show how many units lost the sequence and to which blank types.

### Output Files

All in `results/{project}/tables/`:

| File | Content |
|------|---------|
| `{project}-sample_summary.tsv` | one row per sample (isolation × library unit) × retained sequence, blanks included: sample fields, replicate statistics, blank flags, `remove_blank`, assignment of each database and the best one (`best_bestid`, `best_db`, `best_taxon`) |
| `{project}-bio_sample_summary.tsv` | one row per sediment sample (`bio_sample`, all isolations and libraries pooled) × sequence: pooled replicate statistics, `n_units`, `units_removed_by_blank`, `removed_by`, assignments; input of `cluster_taxa.R` |
| `{project}-sequences.tsv` | one row per sequence: reads in samples and blanks, observations, `pass_reads_across`, `pass_reads_replicates`, `pass_identity`, `retained`, assignment of each database |
| `{project}-pcrs.tsv` | one row per PCR replicate: parsed name, raw reads, reads after each step and the share of raw reads removed by the technical filters (including OBITools), identity, blanks, and retained |
| `{project}-samples.tsv` | one row per sample: reads, share retained, number of sequences, mean replication of the 10 most abundant sequences before (`avg_rep`) and after (`avg_filt_rep`) filtering, and their `overlap` |
| `{project}-blank_report.tsv` | per blank and sequence: reads in the blank, samples flagged, samples where this blank was the only reason, reads removed |
| `{project}-replicates.tsv.gz` | every observation (`pcr_id`, `sequence_id`, `reads`), unfiltered and including blanks, so filters can be re-applied without re-running the pipeline |
| `{project}-filter_params.json` | parameters and counts of the run |
| `{project}-process_motu_info.txt` | log with the counts at each step |

---

## 5. Taxonomic Clustering

### Overview

`cluster_taxa.R` makes one taxon table per database and sediment sample from `{project}-bio_sample_summary.tsv` (blank removal is already applied there).

### Processing Steps

1. Keep sequences whose assignment **in this database** has `bestid ≥ min_identity` (process_motu.R kept a sequence if any database reached it; here the database's own assignment must reach it).
2. Keep assignments at family level or below: family, subfamily, tribe, subtribe, genus, subgenus, section, subsection, series, species group, species subgroup, species, subspecies, varietas, forma. `clade` and `no rank` cannot be placed without the lineage and are dropped; the number of rows per rank is written to the log.
3. Merge the sequences of one taxon within each sediment sample. Replicate reads are added replicate by replicate, then `totread`, `totrep`, `proprep`, `avgpread`, `sdpread` and `weightrep` are recomputed.

Example:
```
sample  sequence  replicate_reads
S1      seq0042   50;0;40;0
S1      seq0089   0;30;0;20
-> Pinus  80;30;40;20   totrep = 4 (not 2)
```

Clustering does not use the taxonomy tree: reads assigned to a genus do not include the reads assigned to species within it.

### Output Files

**Clustered table:** `results/{project}/tables/{project}-{db}-clustered_taxa_table.csv`

`bio_sample`, `lab_sample_ids`, `core`, `depth`, `n_units`, `units`, `taxid`, `taxon`, `obitag_rank`, `n_sequences`, `sequence_ids`, `max_bestid`, `replicates`, `replicate_reads`, `replicate_totals`, `totread`, `totrep`, `nrep`, `proprep`, `avgpread`, `sdpread`, `weightrep`, `removed_by`. Blanks are not in this table (see `sample_summary` and `blank_report`).

**Processing info:** `results/{project}/tables/{project}-{db}-clustered_taxa_info.txt`

**Error log:** `logs/{project}/{project}-{db}-clustered_taxa_table.log`

---

## 6. Complete Filtering Pipeline Summary
```
Raw paired-end reads
    ↓ [obipairing: Pairing & merging]
Paired reads (aligned F+R)
    ↓ [obimultiplex: Demultiplexing by barcode]
Demultiplexed reads (per sample/replicate)
    ↓ [obiuniq: Dereplication]
Unique sequences with counts
    ↓ [obigrep: Count & length filtering]
Filtered unique sequences
    ↓ [obiclean: Denoising & error removal]
Clean unique sequences
    ↓ [obitag: Taxonomic classification, every database]
Classified sequences
    ↓ [process_motu.R: once per project, MergeAndFilter design]
├─ reads_within: counts < 3 in a replicate set to 0
├─ reads_across: ≥ 10 reads per sequence over all sample replicates
├─ reads_replicates: ≥ 3 sample replicates per sequence
├─ min_identity: reached in at least one database
├─ blanks: flag/remove per sample (LB, PB by library; IB by isolation; SB by core + Dn)
├─ per-sample replicate statistics
└─ pooling of isolations/libraries to bio_sample (core_depth)
    ↓ [cluster_taxa.R: per database]
├─ min_identity in this database
├─ ranks at family level or below
└─ sequences merged per taxon, replicate by replicate
Final sample × taxon tables
```

---

## 7. Output Directory Structure
```
results/
└── {project}/
    ├── validation/
    │   ├── {library}.barcode_validation.txt
    │   ├── {project}.all_replicates.tsv
    │   ├── {project}.sample_replicates.tsv
    │   ├── {project}.LB_PB_replicates.tsv
    │   ├── {project}.SB_replicates.tsv
    │   ├── {project}.IB_replicates.tsv
    │   ├── {project}.duplicate_sample_identity.tsv
    │   ├── {project}.sample_name_errors.tsv
    │   ├── {project}.validation_passed.txt
    │   └── {project}.validation_asserted.txt
    ├── sequences/
    │   ├── barcodes-{library}_{length}bp_only.txt
    │   ├── {library}.paired.fastq.gz
    │   ├── {library}.demux_{length}bp.fastq.gz
    │   ├── {library}.demux.fastq.gz
    │   ├── {library}.demux.uniq.fasta.gz
    │   ├── {library}.demux.uniq.filtered.fasta.gz
    │   ├── {library}.demux.uniq.filtered.denoised.fasta.gz
    │   └── {project}-merged.fasta.gz
    ├── classified/
    │   ├── {project}-{db}.classified.fasta
    │   ├── {project}-{db}.classified.no_annot.fasta
    │   ├── {project}-{db}.motu_table.csv
    │   └── {project}-{db}.classification_table.csv
    ├── tables/
    │   ├── {project}-sample_summary.tsv
    │   ├── {project}-bio_sample_summary.tsv
    │   ├── {project}-sequences.tsv
    │   ├── {project}-pcrs.tsv
    │   ├── {project}-samples.tsv
    │   ├── {project}-blank_report.tsv
    │   ├── {project}-replicates.tsv.gz
    │   ├── {project}-filter_params.json
    │   ├── {project}-process_motu_info.txt
    │   ├── {project}-{db}-clustered_taxa_table.csv
    │   └── {project}-{db}-clustered_taxa_info.txt
    └── stats/
        ├── {library}.raw_stats.json
        ├── {library}.pair_stats.json
        ├── {library}.demux_stats.json
        ├── {library}.merged_stats.tsv
        └── {project}.demux_stats_combined.json

logs/
└── {project}/
    ├── {library}.paired.log
    ├── {library}.demux_{length}bp.log
    ├── {project}-process_motu.log
    └── {project}-{db}-clustered_taxa_table.log
```

---

## 8. Citation

If using this pipeline, please cite:

- OBITools: Boyer F, et al. (2016) obitools: a unix-inspired software package for DNA metabarcoding. Molecular Ecology Resources 16:176-182

- Filtering approach: Lammers Y (2021) https://github.com/Y-Lammers/MergeAndFilter

- Your reference databases


---

## 9. Advanced Usage

### Running Specific Steps

Run only validation checks:
```bash
snakemake validate --cores 4
```

Run only classification for one database:
```bash
snakemake results/project/classified/project-PhyloAlps.motu_table.csv --cores 4
```

Run only final tables:
```bash
snakemake results/project/tables/project-PhyloAlps-clustered_taxa_table.csv --cores 4
```

### Reprocessing with Different Parameters

Edit config file and rerun with force:
```bash
snakemake --forcerun process_motu --cores 4
```

### Adding New Samples

1. Add new fastq and barcode files to config
2. Run pipeline (Snakemake will detect new inputs)

### Using Multiple Databases

Add new databases to config:
```yaml
reference_dbs:
  Database1: "path1.fasta"
  Database2: "path2.fasta"
```

Pipeline will automatically classify against all databases.
