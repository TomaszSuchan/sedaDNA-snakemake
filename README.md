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
    reads_within: 3        # <INTEGER> minimum number of reads within a single replicate to count that replicate (default 3)
    reads_across: 10       # <INTEGER> minimum total number of reads across all valid replicates (default 10)
    reads_replicates: 3    # <INTEGER> minimum number of replicates that must meet reads_within threshold (default 3)
  
  # Final filters for taxonomic assignments to remove low confidence assignments:
  tax_filters:
    min_identity: 0.97     # <FLOAT64> minimum taxonomic match identity to keep a sequence assignment (0-1, default 0.97)
  
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

**Rule:** `validate_barcodes`, `split_barcodes`

**Process:**
1. Validates barcode file format (required columns, duplicate checks)
2. Analyzes barcode length distribution
3. Splits barcode files by length (e.g., 8bp, 9bp)
4. Creates OBITools-compatible barcode files for each length

**Output:**
- `results/{project}/sequences/barcodes-{library}_{length}bp_only.txt`
- `results/{project}/validation/{library}.barcode_validation.txt`
- `results/{project}/validation/{project}.sample_replicates.tsv`
- `results/{project}/validation/{project}.pb_ib_per_library.tsv`
- `results/{project}/validation/{project}.sb_per_sampling.tsv`
- `results/{project}/validation/{project}.ib_per_isolation.tsv`
- `results/{project}/validation/{project}.duplicate_sample_identity.tsv`
- `results/{project}/validation/{project}.sample_name_errors.tsv`

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

The `process_motu.R` script applies multi-stage quality filtering, integrating read count statistics to calculate weighted proportions and identify contamination from blank controls.

### Sample Name Parsing

Each sample ID is automatically parsed into components:
```
ZSG-3_025_D1_250129_ZSG3MET1_1
  ↓
core: ZSG-3
depth: 025_D1
sampling_batch: ZSG-3_250129
isolation_batch: 250129
library: ZSG3MET1
replicate: 1
blank_type: SAMPLE
```

**Blank type detection:**
- **LB**: Library blank (isolation_batch starts with "LB")
- **PB**: PCR blank (isolation_batch starts with "PB")
- **IB**: Isolation blank (sampling_batch starts with "IB")
- **SB**: Sampling blank (depth starts with "SB")
- **SAMPLE**: Environmental samples (everything else)

### Proportional Abundance Calculation

For each sequence in each replicate:
```
proportion = reads_in_replicate / total_reads_in_sample
```

This accounts for sequencing depth variation between replicates, preventing samples with more reads from artificially dominating the analysis.

### Weighted Average Proportion

A **weighted average proportion** is calculated across replicates:
```
weighted_avg_proportion = Σ(proportion_i × total_reads_i) / Σ(total_reads_i)
```

This gives appropriate weight to replicates with higher sequencing depth (MergeAndFilter approach).

### Read Count Filtering

**Parameters** (from `config["parameters"]["seq_filters"]`):

1. **reads_within** `<INTEGER>` (default: 3):
   - Minimum number of reads within a single replicate to count that replicate
   - Filters out low-abundance observations within each replicate
   - Only replicates meeting this threshold contribute to total_reads

2. **reads_across** `<INTEGER>` (default: 10):
   - Minimum total number of reads across all valid replicates
   - Ensures sequences have substantial representation overall
   - Applied after within-replicate filtering

3. **reads_replicates** `<INTEGER>` (default: 3):
   - Minimum number of replicates that must meet reads_within threshold
   - Ensures reproducibility across technical replicates
   - Adjust if using different number of replicates (e.g., 4)

**Example with reads_within=3, reads_replicates=3, reads_across=10:**
```
Sequence A in Sample 1:
  Rep1: 5 reads, Rep2: 2 reads, Rep3: 6 reads, Rep4: 4 reads
  → Only Rep1, Rep3, Rep4 counted (≥3 reads each)
  → total_reads = 15, n_replicates_present = 3
  → KEEP (3 reps ≥ threshold, total ≥ 10)

Sequence B in Sample 2:
  Rep1: 4 reads, Rep2: 2 reads, Rep3: 3 reads, Rep4: 1 read
  → Only Rep1, Rep3 counted (≥3 reads each)
  → total_reads = 7, n_replicates_present = 2
  → REMOVE (only 2 reps ≥ threshold, need 3)
```

### Blank Contamination Flagging

Sequences are flagged if present in any blank control:

| Flag | Description | 
|------|-------------|
| `in_LB` | Library blank contamination |
| `in_PB` | PCR blank contamination |
| `in_IB` | Isolation blank contamination |
| `in_SB` | Sampling blank contamination |

**Removal criteria:**
```
remove = not_replicated OR in_LB OR in_PB OR in_IB OR in_SB
```

A sequence is flagged for removal if **ANY** condition is met:
- Fails replication requirements
- Appears in any blank control

### Output Files

**Main table:** `results/{project}/tables/{project}-{db}-classification_table.csv`

**Key columns:**

| Column | Description |
|--------|-------------|
| `core`, `depth` | Sample location identifiers |
| `sampling_batch` | Sampling batch ID |
| `isolation_batch` | Isolation batch ID |
| `library` | Library ID |
| `replicate` | Technical replicate of the sample within the library |
| `blank_type` | Sample type (SAMPLE/LB/PB/IB/SB) |
| `sequence_id` | Unique sequence ID (e.g., seq0001) |
| `total_reads` | Sum of reads from valid replicates after within-replicate filtering |
| `n_replicates_present` | Number of replicates meeting the reads_within threshold for this sequence |
| `weighted_avg_proportion` | Weighted average proportion across replicates (Σ(proportion_i × total_reads_i) / Σ(total_reads_i)) |
| `mean_proportion` | Simple mean of replicate proportions across replicates (unweighted) |
| `replicate_summary` | Semicolon-separated read counts per replicate (human-readable summary) |
| `proportion_summary` | Semicolon-separated proportions per replicate (human-readable summary) |
| `not_replicated` | TRUE if the sequence fails the replication criteria (e.g., fewer than reads_replicates) |
| `in_SB` | TRUE if the sequence is detected in any sampling blank (SB) |
| `in_IB` | TRUE if the sequence is detected in any isolation blank (IB) |
| `in_LB` | TRUE if the sequence is detected in any library blank (LB) |
| `in_PB` | TRUE if the sequence is detected in any PCR blank (PB) |
| `remove` | TRUE if flagged for removal (not_replicated OR in_LB OR in_PB OR in_IB OR in_SB) |
| `obitag_bestid` | Best taxonomic match identity returned by obitag (0–1) |
| `taxid` | NCBI taxonomy identifier for the assigned taxon |
| `obitag_rank` | Taxonomic rank returned by obitag (e.g., species, genus, family) |
| `taxon` | Taxonomic name assigned by obitag |

**Processing info:** `results/{project}/tables/{project}-{db}-classification_info.txt`

Contains statistics on filtering at each stage.

**Error log:** `logs/{project}/{project}-{db}-classification_table.log`

Contains error messages (if any).

---

## 5. Taxonomic Clustering

### Overview

The `cluster_taxa.R` script aggregates sequences by taxonomic assignment, applies identity and rank filters, and clusters by taxon at each sampling location.

### Processing Steps

#### Step 1: Remove Flagged Sequences
```r
data %>% filter(remove == FALSE)
```
Removes sequences flagged during MOTU processing.

#### Step 2: Taxonomic Identity Filtering

**Parameter** (from `config["parameters"]["tax_filters"]`):
- **min_identity** `<FLOAT64>` (default: 1.0): Minimum taxonomic match identity to keep a sequence assignment (0-1 scale)
```r
data %>% filter(obitag_bestid >= min_identity)
```

#### Step 3: Taxonomic Rank Filtering

Keeps only ranks at or above family level:
```r
data %>% filter(obitag_rank %in% c("species", "subgenus", "section", 
                                    "genus", "family", "subfamily", "tribe"))
```

Retained ranks are: species, subgenus, section, genus, subfamily, tribe, family.
Everything higher than family-level is discarded.

#### Step 4: Taxonomic Clustering

Aggregates by unique taxon at each location:
```r
data %>%
  group_by(core, depth, taxid, taxon, obitag_rank) %>%
  summarise(total_reads = sum(total_reads))
```

Multiple sequences representing the same taxon (intraspecific variation) are clustered together. This provides  taxon-level abundance.

**Important:** Clustering does not take taxonomy into consideration. So if you have sequences clustered on the genus level and the species within, you would have to add these togehter to get the total number of reads within this genus. In other words, the sum of reads assigned to higher taxonomic levels do not include reads assigned to lower taxonomic levels within.

**Example:**
```
Before:
core    depth   taxon            sequence  reads
ZSG-3   025_D1  Picea abies      seq0042   145
ZSG-3   025_D1  Picea abies      seq0089   89
ZSG-3   025_D1  Picea abies      seq0156   52

After:
core    depth   taxon            total_reads
ZSG-3   025_D1  Picea abies      286
```

### Output Files

**Clustered table:** `results/{project}/tables/{project}-{db}-clustered_taxa_table.csv`

**Columns:**
- `core`, `depth`: Sample location
- `taxid`: NCBI taxonomy ID
- `taxon`: Taxonomic name
- `obitag_rank`: Taxonomic rank
- `total_reads`: Aggregated abundance across all sequences

**Processing info:** `results/{project}/tables/{project}-{db}-clustered_taxa_info.txt`

Contains statistics on filtering and clustering.

**Error log:** `logs/{project}/{project}-{db}-clustered_taxa_table.log`

Contains error messages (if any).

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
Clean unique sequences: ~3,000
    ↓ [obitag: Taxonomic classification]
Classified sequences
    ↓ [process_motu.R: Replication + Blank filtering]
├─ reads_within filter: Remove low-abundance within replicates
├─ reads_replicates filter: Require minimum replication
├─ reads_across filter: Require minimum total reads
├─ Blank contamination filter: Remove if in LB/PB/IB/SB
└─ Valid sequences
    ↓ [cluster_taxa.R: Identity + Rank filtering]
├─ min_identity filter: Remove low-confidence assignments (<97%)
├─ Rank filter: Remove assignments above family level
└─ Valid taxa
    ↓ [Clustering by taxon]
Final taxa-location combinations
```

---

## 7. Output Directory Structure
```
results/
└── {project}/
    ├── validation/
    │   ├── {library}.barcode_validation.txt
    │   ├── {project}.sample_replicates.tsv
    │   ├── {project}.pb_ib_per_library.tsv
    │   ├── {project}.sb_per_sampling.tsv
    │   ├── {project}.ib_per_isolation.tsv
    │   ├── {project}.duplicate_sample_identity.tsv
    │   └── {project}.sample_name_errors.tsv
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
    │   ├── {project}-{db}-classification_table.csv
    │   ├── {project}-{db}-classification_info.txt
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
    ├── {project}-{db}-classification_table.log
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
