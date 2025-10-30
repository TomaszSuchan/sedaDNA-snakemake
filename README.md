# sedaDNA Metabarcoding Pipeline - Complete Workflow Guide

## Overview

This pipeline processes sedimentary ancient DNA (sedaDNA) metabarcoding data using OBITools, applying multi-stage quality filtering, blank contamination control, and taxonomic classification. The workflow handles multiple projects, libraries, and reference databases, with comprehensive quality control at each step.

---

## 1. Configuration

### Edit the Configuration File

Edit `config/config.yaml` to define your project structure, input files, and parameters:
```yaml
projects:
  project_name:
    libraries:
      library1:
        forward: "path/to/forward.fastq.gz"
        reverse: "path/to/reverse.fastq.gz"
        barcode_file: "path/to/barcodes.csv"
    parameters:
      reference_dbs:
        PhyloAlps: "path/to/phyloalps.fasta"
        EMBL143: "path/to/embl143.fasta"
      seq_filters:
        reads_within: 3      # Minimum reads per replicate
        reads_across: 10     # Minimum total reads
        reads_replicates: 3  # Minimum number of replicates
      tax_filters:
        min_identity: 0.97   # Minimum taxonomic match (97%)
```

### Prepare Barcode Files

OBITools requires separate processing for barcodes of different lengths. The pipeline automatically detects barcode lengths and splits processing accordingly.

**Barcode file format** (`barcodes.csv`):
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
- `results/{project}/barcodes-{library}_{length}bp_only.txt`
- `results/{project}/{library}.barcode_validation.txt`

### Stage 2: Read Pairing and Quality Filtering

**Rule:** `pair_reads`

**Process:**
1. Pairs forward and reverse reads using `obipairing`
2. Keeps only successfully merged reads (`mode="alignment"`)
3. Filters by alignment quality (identity, overlap)

**Parameters:**
- `gap-penalty`: Gap penalty for alignment (default: 2.0)
- `min-identity`: Minimum identity between F/R reads (default: 0.9)
- `min-overlap`: Minimum overlap length (default: 20bp)

**Output:**
- `results/{project}/{library}.paired.fastq.gz`

### Stage 3: Demultiplexing

**Rule:** `demultiplex`

**Process:**
1. Demultiplexes reads by sample barcode using `obimultiplex`
2. Processes each barcode length separately
3. Concatenates results across barcode lengths

**Output:**
- `results/{project}/{library}.demux_{length}bp.fastq.gz` (per length)
- `results/{project}/{library}.demux.fastq.gz` (concatenated)
- `stats/{project}/{library}.demux_stats.json` (read counts per sample)

### Stage 4: Dereplication and Filtering

**Rules:** `dereplicate`, `filter_counts`, `denoise`

**Process:**
1. **Dereplication** (`obiuniq`): Collapses identical sequences, counts abundance
2. **Count filtering** (`obigrep`):
   - Removes sequences with `< min-count` reads (default: 2)
   - Removes sequences `< min-length` bp (default: 10)
3. **Denoising** (`obiclean`):
   - Identifies PCR/sequencing errors using sequence similarity
   - Flags low-abundance variants as potential errors
   - Optional chimera detection

**Output:**
- `results/{project}/{library}.demux.uniq.filtered.denoised.fasta.gz`

### Stage 5: Library Merging and Classification

**Rules:** `merge_all_libraries`, `classify`

**Process:**
1. Merges all libraries within a project
2. Performs taxonomic classification against each reference database using `obitag`
3. Creates MOTU tables and classification tables

**Output:**
- `results-classified/{project}/{project}-merged.fasta.gz`
- `results-classified/{project}/{project}-{db}.classified.fasta`
- `results-classified/{project}/{project}-{db}.motu_table.csv`
- `results-classified/{project}/{project}-{db}.classification_table.csv`

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

**Three filtering thresholds work together:**

1. **`reads_within`** (default: 3):
   - Minimum reads required per replicate
   - Filters out low-abundance observations within each replicate
   - Only replicates meeting this threshold are counted

2. **`reads_across`** (default: 10):
   - Minimum total reads after within-replicate filtering
   - Ensures sequences have substantial representation overall

3. **`reads_replicates`** (default: 3):
   - Minimum number of replicates that must pass `reads_within` threshold
   - Ensures reproducibility across technical replicates

**Example with reads_within=3, reads_replicates=3, reads_across=10:**
```
Sequence A in Sample 1:
  Rep1: 5 reads, Rep2: 2 reads, Rep3: 6 reads, Rep4: 4 reads
  → Only Rep1, Rep3, Rep4 counted (≥3 reads)
  → total_reads = 15, n_replicates_present = 3
  → KEEP (3 reps ≥ threshold, total ≥ 10)

Sequence B in Sample 2:
  Rep1: 4 reads, Rep2: 2 reads, Rep3: 3 reads, Rep4: 1 read
  → Only Rep1, Rep3 counted (≥3 reads)
  → total_reads = 7, n_replicates_present = 2
  → REMOVE (only 2 reps ≥ threshold, need 3)
```

### Blank Contamination Flagging

Sequences are flagged if present in any blank control:

| Flag | Description | Grouping |
|------|-------------|----------|
| `in_LB` | Library blank contamination | library + sequence |
| `in_PB` | PCR blank contamination | isolation_batch + sequence |
| `in_IB` | Isolation blank contamination | isolation_batch + sequence |
| `in_SB` | Sampling blank contamination | sampling_batch + sequence |

**Removal criteria:**
```
remove = not_replicated OR in_LB OR in_PB OR in_IB OR in_SB
```

A sequence is flagged for removal if **ANY** condition is met:
- Fails replication requirements
- Appears in any blank control

### Output Files

**Main table:** `results-tables/{project}/{project}-{db}-combined_classification_table.csv`

**Key columns:**

| Column | Description |
|--------|-------------|
| `core`, `depth` | Sample location identifiers |
| `library`, `replicate` | Sequencing identifiers |
| `blank_type` | Sample type (SAMPLE/LB/PB/IB/SB) |
| `sequence_id` | Unique sequence ID (e.g., seq0001) |
| `total_reads` | Sum of reads from valid replicates |
| `n_replicates_present` | Number of replicates ≥ threshold |
| `weighted_avg_proportion` | Weighted proportion across replicates |
| `replicate_summary` | Semicolon-separated read counts |
| `proportion_summary` | Semicolon-separated proportions |
| `not_replicated` | TRUE if replication criteria not met |
| `in_LB/PB/IB/SB` | TRUE if present in respective blank |
| `remove` | TRUE if flagged for removal |
| `obitag_bestid` | Taxonomic match identity (0-1) |
| `taxid` | NCBI taxonomy ID |
| `obitag_rank` | Taxonomic rank |
| `taxon` | Taxonomic name |

**Processing log:** `results-tables/{project}/{project}-{db}-combined_classification_info.txt`

Contains statistics on filtering at each stage.

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
```r
data %>% filter(obitag_bestid >= min_identity)
```

Filters by taxonomic match quality (default: 1 or 100%):

#### Step 3: Taxonomic Rank Filtering

Keeps only ranks ≥ family level:
```r
data %>% filter(obitag_rank %in% c("species", "subgenus", "section", 
                                    "genus", "family", "subfamily", "tribe"))
```

This only keeps the following ranks: species, subgenus, section, genus, family,
subfamily, tribe. Anything above is considered too broad.

**Rationale:** Family-level is minimum for meaningful ecological interpretation.

#### Step 4: Taxonomic Clustering

Aggregates by unique taxon at each location:
```r
data %>%
  group_by(core, depth, taxid, taxon, obitag_rank) %>%
  summarise(total_reads = sum(total_reads))
```

**Why cluster?**
- Multiple sequences can represent the same taxon (intraspecific variation)
- Different genetic regions may yield different sequences
- Clustering gives true taxon-level abundance

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

**Clustered table:** `results-tables/{project}/{project}-{db}-clustered_taxa_table.csv`

**Columns:**
- `core`, `depth`: Sample location
- `taxid`: NCBI taxonomy ID
- `taxon`: Taxonomic name
- `obitag_rank`: Taxonomic rank
- `total_reads`: Aggregated abundance

**Processing log:** `logs/{project}/{project}-{db}-clustered_taxa_table.log`

---

## 6. Visualization

### Taxa Heatmaps

**Rules:** `plot_taxa_heatmap`, `plot_taxa_heatmap_log`

Generates heatmaps showing taxonomic composition across samples:
- Linear scale: Shows absolute abundance patterns
- Log scale: Better for visualizing rare taxa

**Parameters:**
- `top_n_taxa`: Number of most abundant taxa to display (default: 200)
- `width`, `height`: Plot dimensions in inches

**Output:**
- `final_plots/{project}/{project}-{db}-taxa_heatmap.pdf`
- `final_plots/{project}/{project}-{db}-taxa_heatmap_log.pdf`

---

## 7. Complete Filtering Pipeline Summary
```
Raw reads
    ↓ [Pairing & merging]
Paired reads
    ↓ [Demultiplexing]
Demultiplexed reads (per sample/replicate)
    ↓ [Dereplication & denoising]
Unique sequences: ~3,000
    ↓ [Taxonomic classification]
Classified sequences
    ↓ [process_motu.R: Replication + Blank filtering]
├─ Remove not_replicated: -40%
├─ Remove in_blanks: -20%
└─ Valid sequences: ~2,000 sequence-location combinations
    ↓ [cluster_taxa.R: Identity + Rank filtering]
├─ Filter identity (<97%): -15%
├─ Filter rank (<family): -10%
└─ Valid taxa: ~1,500
    ↓ [Clustering by taxon]
Final taxa-location combinations: ~500
```

---

## 8. Quality Control Checkpoints

### Sequence-Level QC (process_motu.R)
- ✓ Technical replication (minimum N replicates)
- ✓ Read abundance (within and across replicates)
- ✓ Blank contamination (LB, PB, IB, SB)
- ✓ Proportional weighting (accounts for depth variation)

### Taxonomic-Level QC (cluster_taxa.R)
- ✓ Match identity (≥97% similarity to reference)
- ✓ Taxonomic resolution (≥family level)
- ✓ Intraspecific clustering (multiple sequences → single taxon)

---

## 9. Parameter Recommendations

### Read Filtering
- **reads_within: 3** - Conservative, removes sporadic observations
- **reads_across: 10** - Ensures meaningful abundance
- **reads_replicates: 3** - Standard for technical replicates (adjust if using 4+ replicates)

### Taxonomic Filtering
- **min_identity: 0.97-0.99** - High confidence (recommended)
- **min_identity: 0.95-0.97** - Moderate confidence (more permissive)
- **min_identity: <0.95** - Not recommended (too many uncertain matches)

### Blank Controls
Always include:
- Library blanks (LB) - one per library
- PCR blanks (PB) - one per PCR batch
- Isolation blanks (IB) - one per extraction batch
- Sampling blanks (SB) - field controls (if applicable)

---

## 10. Output Directory Structure
```
results/
├── {project}/
│   ├── {library}.paired.fastq.gz
│   ├── {library}.demux.fastq.gz
│   ├── {library}.demux.uniq.filtered.denoised.fasta.gz
│   └── barcodes-{library}_{length}bp_only.txt

results-classified/
├── {project}/
│   ├── {project}-merged.fasta.gz
│   ├── {project}-{db}.classified.fasta
│   ├── {project}-{db}.motu_table.csv
│   └── {project}-{db}.classification_table.csv

results-tables/
├── {project}/
│   ├── {project}-{db}-combined_classification_table.csv
│   └── {project}-{db}-clustered_taxa_table.csv

final_plots/
├── {project}/
│   ├── {project}-{db}-taxa_heatmap.pdf
│   └── {project}-{db}-taxa_heatmap_log.pdf

stats/
├── {project}/
│   ├── {library}.raw_stats.json
│   ├── {library}.pair_stats.json
│   ├── {library}.demux_stats.json
│   └── {project}.demux_stats_combined.json

logs/
└── {project}/
    ├── {library}.paired.log
    ├── {library}.demux_{length}bp.log
    └── {project}-{db}-combined_classification_table.log
```

---

## 11. Troubleshooting

### Issue: Samples not matching between MOTU and JSON
**Symptom:** Warning in log: "X samples in MOTU table not found in demux stats"
**Solution:** Check sample naming consistency between barcode file and demux stats

### Issue: All sequences flagged as remove=TRUE
**Symptom:** No sequences in final output
**Solution:** 
- Check replication parameters are appropriate for your data
- Verify blank samples are correctly labeled
- Review `reads_within`, `reads_across`, `reads_replicates` thresholds

### Issue: Low taxonomic assignment rate
**Symptom:** Many sequences with low `obitag_bestid` or "no rank"
**Solution:**
- Check reference database quality and coverage
- Consider lowering `min_identity` threshold (not below 0.95)
- Verify primers target appropriate region for your database

### Issue: High blank contamination
**Symptom:** Many sequences flagged with in_LB/PB/IB/SB=TRUE
**Solution:**
- Review laboratory procedures
- Check blank controls for unusual contamination
- Consider stricter laboratory protocols for future samples

---

## 12. Citation

If using this pipeline, please cite:
- OBITools: Boyer et al. (2016) Mol Ecol Resour 16:176-182
- MergeAndFilter approach: Lammers et al. (2021) Mol Ecol Resour
- Your reference databases (e.g., PhyloAlps, EMBL, etc.)

---

## 13. Support

For issues or questions:
1. Check the logs in `logs/{project}/` for error messages
2. Review the processing statistics in log files
3. Verify configuration file syntax and file paths
4. Consult the Snakemake and OBITools documentation