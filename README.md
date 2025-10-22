# Complete OBI-tools Workflow Example

This example shows how to use the optimized pipeline with library blanks and isolation blanks.

## 1. Edit the configuration file
Edit `config/config.yaml`, define location of your reads, barcodes, and name of the project.

## 2. Prepare Your Barcode Files
OBItools does not work with barcodes of different lengths, thus in the first part of the pipeline, the barcode file is divided by barcode lengths and the raw files are demultilexed separately for barcodes of different lengths, then results are concatenated.

Example barcode file (note that it does not contain @headers, these are defined in the confg file):
```csv
experiment,sample,sample_tag,forward_primer,reverse_primer
ZSG3,LB_1,AACAAGCC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,LB_2,TGAGAGCT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,LB_3,ACAACCGA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,LB_4,TTACGCCA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_1,GGATAGCA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_2,ACACACAG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_3,TCCAACAC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,PB_4,CAACCTCA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_SB_D1_isol1_1,GTGTAGTC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_SB_D1_isol1_2,TATCTGGC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_SB_D1_isol1_3,GTGTGTGT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_SB_D1_isol1_4,CTGTCAAC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,IB_isol1_1,TGAGTTCCT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,IB_isol1_2,TCTGGTTGA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,IB_isol1_3,TCACGGATA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,IB_isol1_4,TTCCACCTA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_025_D1_isol1_1,AGGAATGAG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_025_D1_isol1_2,ACTGACCTT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_025_D1_isol1_3,ATGAGCCTA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_025_D1_isol1_4,ATCATAGCG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_041_D1_isol1_1,AATTGCCG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_041_D1_isol1_2,ATGCTTGG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_041_D1_isol1_3,ATGGAGGT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_041_D1_isol1_4,TGAGGACA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_057_D1_isol1_1,CCGACCATA,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_057_D1_isol1_2,CAACACCGT,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_057_D1_isol1_3,CTCATACGC,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
ZSG3,ZSG3_057_D1_isol1_4,CAACAGGAG,GGGCAATCCTGAGCCAA,CCATTGAGTCTCTGCACCTATC
```

## 3. Run the Pipeline

```bash
# First, do a dry run
snakemake -n

# Run with 4 cores
snakemake --cores 4

# The pipeline will:
# 1. Automatically detect barcode lengths (8bp, 9bp, etc.)
# 2. Split barcode files by length
# 3. Process each library separately
# 4. Merge all libraries
# 5. Dereplicate, filer by read counts, denoise
# 6. Perform taxonomic classification
```

## 4. Multi-Project Batch Processing (NEW!)

The pipeline now supports processing multiple projects simultaneously with cross-project comparative analysis.

### 4.1 Setup Multi-Project Mode

Create or edit `config/projects.yaml`:

```yaml
mode: "multi"  # Set to "multi" to enable batch processing

global:
  # Shared reference databases
  reference_dbs:
    PhyloAlps: "data/ref_db/Alps_GH_clean_edit.cl.final_20-07-21.fasta"
    # ... add more databases

  # Shared processing parameters
  parameters:
    max-cpu: 8
    # ... other parameters

projects:
  PROJECT1:
    libraries:
      LIB1:
        forward: "path/to/forward.fq.gz"
        reverse: "path/to/reverse.fq.gz"
        barcode_file: "path/to/barcodes.txt"

  PROJECT2:
    libraries:
      LIB2:
        forward: "path/to/forward.fq.gz"
        reverse: "path/to/reverse.fq.gz"
        barcode_file: "path/to/barcodes.txt"
    # Optional: project-specific parameter overrides
    parameters:
      seq_filters:
        reads_within: 5

meta_analysis:
  enabled: true
  comparisons:
    combined_heatmap: true
    diversity_metrics: true
    contamination_report: true
    taxa_overlap: true
```

### 4.2 Run Multi-Project Pipeline

```bash
# Make sure config/projects.yaml has mode: "multi"
snakemake -n  # Dry run to check

# Run with multiple cores
snakemake --cores 8
```

### 4.3 Multi-Project Outputs

When multi-project mode is enabled, additional outputs are generated:

```
meta_analysis/
├── combined_taxa_heatmap.pdf           # Combined heatmap across all projects
├── combined_taxa_heatmap_log.pdf       # Log-scale version
├── taxa_overlap_venn.pdf               # Venn diagram of shared taxa
├── diversity_metrics.csv               # Alpha/beta diversity statistics
└── contamination_report.csv            # Cross-project contamination summary
```

### 4.4 Features

**Cross-Project Comparisons:**
- Combined taxonomic heatmaps with hierarchical clustering
- Taxa overlap analysis (Venn diagrams)
- Diversity metrics (Shannon, Simpson, Bray-Curtis)
- Contamination patterns across projects

**Flexible Configuration:**
- Global settings shared across projects
- Project-specific parameter overrides
- Shared reference databases

**Resource Optimization:**
- Single taxonomy download
- Shared reference databases
- Parallel processing of projects

### 4.5 Switching Between Modes

**Single-project mode** (default, backward compatible):
```yaml
# config/projects.yaml
mode: "single"
```
Or simply use `config/config.yaml` as before.

**Multi-project mode**:
```yaml
# config/projects.yaml
mode: "multi"
```

The pipeline automatically detects the mode and processes accordingly.

## 5. Output Files

### Single-Project Outputs
- `results/{project}/`: Raw processing outputs
- `stats/{project}/`: Read statistics
- `final_tables/{project}/`: Analysis tables
- `final_plots/{project}/`: Visualizations

### Multi-Project Outputs (additional)
- `meta_analysis/`: Cross-project comparisons
