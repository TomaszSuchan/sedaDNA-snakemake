# Multi-Project Batch Processing Guide

## Overview

The sedaDNA-snakemake pipeline now supports processing multiple projects simultaneously with automated cross-project comparative analysis. This feature enables:

- **Batch Processing**: Run multiple projects with a single command
- **Comparative Analysis**: Automatic generation of cross-project comparisons
- **Resource Optimization**: Share reference databases and taxonomies
- **Flexible Configuration**: Project-specific parameter overrides

## Quick Start

### 1. Create Multi-Project Configuration

```bash
# Copy the example configuration
cp config/projects.example.yaml config/projects.yaml

# Edit to add your projects
nano config/projects.yaml
```

### 2. Set Mode to Multi

Ensure the first line in `config/projects.yaml` is:
```yaml
mode: "multi"
```

### 3. Add Your Projects

```yaml
projects:
  MyProject1:
    libraries:
      LIB1:
        forward: "/path/to/forward_R1.fq.gz"
        reverse: "/path/to/forward_R2.fq.gz"
        barcode_file: "/path/to/barcodes.txt"
  
  MyProject2:
    libraries:
      LIB2:
        forward: "/path/to/forward_R1.fq.gz"
        reverse: "/path/to/forward_R2.fq.gz"
        barcode_file: "/path/to/barcodes.txt"
```

### 4. Run the Pipeline

```bash
# Dry run to verify
snakemake -n

# Run with 8 cores
snakemake --cores 8
```

## Configuration Structure

### Global Settings

Settings applied to all projects unless overridden:

```yaml
global:
  reference_dbs:
    Database1: "path/to/db1.fasta"
  
  parameters:
    max-cpu: 8
    # ... other parameters
  
  barcodes:
    matching: strict
```

### Project-Specific Settings

Override global settings for individual projects:

```yaml
projects:
  SpecialProject:
    libraries:
      # ... library definitions
    
    parameters:
      seq_filters:
        reads_within: 5  # Override global value
```

### Meta-Analysis Configuration

Control cross-project comparisons:

```yaml
meta_analysis:
  enabled: true
  
  filters:
    min_projects: 2        # Taxon must appear in 2+ projects
    min_total_reads: 50    # Minimum reads across all projects
    min_identity: 0.95     # Minimum taxonomic identity
  
  plotting:
    color_scheme: "Set2"
    width: 15
    height: 20
```

## Output Files

### Per-Project Outputs

Same as single-project mode:
```
results/PROJECT1/       # Raw processing results
stats/PROJECT1/         # Statistics
final_tables/PROJECT1/  # Analysis tables
final_plots/PROJECT1/   # Visualizations
```

### Meta-Analysis Outputs

New cross-project outputs:
```
meta_analysis/
├── all_projects_taxa.csv               # Combined taxa data
├── combined_taxa_heatmap.pdf           # Combined heatmap (linear)
├── combined_taxa_heatmap_log.pdf       # Combined heatmap (log scale)
├── taxa_overlap_venn.pdf               # Venn diagram of shared taxa
├── diversity_metrics.csv               # Alpha/beta diversity stats
└── contamination_report.csv            # Contamination summary
```

## Meta-Analysis Features

### 1. Combined Taxonomic Heatmap

Visualizes taxa abundance across all samples from all projects:
- Hierarchical clustering of projects and taxa
- Linear and log10 scales
- Color-coded by project
- Configurable number of top taxa

### 2. Taxa Overlap Analysis

Venn diagram showing:
- Unique taxa per project
- Shared taxa between projects
- Core taxa across all projects

### 3. Diversity Metrics

Calculates and compares:
- **Alpha diversity**: Shannon, Simpson, Richness, Evenness
- **Beta diversity**: Bray-Curtis dissimilarity
- Within-project vs. between-project comparisons

### 4. Contamination Report

Summarizes across all projects:
- Total sequences per project
- Contamination rates (% flagged)
- Blank type distribution
- Cross-contamination patterns

## Use Cases

### Use Case 1: Temporal Study

Process samples from multiple time points:
```yaml
projects:
  Site_2021:
    libraries: ...
  Site_2022:
    libraries: ...
  Site_2023:
    libraries: ...
```

**Benefits**: Track taxonomic changes over time

### Use Case 2: Multi-Site Comparison

Compare different geographic locations:
```yaml
projects:
  Alps_Site1:
    libraries: ...
  Alps_Site2:
    libraries: ...
  Carpathians_Site1:
    libraries: ...
```

**Benefits**: Identify biogeographic patterns

### Use Case 3: Method Comparison

Test different extraction or amplification methods:
```yaml
projects:
  Method_A:
    libraries: ...
  Method_B:
    libraries: ...
    parameters:
      seq_filters:
        reads_within: 5  # More stringent
```

**Benefits**: Validate methodological choices

## Best Practices

### 1. Consistent Naming

Use systematic project names:
- `Site_Year` format: `Alps_2023`, `Alps_2024`
- `Method_Replicate` format: `ExtractA_Rep1`, `ExtractA_Rep2`

### 2. Parameter Tuning

- Start with global defaults
- Override only when necessary
- Document reasons for project-specific changes

### 3. Reference Databases

- Use same databases across projects for comparability
- Update all projects together when databases change

### 4. Resource Allocation

- Allocate CPUs based on number of projects
- Consider memory requirements for meta-analysis
- Use `--cores` wisely (2-4 cores per project)

## Troubleshooting

### Issue: "No such file or directory"

**Solution**: Check all file paths in `config/projects.yaml` are absolute paths

### Issue: Meta-analysis not running

**Solution**: Verify `meta_analysis.enabled: true` and `mode: "multi"`

### Issue: Memory errors during meta-analysis

**Solution**: Reduce number of projects or increase `min_total_reads` filter

### Issue: Different results between single and multi mode

**Solution**: Check for project-specific parameter overrides

## Switching Modes

### From Single to Multi

1. Keep existing `config/config.yaml`
2. Create `config/projects.yaml` with `mode: "multi"`
3. Add your project(s) to `projects:` section

### From Multi to Single

Option 1: Set `mode: "single"` in `config/projects.yaml`
Option 2: Remove/rename `config/projects.yaml`

## Performance Considerations

### Parallel Processing

Projects are processed in parallel. Optimize with:
```bash
snakemake --cores 16  # Total cores across all projects
```

### Disk Space

Multi-project mode requires:
- Per-project space: ~5-10x input data size
- Meta-analysis: Additional ~1-2 GB

### Runtime

- Single project: 2-8 hours (depends on data size)
- Multi-project (3 projects): ~3-10 hours (parallel processing)

## Advanced Features

### Conditional Meta-Analysis

Enable only specific comparisons:
```yaml
meta_analysis:
  enabled: true
  comparisons:
    combined_heatmap: true
    diversity_metrics: true
    contamination_report: false  # Skip this
    taxa_overlap: false           # Skip this
```

### Custom Filtering

Fine-tune what appears in meta-analysis:
```yaml
meta_analysis:
  filters:
    min_projects: 3          # Taxon must be in 3+ projects
    min_total_reads: 100     # At least 100 reads total
    min_identity: 0.98       # Very high confidence only
```

## Example Workflow

Complete workflow for 3 projects:

```bash
# 1. Setup
cp config/projects.example.yaml config/projects.yaml

# 2. Edit configuration
nano config/projects.yaml
# - Set mode: "multi"
# - Add your 3 projects
# - Configure meta_analysis settings

# 3. Validate
snakemake -n  # Dry run

# 4. Run pipeline
snakemake --cores 12

# 5. Check outputs
ls meta_analysis/
# - combined_taxa_heatmap.pdf
# - diversity_metrics.csv
# - contamination_report.csv
# - taxa_overlap_venn.pdf

# 6. Analyze results
Rscript -e "read.csv('meta_analysis/diversity_metrics.csv')"
```

## Support

For issues or questions:
1. Check `logs/` directory for error messages
2. Verify configuration files with: `python3 -c "import yaml; yaml.safe_load(open('config/projects.yaml'))"`
3. Report issues with example configuration
