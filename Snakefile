# Snakemake workflow for OBI-tools metabarcoding pipeline

# Load configuration
configfile: "config/config.yaml"

# Extract configuration variables
PROJECT = config["project"]
LIBRARIES = list(config["libraries"].keys())

# Include rules
include: "workflow/rules/demultiplex.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/stats.smk"
#include: "workflow/rules/classify.smk"

# Final outputs
rule all:
    input:
        # Validation reports
        expand("results/{project}/{library}.barcode_validation.txt",
               project=PROJECT,
               library=LIBRARIES),
        # Demultiplexed files
        expand("results/{project}/{library}.demux.fastq.gz",
               project=PROJECT,
               library=LIBRARIES),
        # Pairing, dereplication stats
        expand("stats/{project}/{library}.merged_stats.tsv",
               project=PROJECT,
               library=LIBRARIES),
        # Filtered and denoised sequences
        expand("results/{project}/{library}.demux.uniq.filtered.denoised.fasta.gz",
               project=PROJECT,
               library=LIBRARIES),
        # Taxonomic classification
        data/ncbitaxo.tgz,
        expand("results/{project}/{library}.classif",
               project=PROJECT)