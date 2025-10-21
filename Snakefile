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
include: "workflow/rules/classify.smk"
include: "workflow/rules/analyze.smk"

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
        # Classified fasta files for each database
        expand("results/{project}/{project}-{db}.classified.fasta",
               project=PROJECT,
               db=config["reference_dbs"].keys()),
        # Fasta files without annotations
        expand("results/{project}/{project}-{db}.classified.no_annot.fasta",
               project=PROJECT,
               db=config["reference_dbs"].keys()),
        # MOTU tables
        expand("results/{project}/{project}-{db}.motu_table.csv",
               project=PROJECT,
               db=config["reference_dbs"].keys()),
        # Classification tables
        expand("results/{project}/{project}-{db}.classification_table.csv",
               project=PROJECT,
               db=config["reference_dbs"].keys()),
        # Final processed MOTU table
        "final_tables/{project}/{project}-combined_classification_table.csv".format(project=PROJECT),
        # Final clustered taxa table
        "final_tables/{project}/{project}-clustered_taxa_table.csv".format(project=PROJECT),
        # Taxa heatmap plot
        "final_plots/{project}/{project}-taxa_heatmap.pdf".format(project=PROJECT)