# Snakemake workflow for OBI-tools metabarcoding pipeline

# Load configuration
configfile: "config/config.yaml"

# Extract configuration variables
PROJECTS = list(config["projects"].keys())
PROJECT_LIBRARIES = {
    project: list(config["projects"][project]["libraries"].keys())
    for project in PROJECTS
}
PROJECT_DBS = {
    project: list(config["projects"][project]["parameters"]["reference_dbs"].keys())
    for project in PROJECTS
}

# Include rules
include: "workflow/rules/demultiplex.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/stats.smk"
include: "workflow/rules/classify.smk"
include: "workflow/rules/analyze.smk"

# Final outputs
rule all:
    input:
        # Validation reports - all projects
        [expand("results-sequences/{project}/{library}.barcode_validation.txt",
                project=project,
                library=PROJECT_LIBRARIES[project])
         for project in PROJECTS],
        # Barcode sanity reports - all projects
        [expand([
                "results-sequences/{project}/{project}.sample_replicates.tsv",
                "results-sequences/{project}/{project}.pb_ib_per_library.tsv",
                "results-sequences/{project}/{project}.sb_per_sampling.tsv",
                "results-sequences/{project}/{project}.ib_per_isolation.tsv"
            ],
            project=project)
         for project in PROJECTS],
        # Demultiplexed files - all projects
        [expand("results-sequences/{project}/{library}.demux.fastq.gz",
                project=project,
                library=PROJECT_LIBRARIES[project])
         for project in PROJECTS],
        # Pairing, dereplication stats - all projects
        [expand("stats/{project}/{library}.merged_stats.tsv",
                project=project,
                library=PROJECT_LIBRARIES[project])
         for project in PROJECTS],
        # Classified fasta files for each database
        [expand("results-classified/{project}/{project}-{db}.classified.no_annot.fasta",
                project=project,
                db=PROJECT_DBS[project])
         for project in PROJECTS],
        # MOTU tables - all projects
        [expand("results-classified/{project}/{project}-{db}.motu_table.csv",
                project=project,
                db=PROJECT_DBS[project])
         for project in PROJECTS],
        # Classification tables - all projects
        [expand("results-classified/{project}/{project}-{db}.classification_table.csv",
                project=project,
                db=PROJECT_DBS[project])
         for project in PROJECTS],
        # Final processed MOTU table - all projects, all databases
        [expand("results-tables/{project}/{project}-{db}-classification_table.csv",
                project=project,
                db=PROJECT_DBS[project])
         for project in PROJECTS],
        # Final clustered taxa table - all projects, all databases
        [expand("results-tables/{project}/{project}-{db}-clustered_taxa_table.csv",
                project=project,
                db=PROJECT_DBS[project])
         for project in PROJECTS]