# Snakemake workflow for OBI-tools metabarcoding pipeline

import os
import yaml

# Determine configuration mode
PROJECTS_CONFIG = "config/projects.yaml"
SINGLE_CONFIG = "config/config.yaml"

# Check if multi-project mode is enabled
if os.path.exists(PROJECTS_CONFIG):
    with open(PROJECTS_CONFIG) as f:
        projects_config = yaml.safe_load(f)

    if projects_config.get("mode") == "multi":
        # Multi-project mode
        MULTI_PROJECT = True

        # Load global settings
        config = projects_config["global"]

        # Extract all projects and their libraries
        PROJECTS = list(projects_config["projects"].keys())
        PROJECT_LIBRARIES = {
            proj: list(projects_config["projects"][proj]["libraries"].keys())
            for proj in PROJECTS
        }

        # Merge project-specific parameters with global parameters
        for proj in PROJECTS:
            proj_params = projects_config["projects"][proj].get("parameters", {})
            merged_params = config["parameters"].copy()
            for key, value in proj_params.items():
                if isinstance(value, dict):
                    merged_params.setdefault(key, {}).update(value)
                else:
                    merged_params[key] = value
            projects_config["projects"][proj]["_merged_params"] = merged_params

        # Store for access in rules
        PROJECTS_DATA = projects_config["projects"]
        META_ANALYSIS = projects_config.get("meta_analysis", {})

    else:
        # Single project mode
        MULTI_PROJECT = False
        configfile: SINGLE_CONFIG
        PROJECTS = [config["project"]]
        PROJECT_LIBRARIES = {config["project"]: list(config["libraries"].keys())}
        PROJECTS_DATA = None
        META_ANALYSIS = None
else:
    # Single project mode (backward compatibility)
    MULTI_PROJECT = False
    configfile: SINGLE_CONFIG
    PROJECTS = [config["project"]]
    PROJECT_LIBRARIES = {config["project"]: list(config["libraries"].keys())}
    PROJECTS_DATA = None
    META_ANALYSIS = None

# For backward compatibility with existing rules
PROJECT = PROJECTS[0] if not MULTI_PROJECT else None
LIBRARIES = PROJECT_LIBRARIES[PROJECTS[0]] if not MULTI_PROJECT else None

# Include rules
include: "workflow/rules/demultiplex.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/stats.smk"
include: "workflow/rules/classify.smk"
include: "workflow/rules/analyze.smk"

# Include meta-analysis rules if in multi-project mode
if MULTI_PROJECT and META_ANALYSIS.get("enabled", False):
    include: "workflow/rules/meta_analysis.smk"

# Helper function to generate project-library combinations
def get_project_library_combinations():
    """Generate all project-library combinations for expand()"""
    combinations = []
    for proj in PROJECTS:
        for lib in PROJECT_LIBRARIES[proj]:
            combinations.append({"project": proj, "library": lib})
    return combinations

# Final outputs
rule all:
    input:
        # Per-project outputs
        # Validation reports
        lambda wildcards: [
            f"results/{combo['project']}/{combo['library']}.barcode_validation.txt"
            for combo in get_project_library_combinations()
        ],
        # Demultiplexed files
        lambda wildcards: [
            f"results/{combo['project']}/{combo['library']}.demux.fastq.gz"
            for combo in get_project_library_combinations()
        ],
        # Pairing, dereplication stats
        lambda wildcards: [
            f"stats/{combo['project']}/{combo['library']}.merged_stats.tsv"
            for combo in get_project_library_combinations()
        ],
        # Classified fasta files for each database
        lambda wildcards: [
            f"results/{proj}/{proj}-{db}.classified.fasta"
            for proj in PROJECTS
            for db in config["reference_dbs"].keys()
        ],
        # Fasta files without annotations
        lambda wildcards: [
            f"results/{proj}/{proj}-{db}.classified.no_annot.fasta"
            for proj in PROJECTS
            for db in config["reference_dbs"].keys()
        ],
        # MOTU tables
        lambda wildcards: [
            f"results/{proj}/{proj}-{db}.motu_table.csv"
            for proj in PROJECTS
            for db in config["reference_dbs"].keys()
        ],
        # Classification tables
        lambda wildcards: [
            f"results/{proj}/{proj}-{db}.classification_table.csv"
            for proj in PROJECTS
            for db in config["reference_dbs"].keys()
        ],
        # Final processed MOTU tables
        lambda wildcards: [
            f"final_tables/{proj}/{proj}-combined_classification_table.csv"
            for proj in PROJECTS
        ],
        # Final clustered taxa tables
        lambda wildcards: [
            f"final_tables/{proj}/{proj}-clustered_taxa_table.csv"
            for proj in PROJECTS
        ],
        # Taxa heatmap plots
        lambda wildcards: [
            f"final_plots/{proj}/{proj}-taxa_heatmap.pdf"
            for proj in PROJECTS
        ],
        # Taxa heatmap plots (log)
        lambda wildcards: [
            f"final_plots/{proj}/{proj}-taxa_heatmap_log.pdf"
            for proj in PROJECTS
        ],
        # Meta-analysis outputs (if enabled and multi-project mode)
        lambda wildcards: (
            [
                "meta_analysis/combined_taxa_heatmap.pdf",
                "meta_analysis/combined_taxa_heatmap_log.pdf",
                "meta_analysis/taxa_overlap_venn.pdf",
                "meta_analysis/diversity_metrics.csv",
                "meta_analysis/contamination_report.csv"
            ] if MULTI_PROJECT and META_ANALYSIS.get("enabled", False) else []
        )