"""
Snakemake workflow for OBI-tools metabarcoding pipeline
Multi-project batch processing mode

This workflow processes multiple sedaDNA metabarcoding projects
simultaneously with automated cross-project comparative analysis.
"""

import sys
from pathlib import Path

# Add workflow lib to path
sys.path.insert(0, str(Path("workflow/lib").resolve()))

from config_loader import (
    load_and_validate_config,
    extract_projects_info,
    ConfigurationError,
    SecurityError,
    ConfigKeys
)

# Configuration file path
PROJECTS_CONFIG = "config/projects.yaml"

# Load and validate configuration
try:
    projects_config = load_and_validate_config(
        PROJECTS_CONFIG,
        check_file_existence=False  # Set to True for strict validation
    )
except (ConfigurationError, SecurityError) as e:
    print(f"\n❌ Configuration Error:\n{e}\n", file=sys.stderr)
    sys.exit(1)
except FileNotFoundError as e:
    print(f"\n❌ File Not Found:\n{e}\n", file=sys.stderr)
    sys.exit(1)

# Extract global config
config = projects_config[ConfigKeys.GLOBAL]

# Extract project information
PROJECTS, PROJECT_LIBRARIES, PROJECTS_DATA = extract_projects_info(projects_config)

# Meta-analysis settings
META_ANALYSIS = projects_config.get(ConfigKeys.META_ANALYSIS, {"enabled": False})

# Log configuration
print(f"\n{'='*60}")
print(f"sedaDNA-snakemake Multi-Project Pipeline")
print(f"{'='*60}")
print(f"Projects to process: {len(PROJECTS)}")
for proj in PROJECTS:
    n_libs = len(PROJECT_LIBRARIES[proj])
    print(f"  • {proj}: {n_libs} {'library' if n_libs == 1 else 'libraries'}")
print(f"Reference databases: {len(config[ConfigKeys.REFERENCE_DBS])}")
print(f"Meta-analysis: {'enabled' if META_ANALYSIS.get('enabled') else 'disabled'}")
print(f"{'='*60}\n")

# Include rules
include: "workflow/rules/demultiplex.smk"
include: "workflow/rules/filter.smk"
include: "workflow/rules/stats.smk"
include: "workflow/rules/classify.smk"
include: "workflow/rules/analyze.smk"

# Include meta-analysis rules if enabled
if META_ANALYSIS.get("enabled", False):
    include: "workflow/rules/meta_analysis.smk"


# Helper function to generate project-library combinations
def get_project_library_combinations():
    """Generate all project-library combinations for expand()"""
    combinations = []
    for proj in PROJECTS:
        for lib in PROJECT_LIBRARIES[proj]:
            combinations.append({"project": proj, "library": lib})
    return combinations


# Pre-compute output file lists (performance optimization)
PROJECT_LIBRARY_COMBOS = get_project_library_combinations()

VALIDATION_REPORTS = [
    f"results/{combo['project']}/{combo['library']}.barcode_validation.txt"
    for combo in PROJECT_LIBRARY_COMBOS
]

DEMUX_FILES = [
    f"results/{combo['project']}/{combo['library']}.demux.fastq.gz"
    for combo in PROJECT_LIBRARY_COMBOS
]

STATS_FILES = [
    f"stats/{combo['project']}/{combo['library']}.merged_stats.tsv"
    for combo in PROJECT_LIBRARY_COMBOS
]

CLASSIFIED_FASTA = [
    f"results/{proj}/{proj}-{db}.classified.fasta"
    for proj in PROJECTS
    for db in config[ConfigKeys.REFERENCE_DBS].keys()
]

CLASSIFIED_NO_ANNOT = [
    f"results/{proj}/{proj}-{db}.classified.no_annot.fasta"
    for proj in PROJECTS
    for db in config[ConfigKeys.REFERENCE_DBS].keys()
]

MOTU_TABLES = [
    f"results/{proj}/{proj}-{db}.motu_table.csv"
    for proj in PROJECTS
    for db in config[ConfigKeys.REFERENCE_DBS].keys()
]

CLASSIFICATION_TABLES = [
    f"results/{proj}/{proj}-{db}.classification_table.csv"
    for proj in PROJECTS
    for db in config[ConfigKeys.REFERENCE_DBS].keys()
]

COMBINED_CLASSIFICATION_TABLES = [
    f"final_tables/{proj}/{proj}-combined_classification_table.csv"
    for proj in PROJECTS
]

CLUSTERED_TAXA_TABLES = [
    f"final_tables/{proj}/{proj}-clustered_taxa_table.csv"
    for proj in PROJECTS
]

TAXA_HEATMAPS = [
    f"final_plots/{proj}/{proj}-taxa_heatmap.pdf"
    for proj in PROJECTS
]

TAXA_HEATMAPS_LOG = [
    f"final_plots/{proj}/{proj}-taxa_heatmap_log.pdf"
    for proj in PROJECTS
]

META_ANALYSIS_OUTPUTS = [
    "meta_analysis/combined_taxa_heatmap.pdf",
    "meta_analysis/combined_taxa_heatmap_log.pdf",
    "meta_analysis/taxa_overlap_venn.pdf",
    "meta_analysis/diversity_metrics.csv",
    "meta_analysis/contamination_report.csv"
] if META_ANALYSIS.get("enabled", False) else []


# Final outputs
rule all:
    input:
        VALIDATION_REPORTS,
        DEMUX_FILES,
        STATS_FILES,
        CLASSIFIED_FASTA,
        CLASSIFIED_NO_ANNOT,
        MOTU_TABLES,
        CLASSIFICATION_TABLES,
        COMBINED_CLASSIFICATION_TABLES,
        CLUSTERED_TAXA_TABLES,
        TAXA_HEATMAPS,
        TAXA_HEATMAPS_LOG,
        META_ANALYSIS_OUTPUTS
