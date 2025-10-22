"""
Demultiplexing rules for OBI-tools pipeline

Handles:
- Barcode validation
- Variable-length barcode splitting
- Read pairing
- Demultiplexing by barcode
"""

import pandas as pd
import sys
from pathlib import Path

# Import config loader for constants
sys.path.insert(0, str(Path("workflow/lib").resolve()))
from config_loader import ConfigKeys


def get_library_config(project: str, library: str, key: str):
    """
    Get library configuration value for a specific project and library.

    Args:
        project: Project name
        library: Library name
        key: Configuration key (forward, reverse, barcode_file)

    Returns:
        Configuration value
    """
    return PROJECTS_DATA[project][ConfigKeys.LIBRARIES][library][key]


def get_library_barcode_lengths(project: str, library: str):
    """
    Get barcode lengths present in a specific library.

    Args:
        project: Project name
        library: Library name

    Returns:
        Sorted list of barcode lengths
    """
    barcode_file = get_library_config(project, library, ConfigKeys.BARCODE_FILE)

    if not Path(barcode_file).exists():
        raise FileNotFoundError(
            f"Barcode file not found: {barcode_file}\n"
            f"For project '{project}', library '{library}'"
        )

    try:
        df = pd.read_csv(barcode_file)
        return sorted(df['sample_tag'].str.len().unique().tolist())
    except Exception as e:
        raise ValueError(
            f"Error reading barcode file '{barcode_file}': {e}\n"
            f"For project '{project}', library '{library}'"
        )


# Pre-compute barcode lengths for all project-library combinations
LIBRARY_BARCODE_LENGTHS = {
    (proj, lib): get_library_barcode_lengths(proj, lib)
    for proj in PROJECTS
    for lib in PROJECT_LIBRARIES[proj]
}


# Validate barcode files before processing
rule validate_barcodes:
    input:
        lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, ConfigKeys.BARCODE_FILE)
    output:
        "results/{PROJECT}/{library}.barcode_validation.txt"
    run:
        try:
            df = pd.read_csv(input[0])
        except Exception as e:
            raise ValueError(f"Failed to read barcode file '{input[0]}': {e}")

        # Validation checks
        required_columns = ['experiment', 'sample', 'sample_tag', 'forward_primer', 'reverse_primer']
        missing_columns = [col for col in required_columns if col not in df.columns]

        if missing_columns:
            raise ValueError(
                f"Missing required columns in {input[0]}: {missing_columns}\n"
                f"Required columns: {required_columns}\n"
                f"Found columns: {list(df.columns)}"
            )

        # Check for empty required columns
        for col in required_columns:
            if df[col].isna().any():
                raise ValueError(f"Column '{col}' contains empty values in {input[0]}")

        # Check for duplicate barcodes
        duplicates = df[df.duplicated('sample_tag', keep=False)]
        if not duplicates.empty:
            print(f"⚠️ Warning: Duplicate barcodes found in {input[0]}:")
            print(duplicates[['sample', 'sample_tag']])

        # Check barcode lengths
        barcode_lengths = df['sample_tag'].str.len().value_counts().sort_index()
        print(f"✓ Barcode length distribution for {wildcards.library}:")
        for length, count in barcode_lengths.items():
            print(f"  {int(length)}bp: {int(count)} barcodes")

        # Check for primer consistency
        unique_forward = df['forward_primer'].nunique()
        unique_reverse = df['reverse_primer'].nunique()

        # Write validation results
        with open(output[0], 'w') as f:
            f.write(f"Validation completed for {wildcards.library}\n")
            f.write(f"Project: {wildcards.PROJECT}\n")
            f.write(f"Total barcodes: {len(df)}\n")
            f.write(f"Barcode lengths: {dict((int(k), int(v)) for k, v in barcode_lengths.items())}\n")
            f.write(f"Unique forward primers: {unique_forward}\n")
            f.write(f"Unique reverse primers: {unique_reverse}\n")
            if duplicates.empty:
                f.write("✓ No duplicate barcodes found\n")
            else:
                f.write(f"⚠️  Warning: {len(duplicates)} duplicate barcodes found\n")


# Split barcode files by length (dynamic)
rule split_barcodes:
    input:
        barcodes=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, ConfigKeys.BARCODE_FILE),
        validation="results/{PROJECT}/{library}.barcode_validation.txt"
    output:
        "results/{PROJECT}/barcodes-{library}_{length}bp_only.txt"
    params:
        matching=config[ConfigKeys.BARCODES]["matching"],
        primer_mismatches=config[ConfigKeys.BARCODES]["primer_mismatches"],
        indels=str(config[ConfigKeys.BARCODES]["indels"]).lower()
    run:
        # Read the barcode CSV file
        try:
            df = pd.read_csv(input.barcodes)
        except Exception as e:
            raise ValueError(f"Failed to read barcode file '{input.barcodes}': {e}")

        # Extract the target length from wildcards
        target_length = int(wildcards.length)

        # Create header for OBITools barcode files
        header_lines = [
            f"@param,matching,{params.matching}",
            f"@param,primer_mismatches,{params.primer_mismatches}",
            f"@param,indels,{params.indels}",
            "experiment,sample,sample_tag,forward_primer,reverse_primer"
        ]

        # Filter by barcode length
        df_filtered = df[df['sample_tag'].str.len() == target_length].copy()

        if len(df_filtered) == 0:
            raise ValueError(
                f"No barcodes of length {target_length}bp found in {input.barcodes}\n"
                f"Available lengths: {sorted(df['sample_tag'].str.len().unique())}"
            )

        # Write barcode file for this length
        with open(output[0], 'w') as f:
            f.write('\n'.join(header_lines) + '\n')
            for _, row in df_filtered.iterrows():
                f.write(f"{row['experiment']},{row['sample']},{row['sample_tag']},{row['forward_primer']},{row['reverse_primer']}\n")

        print(f"✓ Created barcode file for {wildcards.library}: {len(df_filtered)} {target_length}bp barcodes")


# Pair reads and keep only merged
rule pair_reads:
    input:
        fvd=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, ConfigKeys.FORWARD),
        rev=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, ConfigKeys.REVERSE)
    output:
        temp("results/{PROJECT}/{library}.paired.fastq.gz")
    params:
        gap_penalty=config[ConfigKeys.PARAMETERS]["obipairing"].get("gap-penalty", 2.0),
        min_identity=config[ConfigKeys.PARAMETERS]["obipairing"].get("min-identity", 0.9),
        min_overlap=config[ConfigKeys.PARAMETERS]["obipairing"].get("min-overlap", 20),
        penalty_scale=config[ConfigKeys.PARAMETERS]["obipairing"].get("penalty-scale", 1.0),
        max_cpu=config[ConfigKeys.PARAMETERS].get("max-cpu", 1)
    log:
        "logs/{PROJECT}/{library}.paired.log"
    shell:
        """
        obipairing \
        -F {input.fvd} -R {input.rev} \
        --gap-penalty {params.gap_penalty} \
        --min-identity {params.min_identity} \
        --min-overlap {params.min_overlap} \
        --penalty-scale {params.penalty_scale} \
        --max-cpu {params.max_cpu} 2> {log} | \
        obigrep \
        -p 'annotations["mode"]=="alignment"' \
        --compress > {output} 2>> {log}
        """


# Demultiplex with different barcode lengths
rule demultiplex:
    input:
        fastq="results/{PROJECT}/{library}.paired.fastq.gz",
        barcodes="results/{PROJECT}/barcodes-{library}_{length}bp_only.txt"
    output:
        "results/{PROJECT}/{library}.demux_{length}bp.fastq.gz"
    params:
        max_cpu=config[ConfigKeys.PARAMETERS].get("max-cpu", 1)
    log:
        "logs/{PROJECT}/{library}.demux_{length}bp.log"
    shell:
        """
        obimultiplex --tag-list {input.barcodes} \
        --max-cpu {params.max_cpu} \
        --compress \
        {input.fastq} > {output} 2> {log}
        """


# Concatenate all demultiplexed files per library
def get_demux_inputs(wildcards):
    """Get all demux files for a library based on its barcode lengths"""
    lengths = LIBRARY_BARCODE_LENGTHS[(wildcards.PROJECT, wildcards.library)]
    return [
        f"results/{wildcards.PROJECT}/{wildcards.library}.demux_{length}bp.fastq.gz"
        for length in lengths
    ]


rule concat_barcodes:
    input:
        get_demux_inputs
    output:
        temp("results/{PROJECT}/{library}.demux.fastq.gz")
    shell:
        """
        cat {input} > {output}
        """
