import pandas as pd
from collections import defaultdict
import os

# Helper function to get library config in both single and multi-project modes
def get_library_config(project, library, key):
    """Get library configuration value for a specific project and library"""
    if MULTI_PROJECT:
        return PROJECTS_DATA[project]["libraries"][library][key]
    else:
        return config["libraries"][library][key]

# Function to get barcode lengths for a specific library
def get_library_barcode_lengths(project, library):
    """Get barcode lengths present in a specific library"""
    barcode_file = get_library_config(project, library, "barcode_file")
    if os.path.exists(barcode_file):
        df = pd.read_csv(barcode_file)
        return sorted(df['sample_tag'].str.len().unique().tolist())
    return []

# Create a dictionary mapping (project, library) tuples to their barcode lengths
if MULTI_PROJECT:
    LIBRARY_BARCODE_LENGTHS = {
        (proj, lib): get_library_barcode_lengths(proj, lib)
        for proj in PROJECTS
        for lib in PROJECT_LIBRARIES[proj]
    }
else:
    LIBRARY_BARCODE_LENGTHS = {
        (PROJECT, lib): get_library_barcode_lengths(PROJECT, lib)
        for lib in LIBRARIES
    }

# Validate barcode files before processing
rule validate_barcodes:
    input:
        lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, "barcode_file")
    output:
        "results/{PROJECT}/{library}.barcode_validation.txt"
    run:
        df = pd.read_csv(input[0])
        
        # Validation checks
        required_columns = ['experiment', 'sample', 'sample_tag', 'forward_primer', 'reverse_primer']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing required columns in {input[0]}: {missing_columns}")
        
        # Check for duplicate barcodes
        duplicates = df[df.duplicated('sample_tag', keep=False)]
        if not duplicates.empty:
            print(f"Warning: Duplicate barcodes found in {input[0]}:")
            print(duplicates[['sample', 'sample_tag']])
        
        # Check barcode lengths
        barcode_lengths = df['sample_tag'].str.len().value_counts().sort_index()
        print(f"Barcode length distribution for {wildcards.library}:")
        for length, count in barcode_lengths.items():
            print(f"  {int(length)}bp: {int(count)} barcodes")
        
        # Check for primer consistency
        unique_forward = df['forward_primer'].nunique()
        unique_reverse = df['reverse_primer'].nunique()
        
        # Write validation results
        with open(output[0], 'w') as f:
            f.write(f"Validation completed for {wildcards.library}\n")
            f.write(f"Total barcodes: {len(df)}\n")
            # Convert dict keys and values to regular Python ints
            f.write(f"Barcode lengths: {dict((int(k), int(v)) for k, v in barcode_lengths.items())}\n")
            f.write(f"Unique forward primers: {unique_forward}\n")
            f.write(f"Unique reverse primers: {unique_reverse}\n")
            if duplicates.empty:
                f.write("No duplicate barcodes found\n")
            else:
                f.write(f"Warning: {len(duplicates)} duplicate barcodes found\n")

# Split barcode files by length (dynamic)
rule split_barcodes:
    input:
        barcodes=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, "barcode_file"),
        validation="results/{PROJECT}/{library}.barcode_validation.txt"
    output:
        "results/{PROJECT}/barcodes-{library}_{length}bp_only.txt"
    params:
        matching=config["barcodes"]["matching"],
        primer_mismatches=config["barcodes"]["primer_mismatches"],
        indels=str(config["barcodes"]["indels"]).lower()
    run:
        # Read the barcode CSV file
        df = pd.read_csv(input.barcodes)
        
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
        
        # Write barcode file for this length
        with open(output[0], 'w') as f:
            f.write('\n'.join(header_lines) + '\n')
            for _, row in df_filtered.iterrows():
                f.write(f"{row['experiment']},{row['sample']},{row['sample_tag']},{row['forward_primer']},{row['reverse_primer']}\n")
        
        print(f"Created barcode file for {wildcards.library} with {len(df_filtered)} {target_length}bp barcodes")

# Pair reads and keep only merged
rule pair_reads:
    input:
        fvd=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, "forward"),
        rev=lambda wildcards: get_library_config(wildcards.PROJECT, wildcards.library, "reverse")
    output:
        temp("results/{PROJECT}/{library}.paired.fastq.gz")
    params:
        gap_penalty = config["parameters"]["obipairing"].get("gap-penalty", 2.0),
        min_identity = config["parameters"]["obipairing"].get("min-identity", 0.9),
        min_overlap = config["parameters"]["obipairing"].get("min-overlap", 20),
        penalty_scale = config["parameters"]["obipairing"].get("penalty-scale", 1.0),
        max_cpu = config["parameters"].get("max-cpu", 1)
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
        max_cpu=config["parameters"].get("max-cpu", 1)
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
    return [f"results/{wildcards.PROJECT}/{wildcards.library}.demux_{length}bp.fastq.gz" for length in lengths]

rule concat_barcodes:
    input:
        get_demux_inputs
    output:
        temp("results/{PROJECT}/{library}.demux.fastq.gz")
    shell:
        """
        cat {input} > {output}
        """
