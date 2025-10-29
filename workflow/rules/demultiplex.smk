import pandas as pd
from collections import defaultdict
import os

# Function to get barcode lengths for a specific library
def get_library_barcode_lengths(project, library):
    """Get barcode lengths present in a specific library"""
    barcode_file = config["projects"][project]["libraries"][library]["barcode_file"]
    if os.path.exists(barcode_file):
        df = pd.read_csv(barcode_file)
        return sorted(df['sample_tag'].str.len().unique().tolist())
    return []

# Create a dictionary mapping libraries to their barcode lengths
LIBRARY_BARCODE_LENGTHS = {
    project: {
        lib: get_library_barcode_lengths(project, lib) 
        for lib in PROJECT_LIBRARIES[project]
    }
    for project in PROJECTS
}

# Validate barcode files before processing
rule validate_barcodes:
    input:
        lambda wildcards: config["projects"][wildcards.project]["libraries"][wildcards.library]["barcode_file"]
    output:
        "results-sequences/{project}/{library}.barcode_validation.txt"
    run:
        try:
            df = pd.read_csv(input[0])
        except Exception as e:
            raise ValueError(f"Failed to read barcode file '{input[0]}': {e}")
        
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
        barcodes=lambda wildcards: config["projects"][wildcards.project]["libraries"][wildcards.library]["barcode_file"],
        validation="results-sequences/{project}/{library}.barcode_validation.txt"
    output:
        "results-sequences/{project}/barcodes-{library}_{length}bp_only.txt"
    params:
        matching=lambda wildcards: config["projects"][wildcards.project]["parameters"]["obimultiplex"]["matching"],
        primer_mismatches=lambda wildcards: config["projects"][wildcards.project]["parameters"]["obimultiplex"]["primer_mismatches"],
        indels=lambda wildcards: str(config["projects"][wildcards.project]["parameters"]["obimultiplex"]["indels"]).lower()
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
        fvd=lambda wildcards: config["projects"][wildcards.project]["libraries"][wildcards.library]["forward"],
        rev=lambda wildcards: config["projects"][wildcards.project]["libraries"][wildcards.library]["reverse"]
    output:
        temp("results-sequences/{project}/{library}.paired.fastq.gz")
    params:
        gap_penalty = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("gap-penalty", 2.0),
        min_identity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("min-identity", 0.9),
        min_overlap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("min-overlap", 20),
        penalty_scale = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("penalty-scale", 1.0),
    log:
        "logs/{project}/{library}.paired.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("time", 60)
    shell:
        """
        obipairing \
        -F {input.fvd} -R {input.rev} \
        --gap-penalty {params.gap_penalty} \
        --min-identity {params.min_identity} \
        --min-overlap {params.min_overlap} \
        --penalty-scale {params.penalty_scale} \
        --max-cpu {threads} 2> {log} | \
        obigrep \
        -p 'annotations["mode"]=="alignment"' \
        --compress > {output} 2>> {log}
        """

# Demultiplex with different barcode lengths
rule demultiplex:
    input:
        fastq="results-sequences/{project}/{library}.paired.fastq.gz",
        barcodes="results-sequences/{project}/barcodes-{library}_{length}bp_only.txt"
    output:
        temp("results-sequences/{project}/{library}.demux_{length}bp.fastq.gz")
    log:
        "logs/{project}/{library}.demux_{length}bp.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime=lambda wildcards: config["projects"][wildcards.project]["parameters"]["obimultiplex"].get("time", 60)
    shell:
        """
        obimultiplex --tag-list {input.barcodes} \
        --max-cpu {threads} \
        --compress \
        {input.fastq} > {output} 2> {log}
        """

# Concatenate all demultiplexed files per library
def get_demux_inputs(wildcards):
    """Get all demux files for a library based on its barcode lengths"""
    lengths = LIBRARY_BARCODE_LENGTHS[wildcards.project][wildcards.library]
    return [f"results-sequences/{wildcards.project}/{wildcards.library}.demux_{length}bp.fastq.gz" for length in lengths]

rule concat_barcodes:
    input:
        get_demux_inputs
    output:
        temp("results-sequences/{project}/{library}.demux.fastq.gz")
    shell:
        """
        cat {input} > {output}
        """