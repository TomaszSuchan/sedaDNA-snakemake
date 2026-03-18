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
        "results/{project}/validation/{library}.barcode_validation.txt"
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

# Project-level sanity checks after all barcode files are validated
rule validate_samples:
    input:
        barcode_files=lambda wildcards: [
            config["projects"][wildcards.project]["libraries"][lib]["barcode_file"]
            for lib in PROJECT_LIBRARIES[wildcards.project]
        ],
        validations=lambda wildcards: expand(
            "results/{project}/validation/{library}.barcode_validation.txt",
            project=wildcards.project,
            library=PROJECT_LIBRARIES[wildcards.project]
        )
    output:
        all_replicates_table="results/{project}/validation/{project}.all_replicates.tsv",
        sample_replicates_table="results/{project}/validation/{project}.sample_replicates.tsv",
        lb_pb_library_table="results/{project}/validation/{project}.LB_PB_replicttes.tsv",
        sb_sampling_table="results/{project}/validation/{project}.SB_replicates.tsv",
        ib_isolation_table="results/{project}/validation/{project}.IB_replicates.tsv",
        duplicate_identity_table="results/{project}/validation/{project}.duplicate_sample_identity.tsv",
        sample_name_errors_table="results/{project}/validation/{project}.sample_name_errors.tsv"
    run:
        import re

        def detect_blank_type(sample_name):
            if sample_name.startswith("LB_"):
                return "LB"
            if sample_name.startswith("PB_"):
                return "PB"
            if sample_name.startswith("IB_"):
                return "IB"
            if "_SB_" in sample_name:
                return "SB"
            return "SAMPLE"

        def extract_replicate(sample_name):
            match = re.search(r"_(\d+)$", sample_name)
            return int(match.group(1)) if match else None

        def sample_no_replicate(sample_name):
            return re.sub(r"_\d+$", "", sample_name)

        def parse_sample_name(sample_name, library_config):
            parts = str(sample_name).split("_")
            blank_type = detect_blank_type(sample_name)
            errors = []
            core = None
            depth = None
            sampling_batch = None
            isolation_batch = None
            library = None
            replicate = None
            sampling_id = None

            if blank_type in ("LB", "PB"):
                if len(parts) != 3:
                    errors.append("expected_format_LB_or_PB")
                else:
                    _, library, rep = parts
                    if rep.isdigit():
                        replicate = int(rep)
                    else:
                        errors.append("non_numeric_replicate")
            elif blank_type == "IB":
                if len(parts) != 4:
                    errors.append("expected_format_IB")
                else:
                    _, isolation_batch, library, rep = parts
                    if rep.isdigit():
                        replicate = int(rep)
                    else:
                        errors.append("non_numeric_replicate")
            elif blank_type == "SB":
                if len(parts) != 6 or parts[1] != "SB":
                    errors.append("expected_format_SB")
                else:
                    core, _, depth, isolation_batch, library, rep = parts
                    sampling_batch = "SB"
                    sampling_id = f"{core}_SB_{depth}_{isolation_batch}"
                    if rep.isdigit():
                        replicate = int(rep)
                    else:
                        errors.append("non_numeric_replicate")
            else:
                if len(parts) != 6:
                    errors.append("expected_format_SAMPLE")
                else:
                    core, depth, sampling_batch, isolation_batch, library, rep = parts
                    if rep.isdigit():
                        replicate = int(rep)
                    else:
                        errors.append("non_numeric_replicate")

            return pd.Series({
                "blank_type": blank_type,
                "core": core,
                "depth": depth,
                "sampling_batch": sampling_batch,
                "isolation_batch": isolation_batch,
                "library_from_sample": library,
                "replicate": replicate,
                "sampling_id": sampling_id,
                "sample_name_error": ";".join(errors)
            })

        rows = []
        for lib in PROJECT_LIBRARIES[wildcards.project]:
            barcode_file = config["projects"][wildcards.project]["libraries"][lib]["barcode_file"]
            df = pd.read_csv(barcode_file)
            df["library_config"] = lib
            rows.append(df)

        combined = pd.concat(rows, ignore_index=True)
        parsed = combined.apply(
            lambda row: parse_sample_name(row["sample"], row["library_config"]),
            axis=1
        )
        combined = pd.concat([combined, parsed], axis=1)
        combined["replicate"] = combined["replicate"].fillna(combined["sample"].map(extract_replicate))
        combined["sample_no_replicate"] = combined["sample"].map(sample_no_replicate)
        combined_valid = combined[combined["sample_name_error"] == ""].copy()

        combined_valid["sampling_unit"] = combined_valid.apply(
            lambda row: (
                f"{row['core']}_{row['sampling_batch']}"
                if row["blank_type"] == "SAMPLE"
                else (
                    f"{row['core']}_{row['depth']}"
                    if row["blank_type"] == "SB"
                    else None
                )
            ),
            axis=1
        )

        all_replicates_table = (
            combined_valid
            .groupby(["sample_no_replicate"], dropna=False)
            .agg(
                n_replicates=("replicate", "nunique")
            )
            .reset_index()
            .sort_values(["sample_no_replicate"])
            .rename(columns={"sample_no_replicate": "sample"})
        )
        all_replicates_table.to_csv(output.all_replicates_table, sep="\t", index=False)

        sample_replicates = (
            combined_valid[combined_valid["blank_type"] == "SAMPLE"]
            .assign(sample_no_library_replicate=lambda df: df["sample"].str.replace(r"_[^_]+_\\d+$", "", regex=True))
            .groupby(["sample_no_library_replicate"], dropna=False)
            .agg(n_replicates=("sample", "size"))
            .reset_index()
            .sort_values(["sample_no_library_replicate"])
            .rename(columns={"sample_no_library_replicate": "sample"})
        )
        sample_replicates.to_csv(output.sample_replicates_table, sep="\t", index=False)

        libraries = (
            combined_valid[["library_from_sample"]]
            .dropna()
            .drop_duplicates()
            .rename(columns={"library_from_sample": "library"})
            .sort_values("library")
        )

        lb_counts = (
            combined_valid[combined_valid["blank_type"] == "LB"]
            .groupby("library_from_sample", dropna=False)
            .agg(n_LB_replicates=("replicate", "nunique"))
            .reset_index()
            .rename(columns={"library_from_sample": "library"})
        )
        pb_counts = (
            combined_valid[combined_valid["blank_type"] == "PB"]
            .groupby("library_from_sample", dropna=False)
            .agg(n_PB_replicates=("replicate", "nunique"))
            .reset_index()
            .rename(columns={"library_from_sample": "library"})
        )
        lb_pb_per_library = (
            libraries
            .merge(lb_counts, on="library", how="left")
            .merge(pb_counts, on="library", how="left")
            .fillna(0)
        )
        lb_pb_per_library["n_LB_replicates"] = lb_pb_per_library["n_LB_replicates"].astype(int)
        lb_pb_per_library["n_PB_replicates"] = lb_pb_per_library["n_PB_replicates"].astype(int)
        lb_pb_per_library.to_csv(output.lb_pb_library_table, sep="\t", index=False)

        sampling_base = (
            combined_valid[combined_valid["blank_type"].isin(["SAMPLE", "SB"])]
            .loc[:, ["sampling_unit"]]
            .dropna()
            .drop_duplicates()
            .rename(columns={"sampling_unit": "sampling"})
        )
        sb_counts = (
            combined_valid[combined_valid["blank_type"] == "SB"]
            .groupby(["sampling_unit"], dropna=False)
            .agg(n_SB_replicates=("sample", "size"))
            .reset_index()
            .rename(columns={"sampling_unit": "sampling"})
        )
        sb_per_sampling = sampling_base.merge(sb_counts, on=["sampling"], how="left").fillna(0).sort_values(["sampling"])
        sb_per_sampling["n_SB_replicates"] = sb_per_sampling["n_SB_replicates"].astype(int)
        sb_per_sampling = sb_per_sampling[["sampling", "n_SB_replicates"]]
        sb_per_sampling.to_csv(output.sb_sampling_table, sep="\t", index=False)

        isolation_base = (
            combined_valid[combined_valid["blank_type"].isin(["SAMPLE", "SB", "IB"])]
            .loc[:, ["isolation_batch"]]
            .dropna()
            .drop_duplicates()
            .rename(columns={"isolation_batch": "isolation"})
        )
        ib_counts = (
            combined_valid[combined_valid["blank_type"] == "IB"]
            .groupby(["isolation_batch"], dropna=False)
            .agg(n_IB_replicates=("sample", "size"))
            .reset_index()
            .rename(columns={"isolation_batch": "isolation"})
        )
        ib_per_isolation = (
            isolation_base
            .merge(ib_counts, on=["isolation"], how="left")
            .fillna(0)
            .sort_values(["isolation"])
        )
        ib_per_isolation["n_IB_replicates"] = ib_per_isolation["n_IB_replicates"].astype(int)
        ib_per_isolation = ib_per_isolation[["isolation", "n_IB_replicates"]]
        ib_per_isolation.to_csv(output.ib_isolation_table, sep="\t", index=False)

        duplicate_identity = (
            combined_valid
            .groupby(
                [
                    "blank_type",
                    "core",
                    "depth",
                    "sampling_batch",
                    "isolation_batch",
                    "library_from_sample",
                    "replicate"
                ],
                dropna=False
            )
            .agg(
                n_rows=("sample", "size"),
                sample_names=("sample", lambda x: ";".join(sorted(set([str(v) for v in x])))),
                sample_tags=("sample_tag", lambda x: ";".join(sorted(set([str(v) for v in x]))))
            )
            .reset_index()
        )
        duplicate_identity = duplicate_identity[duplicate_identity["n_rows"] > 1].sort_values(
            ["library_from_sample", "blank_type", "isolation_batch", "core", "depth", "replicate"]
        )
        duplicate_identity.to_csv(output.duplicate_identity_table, sep="\t", index=False)

        sample_name_errors = combined[combined["sample_name_error"] != ""][
            ["library_config", "sample", "sample_tag", "sample_name_error"]
        ].copy()
        sample_name_errors = sample_name_errors.sort_values(["library_config", "sample"])
        sample_name_errors.to_csv(output.sample_name_errors_table, sep="\t", index=False)

        # Stop the workflow if critical naming/identity issues are detected.
        duplicate_count = len(duplicate_identity)
        name_error_count = len(sample_name_errors)
        if duplicate_count > 0 or name_error_count > 0:
            raise ValueError(
                "validate_samples failed: "
                f"{duplicate_count} duplicate sample identities found "
                f"(see {output.duplicate_identity_table}), "
                f"{name_error_count} sample name errors found "
                f"(see {output.sample_name_errors_table}). "
                "Fix barcode sample names and rerun."
            )

# Split barcode files by length (dynamic)
rule split_barcodes:
    input:
        barcodes=lambda wildcards: config["projects"][wildcards.project]["libraries"][wildcards.library]["barcode_file"],
        validation="results/{project}/validation/{library}.barcode_validation.txt",
        sample_validation="results/{project}/validation/{project}.all_replicates.tsv"
    output:
        "results/{project}/sequences/barcodes-{library}_{length}bp_only.txt"
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
        temp("results/{project}/sequences/{library}.paired.fastq.gz")
    params:
        gap_penalty = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("gap-penalty", 2.0),
        min_identity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("min-identity", 0.9),
        min_overlap = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("min-overlap", 20),
        penalty_scale = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("penalty-scale", 1.0),
    log:
        "logs/{project}/{library}.paired.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obipairing"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
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
        fastq="results/{project}/sequences/{library}.paired.fastq.gz",
        barcodes="results/{project}/sequences/barcodes-{library}_{length}bp_only.txt"
    output:
        temp("results/{project}/sequences/{library}.demux_{length}bp.fastq.gz")
    log:
        "logs/{project}/{library}.demux_{length}bp.log"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obimultiplex"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
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
    return [f"results/{wildcards.project}/sequences/{wildcards.library}.demux_{length}bp.fastq.gz" for length in lengths]

rule concat_barcodes:
    input:
        get_demux_inputs
    output:
        temp("results/{project}/sequences/{library}.demux.fastq.gz")
    shell:
        """
        cat {input} > {output}
        """