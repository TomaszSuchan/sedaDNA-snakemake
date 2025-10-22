# Meta-analysis rules for cross-project comparisons
# These rules aggregate and compare results across multiple projects

# Combine all clustered taxa tables from all projects
rule combine_all_taxa:
    input:
        lambda wildcards: [
            f"final_tables/{proj}/{proj}-clustered_taxa_table.csv"
            for proj in PROJECTS
        ]
    output:
        "meta_analysis/all_projects_taxa.csv"
    run:
        import pandas as pd

        all_data = []
        for proj, input_file in zip(PROJECTS, input):
            df = pd.read_csv(input_file)
            df['project'] = proj
            all_data.append(df)

        combined = pd.concat(all_data, ignore_index=True)
        combined.to_csv(output[0], index=False)


# Generate combined heatmap across all projects
rule plot_combined_heatmap:
    input:
        "meta_analysis/all_projects_taxa.csv"
    output:
        linear="meta_analysis/combined_taxa_heatmap.pdf",
        log="meta_analysis/combined_taxa_heatmap_log.pdf"
    params:
        min_projects=META_ANALYSIS["filters"]["min_projects"],
        min_total_reads=META_ANALYSIS["filters"]["min_total_reads"],
        min_identity=META_ANALYSIS["filters"]["min_identity"],
        top_n_taxa=META_ANALYSIS["plotting"].get("top_n_taxa", 100),
        color_scheme=META_ANALYSIS["plotting"].get("color_scheme", "Set2"),
        width=META_ANALYSIS["plotting"].get("width", 15),
        height=META_ANALYSIS["plotting"].get("height", 20),
        cluster_projects=str(META_ANALYSIS["plotting"].get("cluster_projects", True)).upper(),
        cluster_taxa=str(META_ANALYSIS["plotting"].get("cluster_taxa", True)).upper()
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_combined_heatmap.R"


# Generate taxa overlap Venn diagram
rule taxa_overlap_venn:
    input:
        "meta_analysis/all_projects_taxa.csv"
    output:
        "meta_analysis/taxa_overlap_venn.pdf"
    params:
        min_identity=META_ANALYSIS["filters"]["min_identity"],
        min_reads=META_ANALYSIS["filters"]["min_total_reads"]
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/taxa_overlap.R"


# Calculate diversity metrics across projects
rule diversity_metrics:
    input:
        "meta_analysis/all_projects_taxa.csv"
    output:
        "meta_analysis/diversity_metrics.csv"
    params:
        min_identity=META_ANALYSIS["filters"]["min_identity"]
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/diversity_metrics.R"


# Generate cross-project contamination report
rule contamination_report:
    input:
        lambda wildcards: [
            f"final_tables/{proj}/{proj}-combined_classification_table.csv"
            for proj in PROJECTS
        ]
    output:
        "meta_analysis/contamination_report.csv"
    run:
        import pandas as pd

        contamination_summary = []

        for proj, input_file in zip(PROJECTS, input):
            df = pd.read_csv(input_file)

            # Calculate contamination statistics
            total_sequences = len(df)
            flagged = df['remove'].sum() if 'remove' in df.columns else 0
            in_lb = df['in_LB'].sum() if 'in_LB' in df.columns else 0
            in_ib = df['in_IB'].sum() if 'in_IB' in df.columns else 0
            in_sb = df['in_SB'].sum() if 'in_SB' in df.columns else 0

            # Count blanks by type
            blank_counts = df['blank_type'].value_counts().to_dict() if 'blank_type' in df.columns else {}

            contamination_summary.append({
                'project': proj,
                'total_sequences': total_sequences,
                'flagged_sequences': flagged,
                'flagged_percent': 100 * flagged / total_sequences if total_sequences > 0 else 0,
                'in_library_blank': in_lb,
                'in_isolation_blank': in_ib,
                'in_sampling_blank': in_sb,
                'library_blanks': blank_counts.get('LB', 0),
                'pcr_blanks': blank_counts.get('PB', 0),
                'isolation_blanks': blank_counts.get('IB', 0),
                'sampling_blanks': blank_counts.get('SB', 0),
                'samples': blank_counts.get('SAMPLE', 0)
            })

        report = pd.DataFrame(contamination_summary)
        report.to_csv(output[0], index=False)
