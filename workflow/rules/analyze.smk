rule process_motu:
    input:
        motu_tables = expand(
            "results/{project}/{project}-{db}.motu_table.csv",
            project="{project}",
            db=config["reference_dbs"].keys()
        ),
        classification_tables = expand(
            "results/{project}/{project}-{db}.classification_table.csv",
            project="{project}",
            db=config["reference_dbs"].keys()
        )
    output:
        "final_tables/{project}/{project}-combined_classification_table.csv"
    params:
        db_prefixes = ",".join(config["reference_dbs"].keys()),
        reads_within = config["parameters"]["seq_filters"].get("reads_within", 3),
        reads_across = config["parameters"]["seq_filters"].get("reads_across", 10),
        reads_replicates = config["parameters"]["seq_filters"].get("reads_replicates", 3)
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript workflow/scripts/process_motu.R \
            {input.motu_tables[0]} \
            {input.classification_tables} \
            "{params.db_prefixes}" \
            {params.reads_within} \
            {params.reads_across} \
            {params.reads_replicates} \
            {output}
        """

rule cluster_taxa:
    input:
        "final_tables/{project}/{project}-combined_classification_table.csv"
    output:
        "final_tables/{project}/{project}-clustered_taxa_table.csv"
    params:
        min_identity = config["parameters"]["tax_filters"].get("min_identity", 1.0),
        db_prefix = config["parameters"]["tax_filters"]["db_prefix"]
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript workflow/scripts/cluster_taxa.R \
            {input} \
            {params.min_identity} \
            {params.db_prefix} \
            {output}
        """

rule plot_taxa_heatmap_log:
    input:
        "final_tables/{project}/{project}-clustered_taxa_table.csv"
    output:
        "final_plots/{project}/{project}-taxa_heatmap_log.pdf"
    params:
        log_transform = TRUE,
        top_n_taxa = config["parameters"]["plotting"].get("top_n_taxa", 50)
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript workflow/scripts/plot_taxa_heatmap.R \
            {params.log_transform} {params.top_n_taxa} \
            {input} \
            {output}
        """

rule plot_taxa_heatmap:
    input:
        "final_tables/{project}/{project}-clustered_taxa_table.csv"
    output:
        "final_plots/{project}/{project}-taxa_heatmap.pdf"
    params:
        log_transform = FALSE,
        top_n_taxa = config["parameters"]["plotting"].get("top_n_taxa", 50),
        width = config["parameters"]["plotting"].get("width", 10),
        height = config["parameters"]["plotting"].get("height", 8)
    conda:
        "../envs/r.yaml"
    shell:
        """
        Rscript workflow/scripts/plot_taxa_heatmap.R \
            {params.log_transform} {params.top_n_taxa} \
            {input} \
            {output}
        """