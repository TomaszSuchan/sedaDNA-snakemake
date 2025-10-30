rule process_motu:
    input:
        motu_table = "results-classified/{project}/{project}-{db}.motu_table.csv",
        classification_table = "results-classified/{project}/{project}-{db}.classification_table.csv",
        demultiplex_stats = "stats/{project}/{project}.demux_stats_combined.json"
    output:
        "results-tables/{project}/{project}-{db}-combined_classification_table.csv"
    log:
        "logs/{project}/{project}-{db}-combined_classification_table.log"
    params:
        reads_within = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_within", 3),
        reads_across = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_across", 10),
        reads_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_replicates", 3)
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/process_motu.R \
            {input.motu_table} \
            {input.classification_table} \
            {input.demultiplex_stats} \
            {wildcards.db} \
            {params.reads_within} \
            {params.reads_across} \
            {params.reads_replicates} \
            {output} 2> {log}
        """

rule cluster_taxa:
    input:
        "results-tables/{project}/{project}-{db}-combined_classification_table.csv"
    output:
        "results-tables/{project}/{project}-{db}-clustered_taxa_table.csv"
    log:
        "logs/{project}/{project}-{db}-clustered_taxa_table.log"
    params:
        min_identity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["tax_filters"].get("min_identity", 1.0)
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/cluster_taxa.R \
            {input} \
            {params.min_identity} \
            {output} 2> {log}
        """

rule plot_taxa_heatmap_log:
    input:
        "results-tables/{project}/{project}-{db}-clustered_taxa_table.csv"
    output:
        "results-plots/{project}/{project}-{db}-taxa_heatmap_log.pdf"
    log:
        "logs/{project}/{project}-{db}-taxa_heatmap_log.log"
    params:
        log_transform = "TRUE",
        top_n_taxa = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("top_n_taxa", 50),
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("height", 8)
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/plot_taxa_heatmap.R \
            {params.log_transform} {params.top_n_taxa} {params.width} {params.height} \
            {input} \
            {output} 2> {log}
        """

rule plot_taxa_heatmap:
    input:
        "results-tables/{project}/{project}-{db}-clustered_taxa_table.csv"
    output:
        "results-plots/{project}/{project}-{db}-taxa_heatmap.pdf"
    log:
        "logs/{project}/{project}-{db}-taxa_heatmap.log"
    params:
        log_transform = "FALSE",
        top_n_taxa = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("top_n_taxa", 50),
        width = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("width", 10),
        height = lambda wildcards: config["projects"][wildcards.project]["parameters"]["plotting"].get("height", 8)
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/plot_taxa_heatmap.R \
            {params.log_transform} {params.top_n_taxa} {params.width} {params.height} \
            {input} \
            {output} 2> {log}
        """