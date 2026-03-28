rule process_motu:
    input:
        motu_table = "results/{project}/classified/{project}-{db}.motu_table.csv",
        classification_table = "results/{project}/classified/{project}-{db}.classification_table.csv",
        demultiplex_stats = "results/{project}/stats/{project}.demux_stats_combined.json"
    output:
        table = "results/{project}/tables/{project}-{db}-classification_table.csv"
    log:
        stderr = "logs/{project}/{project}-{db}-classification_table.log",
        stdout = "results/{project}/tables/{project}-{db}-classification_info.txt"
    params:
        reads_within = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_within", 3),
        reads_across = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_across", 10),
        reads_replicates = lambda wildcards: config["projects"][wildcards.project]["parameters"]["seq_filters"].get("reads_replicates", 3)
    conda:
        "../envs/r.yaml"
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["process_motu"].get("time", 240),
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
            {output.table} > {log.stdout} 2> {log.stderr}
        """

rule cluster_taxa:
    input:
        "results/{project}/tables/{project}-{db}-classification_table.csv"
    output:
        table = "results/{project}/tables/{project}-{db}-clustered_taxa_table.csv"
    log:
        stderr = "logs/{project}/{project}-{db}-clustered_taxa_table.log",
        stdout = "results/{project}/tables/{project}-{db}-clustered_taxa_info.txt"
    params:
        min_identity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["tax_filters"].get("min_identity", 1.0)
    conda:
        "../envs/r.yaml"
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["cluster_taxa"].get("time", 120),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/cluster_taxa.R \
            {input} \
            {params.min_identity} \
            {output.table} > {log.stdout} 2> {log.stderr}
        """