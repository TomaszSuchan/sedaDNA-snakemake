def project_param(wildcards, section, key, default):
    return config["projects"][wildcards.project]["parameters"].get(section, {}).get(key, default)


rule process_motu:
    # One run per project (MergeAndFilter design): the read/replicate filters do
    # not depend on the reference database, so every database gets identical
    # filtering. The MOTU tables of all databases hold the same counts (the
    # sequence ids are assigned before classification), so the first one is used.
    input:
        motu_table = lambda wildcards: "results/{project}/classified/{project}-%s.motu_table.csv"
                                       % PROJECT_DBS[wildcards.project][0],
        classification_tables = lambda wildcards: expand(
            "results/{project}/classified/{project}-{db}.classification_table.csv",
            project=wildcards.project, db=PROJECT_DBS[wildcards.project]),
        demultiplex_stats = "results/{project}/stats/{project}.demux_stats_combined.json"
    output:
        replicates     = "results/{project}/tables/{project}-replicates.tsv.gz",
        pcrs           = "results/{project}/tables/{project}-pcrs.tsv",
        samples        = "results/{project}/tables/{project}-samples.tsv",
        sequences      = "results/{project}/tables/{project}-sequences.tsv",
        sample_summary = "results/{project}/tables/{project}-sample_summary.tsv",
        bio_sample_summary = "results/{project}/tables/{project}-bio_sample_summary.tsv",
        blank_report   = "results/{project}/tables/{project}-blank_report.tsv",
        params_json    = "results/{project}/tables/{project}-filter_params.json"
    log:
        stderr = "logs/{project}/{project}-process_motu.log",
        stdout = "results/{project}/tables/{project}-process_motu_info.txt"
    params:
        prefix = "results/{project}/tables/{project}",
        db_tables = lambda wildcards: " ".join(
            "%s=results/%s/classified/%s-%s.classification_table.csv"
            % (db, wildcards.project, wildcards.project, db)
            for db in PROJECT_DBS[wildcards.project]),
        reads_within = lambda wildcards: project_param(wildcards, "seq_filters", "reads_within", 3),
        reads_across = lambda wildcards: project_param(wildcards, "seq_filters", "reads_across", 10),
        reads_replicates = lambda wildcards: project_param(wildcards, "seq_filters", "reads_replicates", 3),
        min_identity = lambda wildcards: project_param(wildcards, "tax_filters", "min_identity", 1.0)
    conda:
        "../envs/r.yaml"
    resources:
        runtime = lambda wildcards: project_param(wildcards, "process_motu", "time", 240),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/process_motu.R \
            {input.motu_table} \
            {input.demultiplex_stats} \
            {params.reads_within} \
            {params.reads_across} \
            {params.reads_replicates} \
            {params.min_identity} \
            {params.prefix} \
            {params.db_tables} > {log.stdout} 2> {log.stderr}
        """

rule cluster_taxa:
    # one table per database, per sediment sample (bio_sample = core_depth)
    input:
        "results/{project}/tables/{project}-bio_sample_summary.tsv"
    output:
        table = "results/{project}/tables/{project}-{db}-clustered_taxa_table.csv"
    log:
        stderr = "logs/{project}/{project}-{db}-clustered_taxa_table.log",
        stdout = "results/{project}/tables/{project}-{db}-clustered_taxa_info.txt"
    params:
        min_identity = lambda wildcards: project_param(wildcards, "tax_filters", "min_identity", 1.0)
    conda:
        "../envs/r.yaml"
    resources:
        runtime = lambda wildcards: project_param(wildcards, "cluster_taxa", "time", 120),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        Rscript workflow/scripts/cluster_taxa.R \
            {input} \
            {wildcards.db} \
            {params.min_identity} \
            {output.table} > {log.stdout} 2> {log.stderr}
        """
