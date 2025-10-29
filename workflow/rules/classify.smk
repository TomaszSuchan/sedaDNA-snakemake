rule download_ncbitaxo:
    output:
        "data/ncbitaxo.tgz"
    shell:
        """
        obitaxonomy --download-ncbi --out data/ncbitaxo.tgz
        """

rule merge_all_libraries:
    input:
        lambda wildcards: expand("results-sequences/{project}/{library}.demux.uniq.filtered.denoised.fasta.gz",
                                project=wildcards.project,
                                library=PROJECT_LIBRARIES[wildcards.project])
    output:
        "results-sequences/{project}/{project}-merged.fasta.gz"
    shell:
        """
        cat {input} > {output}
        """

rule classify:
    input:
        fasta = "results-sequences/{project}/{project}-merged.fasta.gz",
        taxonomy = "data/ncbitaxo.tgz"
    output:
        temp("results-classified/{project}/{project}-{db}.classified.fasta")
    params:
        db = lambda wildcards: config["projects"][wildcards.project]["parameters"]["reference_dbs"][wildcards.db]
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obitag"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obitag --max-cpu {threads} \
        -t {input.taxonomy} \
        -R {params.db} \
        {input.fasta} \
        > {output}
        """

rule remove_annotations:
    input:
        "results-classified/{project}/{project}-{db}.classified.fasta"
    output:
        "results-classified/{project}/{project}-{db}.classified.no_annot.fasta"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obiannotate  --max-cpu {threads} \
             --delete-tag=obiclean_head \
             --delete-tag=obiclean_headcount \
             --delete-tag=obiclean_internalcount \
             --delete-tag=obiclean_samplecount \
             --delete-tag=obiclean_singletoncount \
             {input} | \
             obiannotate --number | \
             obiannotate --set-id 'sprintf("seq%04d",annotations.seq_number)' \
             > {output}
        """

rule export_motu_tables:
    input:
        "results-classified/{project}/{project}-{db}.classified.no_annot.fasta"
    output:
        "results-classified/{project}/{project}-{db}.motu_table.csv"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obimatrix --max-cpu {threads} \
        --map obiclean_weight \
          {input} \
          > {output}
        """

rule export_classification_tables:
    input:
        "results-classified/{project}/{project}-{db}.classified.no_annot.fasta"
    output:
        "results-classified/{project}/{project}-{db}.classification_table.csv"
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obicsv --max-cpu {threads} \
        --auto -i -s \
        {input} \
        > {output}
        """