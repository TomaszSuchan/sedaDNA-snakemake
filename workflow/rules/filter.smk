rule dereplicate:
    input:
        "results-sequences/{project}/{library}.demux.fastq.gz"
    output:
        temp("results-sequences/{project}/{library}.demux.uniq.fasta.gz")
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obiuniq"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obiuniq --max-cpu {threads} -m sample {input} > {output}
        """

rule filter_annotations:
    input:
        "results-sequences/{project}/{library}.demux.uniq.fasta.gz"
    output:
        temp("results-sequences/{project}/{library}.demux.uniq.counts.fasta.gz")
    params:
        annotation_db = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("annotation-db", ""),
        min_identity = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("min-identity", 0.9)
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obiannotate --max-cpu {threads} -k count -k merged_sample {input} > {output}
        """

rule filter_counts:
    input:
        "results-sequences/{project}/{library}.demux.uniq.counts.fasta.gz"
    output:
        temp("results-sequences/{project}/{library}.demux.uniq.filtered.fasta.gz")
    params:
        min_count = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("min-count", 2),
        min_length = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("min-length", 10)
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obigrep --max-cpu {threads} --min-length {params.min_length} \
        --min-count {params.min_count} {input} > {output}
        """

rule denoise:
    input:
        "results-sequences/{project}/{library}.demux.uniq.filtered.fasta.gz"
    output:
        "results-sequences/{project}/{library}.demux.uniq.filtered.denoised.fasta.gz"
    log:
        "logs/{project}/{library}.denoise.log"
    params:
        ratio = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obiclean"].get("ratio", 0.05),
        distance = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obiclean"].get("distance", 1),
        chimera_detection_flag = lambda wildcards: "--detect-chimera" if config["projects"][wildcards.project]["parameters"]["obiclean"].get("detect_chimera", False) else ""
    threads: lambda wildcards: config["projects"][wildcards.project]["parameters"].get("max-cpu", 1)
    resources:
        runtime = lambda wildcards: config["projects"][wildcards.project]["parameters"]["obiclean"].get("time", 60),
        mem_mb = config["parameters"]["max-cpu"] * config["parameters"]["mem-per-cpu"]
    shell:
        """
        obiclean --max-cpu {threads} \
        --ratio {params.ratio} \
        --head \
        --distance {params.distance} \
        {params.chimera_detection_flag} \
        --compress \
        {input} > {output} 2> {log}
        """