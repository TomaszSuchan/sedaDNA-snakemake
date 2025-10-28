rule dereplicate:
    input:
        "results-sequences/{project}/{library}.demux.fastq.gz"
    output:
        temp("results-sequences/{project}/{library}.demux.uniq.fasta.gz")
    shell:
        """
        obiuniq -m sample {input} | \
        obiannotate -k count -k merged_sample --compress > \
        {output}
        """

rule filter_counts:
    input:
        "results-sequences/{project}/{library}.demux.uniq.fasta.gz"
    output:
        temp("results-sequences/{project}/{library}.demux.uniq.filtered.fasta.gz")
    params:
        min_count = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("min-count", 2),
        min_length = lambda wildcards: config["projects"][wildcards.project]["parameters"]["filtering"].get("min-length", 10)
    shell:
        """
        obigrep --min-length {params.min_length} \
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
    shell:
        """
        obiclean --ratio {params.ratio} \
        --head \
        --distance {params.distance} \
        {params.chimera_detection_flag} \
        --compress \
        {input} > {output} 2> {log}
        """