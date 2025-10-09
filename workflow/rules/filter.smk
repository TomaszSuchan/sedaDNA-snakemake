rule dereplicate:
    input:
        "results/{PROJECT}/{library}.demux.fastq.gz"
    output:
        temp("results/{PROJECT}/{library}.demux.uniq.fasta.gz")
    shell:
        """
        obiuniq -m sample {input} | \
        obiannotate -k count -k merged_sample --compress > \
        {output}
        """

rule filter_counts:
    input:
        "results/{PROJECT}/{library}.demux.uniq.fasta.gz"
    output:
        temp("results/{PROJECT}/{library}.demux.uniq.filtered.fasta.gz")
    params:
        min_count = config["parameters"]["filtering"].get("min-count", 2),
        min_length = config["parameters"]["filtering"].get("min-length", 10)
    shell:
        """
        obigrep --min-length {params.min_length} \
        --min-count {params.min_count} {input} > {output}
        """

rule denoise:
    input:
        "results/{PROJECT}/{library}.demux.uniq.filtered.fasta.gz"
    output:
        "results/{PROJECT}/{library}.demux.uniq.filtered.denoised.fasta.gz"
    log:
        "logs/{PROJECT}/{library}.denoise.log"
    params:
        ratio = config["parameters"]["obiclean"].get("ratio", 0.05),
        distance = config["parameters"]["obiclean"].get("distance", 1),
        chimera_detection_flag=lambda wildcards: "--detect-chimera" if config["parameters"]["obiclean"].get("detect_chimera", False) else "",

    shell:
        """
        obiclean --ratio {params.ratio} \
        --head \
        --distance {params.distance} \
        {params.chimera_detection_flag} \
        --compress \
        {input} > {output} 2> {log}
        """
