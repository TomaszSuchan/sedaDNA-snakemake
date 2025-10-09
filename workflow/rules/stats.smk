rule raw_stats:
    input:
        lambda wildcards: config["libraries"][wildcards.library]["forward"]
    output:
        "stats/{PROJECT}/{library}.raw_stats.json"
    shell:
        """
        obisummary {input} > {output}
        """

rule pair_stats:
    input:
        "results/{PROJECT}/{library}.paired.fastq.gz"
    output:
        "stats/{PROJECT}/{library}.pair_stats.json"
    shell:
        """
        obisummary {input} > {output}
        """

rule demux_stats:
    input:
        "results/{PROJECT}/{library}.demux.fastq.gz"
    output:
        "stats/{PROJECT}/{library}.demux_stats.json"
    shell:
        """
        obisummary {input} > {output}
        """

rule merge_json_read_counts:
    input:
        raw="stats/{PROJECT}/{library}.raw_stats.json",
        pair="stats/{PROJECT}/{library}.pair_stats.json",
        demux="stats/{PROJECT}/{library}.demux_stats.json"
    output:
        "stats/{PROJECT}/{library}.merged_stats.tsv"
    run:
        import json, yaml, pandas as pd
        
        def load_stats(path):
            with open(path) as f:
                text = f.read()
            try:
                return json.loads(text)
            except json.JSONDecodeError:
                return yaml.safe_load(text)
        
        raw_data = load_stats(input.raw)
        pair_data = load_stats(input.pair)
        demux_data = load_stats(input.demux)
        
        raw_reads = raw_data.get("count", {}).get("reads", "NA")
        pair_reads = pair_data.get("count", {}).get("reads", "NA")
        demux_reads = demux_data.get("count", {}).get("reads", "NA")
        
        df = pd.DataFrame(
            [{
                "library": wildcards.library,
                "raw_reads": raw_reads,
                "pair_reads": pair_reads,
                "demux_reads": demux_reads,
            }]
        )
        
        df.to_csv(output[0], sep="\t", index=False)