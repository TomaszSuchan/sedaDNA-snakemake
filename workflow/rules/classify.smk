rule download_ncbitaxo:
	output:
		"data/ncbitaxo.tgz"
	shell:
		"""
		obitaxonomy --download-ncbi --out data/ncbitaxo.tgz
		"""

rule classify:
	input:
		fasta = "results/{PROJECT}/{library}.demux.uniq.filtered.denoised.fasta.gz",
		taxonomy = "data/ncbitaxo.tgz"
	output:
		"results/{PROJECT}/{library}.demux.uniq.filtered.denoised.classified.fasta"
	params:
		reference_db = config["reference_db"]
	shell:
		"""
		# run the classification steps
		obitag -t {input.taxonomy} \
        -R {params.reference_db}\
        {input.fasta} \
        > {output}
		"""