rule seqkit_seq:
	input:
		fastq=config['input_dir']
	output:
		trimed=temp(f"{config['res_dir']}/reads.trimed.fastq.gz")
	threads: 2
	params:
		config['seqkit_seq_options']
	shell:
		"seqkit seq -j {threads} {params} -m {wildcards.trim} {input} -o {output} 2> {log}"

rule minimap2:
	input:
		reads=rules.seqkit_seq.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		sam=temp(f"{{path_res}}/reads.trimed.sam")
	threads: 4
	params:
		config['minimap2_options']
	shell:
		"minimap2 -t {threads} {params} {input.ref} {input.reads} -o {output} 2> {log}"

rule samtools_view:
	input:
		sam=rules.minimap2.output
	output:
		aligned=temp(f"{{path_res}}/reads.trimed.aligned.bam")
	threads: 2
	params:
		config['samtools_view_options']
	shell:
		"samtools view -@ {threads} {params} {input} -o {output} 2> {log}"

rule samtools_sort:
	input:
		aligned=rules.samtools_view.output
	output:
		sorted=f"{{path_res}}/reads.trimed.aligned.sorted.bam"
	threads: 2
	params:
		config['samtools_sort_options']
	shell:
		"samtools sort -@ {threads} {params} {input} -o {output} 2> {log}"

rule samtools_index:
	input:
		sorted=rules.samtools_sort.output
	output:
		index=f"{{path_res}}/reads.trimed.aligned.sorted.bam.bai"
	threads: 2
	params:
		config['samtools_index_options']
	shell:
		"samtools index -@ {threads} {params} {input} 2> {log}"

include: "variants_calling.smk"

rule medaka_stitch:
	input:
		hdf=rules.medaka_consensus.output,
		draft=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		fasta=f"{{path_res}}/reads.trimed.fasta"
	threads: 2
	params:
		config['medaka_stitch_options'],
	shell:
		"medaka stitch --threads {threads} {params} {input.hdf} {input.draft} {output} 2> {log}"