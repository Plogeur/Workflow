# Indexation of fasta ref
rule bwa_index:
	input:
		ref="data/ref.fa"
	output:
		"results/indexed.bwt"
	params:
		prefix="results/indexed"
	log:
		"logs/bwa_index.log"
	shell:
		"bwa index -p {params.prefix} {input.ref}"

# Map fastq with indexed ref
rule bwa_mapping:
    input:
        rules.bwa_index.output,
        ref="data/sample.fastq"
    output:
        "results/indexed.trimed.sam"
    params:
        prefix="results/indexed"
    threads: 2
    log:
        "logs/bwa_mapping.log"
    shell:
        "bwa mem -t {threads} {params.prefix} {input.ref} > {output}"

# Convert SAM to BAM with samtools view
rule samtools_view:
	input:
		rules.bwa_mapping.output
	output:
		"results/sample.trimed.aligned.bam"
	threads: 2
	params:
		config['samtools_view_options']
	shell:
		"samtools view -@ {threads} {params} {input} -o {output}"

# Sort aligned BAM with samtools sort
rule samtools_sort:
	input:
		rules.samtools_view.output
	output:
		"results/sample.trimed.aligned.sorted.bam"
	threads: 2
	params:
		config['samtools_sort_options']
	shell:
		"samtools sort -@ {threads} {params} {input} -o {output}"

# Index sorted BAM with samtools index
rule samtools_index:
	input:
		rules.samtools_sort.output
	output:
		"results/sample.trimed.aligned.sorted.bam.bai"
	threads: 2
	params:
		config['samtools_index_options']
	shell:
		"samtools index -@ {threads} {params} {input}"
