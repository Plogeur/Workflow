# Read's quality control
rule fastqc:
	input:
		fastqc=config['sample']
	output:
		html=f"{config['res_dir']}/sample_fastqc.html",
		zip=f"{config['res_dir']}/sample_fastqc.zip"
	threads: 2
	shell:
		"fastqc {input} --threads {threads} --outdir results"

# Sequence trimming
rule seqkit_seq:
	input:
		fastq=config['sample']
	output:
		trimed=temp(f"{config['res_dir']}/sample.trimed.fastq.gz")
	threads: 2
	params:
		config['seqkit_seq_options']
	shell:
		"seqkit seq -j {threads} {params} -m 50 {input} -o {output}"

# Mapping with bwa
rule bwa:
	input:
		reads=rules.seqkit_seq.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		sam=temp(f"{config['res_dir']}/sample.trimed.sam")
	threads: 4
	shell:
		"bwa mem -t {threads} {params} -R {input.ref} {input.reads} -o {output}"

# Convert SAM to BAM with samtools view
rule samtools_view:
	input:
		sam=rules.bwa.output
	output:
		aligned=temp(f"{config['res_dir']}/sample.trimed.aligned.bam")
	threads: 2
	params:
		config['samtools_view_options']
	shell:
		"samtools view -@ {threads} {params} {input} -o {output}"

# Sort aligned BAM with samtools sort
rule samtools_sort:
	input:
		aligned=rules.samtools_view.output
	output:
		sorted=f"{config['res_dir']}/sample.trimed.aligned.sorted.bam"
	threads: 2
	params:
		config['samtools_sort_options']
	shell:
		"samtools sort -@ {threads} {params} {input} -o {output}"

# Index sorted BAM with samtools index
rule samtools_index:
	input:
		sorted=rules.samtools_sort.output
	output:
		index=f"{config['res_dir']}/sample.trimed.aligned.sorted.bam.bai"
	threads: 2
	params:
		config['samtools_index_options']
	shell:
		"samtools index -@ {threads} {params} {input}"