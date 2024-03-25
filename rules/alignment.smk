# Quality controle
rule fastqc:
	input:
		fastqc=config['sample']
	output:
		html="results/sample_fastqc.html",
		zip="results/sample_fastqc.zip"
	threads: 2
	shell:
		"fastqc {input} --threads {threads} --outdir results"

# Preprocessing reads 
rule reads_trimming:
    input:
        reads="data/sample.fastq"
    output:
        trimmed="results/sample_trimmed.fastq.gz",
        failed="results/sample_failed.fastq.gz",
        report="results/sample_adapters_removal_report.html"
    threads:2
    log:
        "logs/fastp_sample.log"
    shell:
        '''
        FASTP -i {input.reads} -o {output.trimmed} --failed_out {output.failed} \
        --trim_poly_x --overrepresentation_analysis -h {output.report} -w {threads} 2> {log}
        '''

# Indexation of fasta ref
rule bwa_index:
	input:
		ref="data/ref.fa"
	output:
		"results/sample.indexed.bwt"
	params:
		prefix="results/sample.indexed"
	log:
		"logs/bwa_index.log"
	shell:
		"bwa index -p {params.prefix} {input.ref} 2> {log}"

# Map fastq with indexed ref
rule bwa_mapping:
    input:
        rules.bwa_index.output,
        ref=rules.reads_trimming.output
    output:
        "results/sample.sam"
    params:
        prefix="results/sample.indexed"
    threads: 2
    log:
        "logs/bwa_mapping.log"
    shell:
        "bwa mem -t {threads} {params.prefix} {input.ref} > {output} 2> {log}"

# Correcting mate-pair information in paired-end sequencing data
rule fixmate:
	input:
		rules.bwa_mapping.output
	output:
		"results/sample.fixed.sam",
	log:
		"logs/fixmate.log"
	threads: 2
	shell:
		"""
		samtools fixmate --threads {threads} --output-fmt sam {input} {output} 2> {log}
		"""

# Converte SAM to BAM
rule sam2bam:
	input:
		rules.fixmate.output
	output:
		"results/sample.bam",
	log:
		"logs/sam2bam.log"
	shell:
		"GATK --java-options SamFormatConverter --INPUT {input} --OUTPUT {output} 2> {log}"

# Clean SAM file
rule cleansam:
	input:
		genome="data/ref.fa",
		sam=rules.sam2bam.output
	output:
		bam="results/sample.clean.bam"
	log:
		"logs/cleansam.log"
	shell:
		"GATK --java-options CleanSam -R {input.genome} -I {input.sam} -O {output.bam} 2> {log}"