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
        ref="data/sample.fastq"
    output:
        "results/sample.sam"
    params:
        prefix="results/sample.indexed"
    threads: 2
    log:
        "logs/bwa_mapping.log"
    shell:
        "bwa mem -t {threads} {params.prefix} {input.ref} > {output} 2> {log}"

# Rule for samtools fixmate
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

# Rule for GATK SamFormatConverter
rule sam2bam:
	input:
		rules.fixmate.output
	output:
		"results/sample.bam",
	log:
		"logs/sam2bam.log"
	shell:
		"GATK --java-options SamFormatConverter --INPUT {input} --OUTPUT {output} 2> {log}"

# Rule for GATK CleanSam
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