"""
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

rule reads_trimming:
    input:
        reads="{reads}"
    output:
        trimmed="{SAMPLE_NAME}_trimmed1.fastq.gz",
        failed="{SAMPLE_NAME}_failed.fastq.gz",
        report="{SAMPLE_NAME}_adapters_removal_report.html"
    params:
        threads="{THREADS}"
    log:
        "logs/fastp_{SAMPLE_NAME}.log"
    shell:
        '''
        $FASTP -i {input.reads} -o {output.trimmed} --failed_out {output.failed} \
        --trim_poly_x --overrepresentation_analysis -h {output.report} -w {params.threads}
        '''
"""

#bwa index -p indexed data/ref.fa
#bwa mem -t 2 indexed data/sample.fastq > results/sample.trimed.sam

rule bwa_index:
	input:
		ref=f"{config['ref']}"
	output:
		index=touch("results/indexed")
	log:
		"logs/bwa_index.log"
	shell:
		"bwa index -p {output.index} {input.ref}"

rule bwa_mapping:
	input:
		indexed="results/indexed",
		fastq="data/sample.fastq"
	output:
		sam=f"{config['res_dir']}/sample.trimed.sam"
	threads: 2
	log:
		"logs/bwa_mapping.log"
	shell:
		"bwa mem -t {params.threads} {input.indexed} {input.fastq} > {output.sam}"

# Convert SAM to BAM with samtools view
rule samtools_view:
	input:
		sam=rules.bwa_mapping.output
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