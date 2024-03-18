rule bcftools_filter:
	input:
		vcf=rules.sniffles.output
	output:
		filter=protected(f"{{path_res}}/reads.trimed.sv.filtered.vcf")
	threads: 2
	params:
		config['bcftools_filter_options']
	shell:
		"bcftools filter --threads {threads} {params} {input} -o {output} 2> {log}"

rule bgzip:
	input:
		vcf=rules.bcftools_filter.output
	output:
		gz=f"{{path_res}}/reads.trimed.sv.filtered.vcf.gz"
	threads: 2
	params:
		config['bgzip_options']
	shell:
		"bgzip -@ {threads} {input} 2> {log}"

rule tabix:
	input:
		vcf_gz=rules.bgzip.output
	output:
		index=f"{{path_res}}/reads.trimed.sv.filtered.vcf.gz.tbi"
	threads: 2
	params:
		config['tabix_options']
	shell:
		"tabix {input} 2> {log}"

rule combisv:
	input:
		sniffles=rules.sniffles.output,
		cutesv=rules.cutesv.output,
		nanosv=rules.nanosv.output,
		svim=rules.svim.output
	output:
		vcf=f"{{path_res}}/reads.trimed.combisv.vcf",
		sniffles=f"{{path_res}}/reads.trimed.combisv_sniffles.vcf",
		cutesv=f"{{path_res}}/reads.trimed.combisv_cutesv.vcf",
		nanosv=f"{{path_res}}/reads.trimed.combisv_nanosv.vcf",
		svim=f"{{path_res}}/reads.trimed.combisv_svim.vcf"
	log:
		f"{{path_res}}/logs/combisv/reads.trimed.log"
	threads: 1
	params:
		config['combisv_options']
	shell:
		"""
		perl combiSV2.2.pl -sniffles {input.sniffles} -cutesv {input.cutesv} -svim {input.svim} -nanosv {input.nanosv} -o {wildcards.sample}.trimed{wildcards.trim}.combisv.vcf 2> {log}
		mv {wildcards.sample}.trimed{wildcards.trim}.combisv.vcf {output.vcf}
		mv Sniffles_{wildcards.sample}.trimed{wildcards.trim}.combisv.vcf {output.sniffles}
		mv cuteSV_{wildcards.sample}.trimed{wildcards.trim}.combisv.vcf {output.cutesv}
		mv SVIM_{wildcards.sample}.trimed{wildcards.trim}.combisv.vcf {output.svim}
		mv NanoSV_{wildcards.sample}.trimed{wildcards.trim}.combisv.vcf {output.nanosv}
		rm simplified_{wildcards.sample}.trimed{wildcards.trim}.combisv.vcf
		"""