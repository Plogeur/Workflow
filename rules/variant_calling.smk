rule sniffles:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/reads.trimed.sv_sniffles.vcf"
	threads: 4
	params:
		config['sniffles_options']
	shell:
		"""
		sniffles --threads {threads} {params} --reference {input.ref} --sample-id {wildcards.sample}.trimed{wildcards.trim} --input {input.sorted} --vcf {output.vcf} > {log.sniffles} 2>> {log.sniffles}
		"""

rule cutesv:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/reads.trimed.sv_cutesv.vcf"
	threads: 4
	params:
		config['cutesv_options']
	shell:
		"""
		mkdir -p {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}
		cuteSV --threads {threads} {params} --sample {wildcards.sample}.trimed{wildcards.trim} {input.sorted} {input.ref} {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim} 2>> {log.cutesv}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}
		"""

rule svim:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/reads.trimed.sv_svim.vcf"
	threads: 4
	params:
		config['svim_options']
	shell:
		"""
		svim alignment {params} --sample {wildcards.sample}.trimed{wildcards.trim} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim {input.sorted} {input.ref} 2> {log.svim}
		bcftools filter --threads 2 --exclude "SUPPORT<10" --output {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim/variants.vcf 2>> {log.svim}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_svim
		"""

rule nanovar:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		ref=f"{config['ref_dir']}/{config['ref_name']}"
	output:
		vcf=f"{{path_res}}/reads.trimed.sv_nanovar.vcf"
	threads: 4
	params:
		config['nanovar_options']
	shell:
		"""
		nanovar --threads {threads} {params} {input.sorted} {input.ref} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar 2> {log.nanovar}
		echo {wildcards.sample}.trimed{wildcards.trim} > {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt
		bcftools reheader --samples {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt --output {output.vcf} {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar/{wildcards.sample}.trimed{wildcards.trim}.aligned.sorted.nanovar.pass.vcf 2>> {log.nanovar}
		rm -rf {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_workdir_nanovar {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanovar.txt
		"""

rule nanosv:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output,
		bed=f"{{path_res}}/reads.trimed.bedgraph"
	output:
		vcf=f"{{path_res}}/reads.trimed.sv_nanosv.vcf"
	threads: 8
	params:
		config['nanosv_options']
	shell:
		"""
		NanoSV --threads {threads} {params} --bed {input.bed} --output {output.vcf} {input.sorted} 2> {log.nanosv}
		bcftools sort --output {output.vcf} {output.vcf} 2>> {log.nanosv} 
		bcftools filter --threads 2 --exclude "SVLEN<10" --output {output.vcf}.tmp {output.vcf} 2>> {log.nanosv}
		bcftools filter --threads 2 --exclude "DV<10" --output {output.vcf} {output.vcf}.tmp 2>> {log.nanosv}
		echo {wildcards.sample}.trimed{wildcards.trim} > {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt
		bcftools reheader --threads 2 --samples {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt --output {output.vcf} {output.vcf}
		rm {output.vcf}.tmp {wildcards.path_res}/{wildcards.sample}.trimed{wildcards.trim}_name_nanosv.txt
		"""

rule medaka_consensus:
	input:
		sorted=rules.samtools_sort.output,
		index=rules.samtools_index.output
	output:
		hdf=f"{{path_res}}/reads.trimed.hdf"
	threads: 2
	params:
		config['medaka_consensus_options']
	shell:
		"medaka consensus --threads {threads} {params} {input.sorted} {output} 2> {log}"

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

rule medaka_snp:
	input:
		ref=f"{config['ref_dir']}/{config['ref_name']}",
		hdf=rules.medaka_consensus.output
	output:
		vcf=f"{{path_res}}/reads.trimed.snp.vcf"
	threads: 4
	params:
		config['medaka_snp_options']
	shell:
		"medaka snp {params} {input} {output} 2> {log}"