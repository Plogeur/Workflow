#Filter SNPs
rule filter_snps:
    input:
        vcf=rules.index_featurefile.output
    output:
        snps_vcf=temp("results/sample.snps.raw.vcf")
    threads: 4
    params:
        config['filter_snps_options']
    shell:
        """
        $GATK SelectVariants -t {threads} {params} --output {output.snps_vcf} -V {input.vcf}
        """

#Filter SNPs 2 (the return)
rule filter_snps_2:
    input:
        reference="data/ref.fa",
        indels_vcf=rules.filter_snps.output
    output:
        filtered_vcf=temp("results/sample.indels.filtered.vcf")
    shell:
        """
        $GATK VariantFiltration -R {input.reference} -V {input.indels_vcf} \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "Broad_indel_Filter" -O {output.filtered_vcf}
        """

#Filter InDels
rule filter_indel:
    input:
        vcf=rules.index_featurefile.output
    output:
        snps_vcf=temp("results/sample.snps.raw.vcf")
    threads: 4
    shell:
        """
        $GATK SelectVariants -t {threads} --select-type-to-exclude SNP --output {output.snps_vcf} -V {input.vcf}
        """

#Merge SNPs and InDels
rule merge_vcfs:
    input:
        snps_filtered_vcf=rules.filter_snps_2.output,
        indels_filtered_vcf=rules.filter_indel.output
    output:
        variants_filtered_tmp_vcf=temp("results/sample.variants.filtered_tmp.vcf")
    shell:
        """
        $GATK MergeVcfs -I {input.snps_filtered_vcf} -I {input.indels_filtered_vcf} -O {output.variants_filtered_tmp_vcf}
        """

#Keep only filtered variants
rule keep_filtered_variants:
    input:
        reference="data/ref.fa",
        variants_filtered_tmp_vcf=rules.merge_vcfs.output
    output:
        variants_filtered_vcf="results/sample.variants.filtered.vcf"
    shell:
        """
        $GATK SelectVariants -R {input.reference} --variant {input.variants_filtered_tmp_vcf} --exclude-filtered -O {output.variants_filtered_vcf}
        """

#TODO : REPARE THE 2 LAST RULES
rule bcftools_filter:
    input:
        vcf=rules.keep_filtered_variants.output
    output:
        filter=protected("results/reads.trimed.sv.filtered.vcf")
    threads: 2
    params:
        config['bcftools_filter_options']
    shell:
        "bcftools filter --threads {threads} {params} {input} -o {output}"

rule vcat:
    input:
        sorted=rules.samtools_sort.output,
        index=rules.samtools_index.output,
        ref="data/ref.fa"
    output:
        vcf="results/reads.trimed.sv_cutesv.vcf"
    threads: 4
    params:
        config['cutesv_options']
    shell:
        "vcat --threads {threads} {params} --mapped_reads {input.sorted} --vcf {output.vcf}"