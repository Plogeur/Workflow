#Filter SNPs
rule filter_snps:
    input:
        vcf=rules.index_featurefile.output
    output:
        snps_vcf=temp("results/sample.snps.raw.vcf")
    threads: 4
    params:
        config['filter_snps_options']
    log:
        "logs/filter_snps.log"
    shell:
        "GATK SelectVariants -t {threads} {params} --output {output.snps_vcf} -V {input.vcf}"

#Filter SNPs 2
rule filter_snps_2:
    input:
        reference="data/ref.fa",
        indels_vcf=rules.filter_snps.output
    output:
        filtered_vcf=temp("results/sample.indels.filtered.vcf")
    log:
        "logs/filter_snps_2.log"
    shell:
    """
    GATK VariantFiltration -R {input.reference} \   # Reference genome file and input VCF containing indels
        -V {input.indels_vcf} --filter-expression \ # Filtering expression for variant filtration
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum \ # Name for the filter applied to variants that fail the expression
        < -20.0" --filter-name "Broad_indel_Filter" \
        -O {output.filtered_vcf}                    # Output file containing filtered VCF data
    """
"""
QD < 2.0: Filters variants with a low Quality by Depth (QD) value (< 2.0).
FS > 200.0: Filters variants with a high Fisher Strand bias (FS) value (> 200.0).
ReadPosRankSum < -20.0: Filters variants with a low Read Position Rank Sum score (< -20.0).
"""
#Filter InDels
rule filter_indel:
    input:
        vcf=rules.index_featurefile.output
    output:
        snps_vcf=temp("results/sample.snps.raw.vcf")
    threads: 4
    log:
        "logs/filter_indel.log"
    shell:
        "GATK SelectVariants -t {threads} --select-type-to-exclude SNP --output {output.snps_vcf} -V {input.vcf}"

#Merge SNPs and InDels
rule merge_vcfs:
    input:
        snps_filtered_vcf=rules.filter_snps_2.output,
        indels_filtered_vcf=rules.filter_indel.output
    output:
        temp("results/sample.variants.filtered_tmp.vcf")
    log:
        "logs/merge_vcfs.log"
    shell:
        "GATK MergeVcfs -I {input.snps_filtered_vcf} -I {input.indels_filtered_vcf} -O {output}"

#Keep only filtered variants
rule keep_filtered_variants:
    input:
        reference="data/ref.fa",
        variants_filtered_tmp_vcf=rules.merge_vcfs.output
    output:
        "results/sample.variants.filtered.vcf"
    log:
        "logs/keep_filtered_variants.log"
    shell:
        "GATK SelectVariants -R {input.reference} --variant {input.variants_filtered_tmp_vcf} --exclude-filtered -O {output}"