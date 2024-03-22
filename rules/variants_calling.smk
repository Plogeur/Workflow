rule fgbio:
    input:
        sorted=rules.samtools_index.output,
    output:
        vcf=f"{config['res_dir']}/sample.trimed.fgbio.vcf"
    threads: 4
    params:
        config['fgbio_options']
    shell:
        "fgbio --threads {threads} {params} --mapped_reads {input.sorted} --vcf {output.vcf}"

rule gatk_calling:
    input:
        reference=f"{config['ref_dir']}/{config['ref_name']}",
        bam_md_clipped=rules.fgbio.output
    output:
        gvcf=f"{config['res_dir']}/sample.complete.raw.g.vcf"
    threads: 4
    params:
        confidence_threshold="30.0",
        ploidy="{PLOIDY}"
    shell:
        """
        $GATK HaplotypeCaller -t {threads} -R {input.reference} -I {input.bam_md_clipped} -ERC GVCF --output {output.gvcf} --standard-min-confidence-threshold-for-calling {params.confidence_threshold} --dont-use-soft-clipped-bases true --sample-ploidy {params.ploidy}
        """

rule genotype_gvcf:
    input:
        reference=f"{config['ref_dir']}/{config['ref_name']}",
        gvcf=rules.gatk_calling.output
    output:
        vcf=f"{config['res_dir']}/sample.complete.raw.vcf"
    threads: 4
    shell:
        """
        $GATK GenotypeGVCFs -t {threads} -R {input.reference} -V {input.gvcf} -G StandardAnnotation -O {output.vcf}
        """

rule index_featurefile:
    input:
        vcf=rules.genotype_gvcf.output
    run:
        shell("$GATK IndexFeatureFile -I {input.vcf}")