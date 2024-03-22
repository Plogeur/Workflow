#ClipOverlaps
rule fgbio:
    input:
        reference=f"{config['ref_dir']}/{config['ref_name']}",
        sorted=rules.samtools_index.output,
    output:
        vcf=temp(f"{config['res_dir']}/sample.trimed.fgbio.vcf")
    threads: 4
    params:
        config['fgbio_options']
    shell:
        "fgbio ClipBam --threads {threads} -i {input.sorted} -r {reference} -o {output} {params} sample_fgbio_metrics.txt"

#Perform variant calling with HaplotypeCaller
rule gatk_calling:
    input:
        reference=f"{config['ref_dir']}/{config['ref_name']}",
        bam_md_clipped=rules.fgbio.output
    output:
        gvcf=temp(f"{config['res_dir']}/sample.complete.raw.g.vcf")
    threads: 4
    params:
        config['gatk_calling_options']
    shell:
        """
        $GATK HaplotypeCaller -t {threads} -R {input.reference} -I {input.bam_md_clipped} -ERC GVCF --output {output.gvcf} {gatk_calling_options}"
        """

#Genotype gVCF
rule genotype_gvcf:
    input:
        reference=f"{config['ref_dir']}/{config['ref_name']}",
        gvcf=rules.gatk_calling.output
    output:
        vcf=temp(f"{config['res_dir']}/sample.complete.raw.vcf")
    threads: 4
    shell:
        """
        $GATK GenotypeGVCFs -t {threads} -R {input.reference} -V {input.gvcf} -G StandardAnnotation -O {output.vcf}
        """
        
#Index featurefile
rule index_featurefile:
    input:
        gvcf=rules.genotype_gvcf.output
    output:
        index_vcf=f"{config['res_dir']}/sample.complete.indexed.raw.vcf"
    shell:
        "$GATK IndexFeatureFile -I {input.gvcf} -o {output.index_vcf}"