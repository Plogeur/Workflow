#ClipOverlaps
rule fgbio:
    input:
        reference="data/ref.fa",
        sorted=rules.samtools_index.output,
    output:
        temp("results/sample.trimed.fgbio.vcf")
    threads: 4
    params:
        config['fgbio_options']
    shell:
        "fgbio ClipBam -i {input.sorted} -r {input.reference} -o {output} {params} sample_fgbio_metrics.txt"

#Perform variant calling with HaplotypeCaller
rule gatk_calling:
    input:
        reference="data/ref.fa",
        bam_md_clipped=rules.fgbio.output
    output:
        gvcf=temp("results/sample.complete.raw.g.vcf")
    threads: 4
    params:
        config['gatk_calling_options']
    shell:
        """
        $GATK HaplotypeCaller -t {threads} -R {input.reference} -I {input.bam_md_clipped} -ERC GVCF --output {output.gvcf} {params}"
        """

#Genotype gVCF
rule genotype_gvcf:
    input:
        reference="data/ref.fa",
        gvcf=rules.gatk_calling.output
    output:
        vcf=temp("results/sample.complete.raw.vcf")
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
        index_vcf="results/sample.complete.indexed.raw.vcf"
    shell:
        "$GATK IndexFeatureFile -I {input.gvcf} -o {output.index_vcf}"