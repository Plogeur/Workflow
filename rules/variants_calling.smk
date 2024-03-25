"""
#ClipOverlaps
rule fgbio:
    input:
        reference="data/ref.fa",
        sorted=rules.samtools_index.output,
    output:
        temp("results/sample.trimed.fgbio.vcf")
    threads: 2
    params:
        config['fgbio_options']
    shell:
        "fgbio ClipBam -i {input.sorted} -r {input.reference} -o {output} {params} sample_fgbio_metrics.txt"
"""

# Perform variant calling with HaplotypeCaller
# gVCF : variants génétiques sur la couverture à travers tout le génome
rule gatk_calling:
    input:
        reference="data/ref.fa",
        bam_md_clipped=rules.cleansam.output
    output:
        gvcf="results/sample.complete.raw.gvcf"
    threads: 4
    params:
        config['gatk_calling_options']
    log:
        "logs/gatk_calling.log"
    shell:
        """
        GATK HaplotypeCaller -t {threads} -R {input.reference} -I {input.bam_md_clipped} -ERC GVCF --output {output.gvcf} {params} 2> {log}"
        """

# Genotype VCF
# VCF : stocke uniquement des informations sur les variants génétiques détectés.
rule genotype_gvcf:
    input:
        reference="data/ref.fa",
        gvcf=rules.gatk_calling.output
    output:
        vcf=temp("results/sample.complete.raw.vcf")
    threads: 4
    log:
        "logs/genotype_gvcf.log"
    shell:
        "GATK GenotypeGVCFs -t {threads} -R {input.reference} -V {input.gvcf} -G StandardAnnotation -O {output.vcf} 2> {log}"

# Index featurefile
rule index_featurefile:
    input:
        gvcf=rules.genotype_gvcf.output
    output:
        index_vcf="results/sample.complete.indexed.raw.vcf"
    log:
        "logs/index_featurefile.log"
    shell:
        "GATK IndexFeatureFile -I {input.gvcf} -o {output.index_vcf} 2> {log}"