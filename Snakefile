workdir: "./"
configfile: "config/config.yaml"

include: "rules/alignment.smk",
         "rules/variants_calling.smk",
         "rules/variants_analysis.smk"

# usage : 
# snakemake --cores <nb_core_max>

# Rule to run all other rules
"""
rule all:
    input:
        f"{config['res_dir']}/sample.trimed.sam"
        #f"{config['res_dir']}/sample.trimed.aligned.sorted.bam.bai"
rule all:
    input:
        f"{config['res_dir']}/sample.trimed.sam"
"""

rule all:
    input: 
        expand("{sample}.{genome}.bam", genome=GENOMES, sample=SAMPLES)

