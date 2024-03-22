workdir: "./"
configfile: "config/config.yaml"

include: "rules/alignment.smk",
         "rules/variants_calling.smk",
         "rules/variants_analysis.smk"

# usage : 
# snakemake --cores <nb_core_max>

rule all:
    input:
        #f"{config['res_dir']}/sample_fastqc.html",
        f"{config['res_dir']}/reads.trimed.aligned.sorted.bam.bai"
