workdir: "./"
configfile: "config/config.yaml"

include: "rules/alignment.smk"
include: "rules/variants_calling.smk"
#include: "rules/variants_analysis.smk"

# usage : 
# snakemake --cores <nb_core_max>

rule all:
    input:
        "results/sample.fixed.sam"