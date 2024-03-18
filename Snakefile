workdir: "./"
configfile: "config/config.yaml"
include: "rules/variants_calling.smk"

# usage : 
# snakemake --cores <nb_core_max>

rule all:
	input:
		expand(f"{config['res_dir']}/reads.seqstats", config=config),
        expand(f"{config['res_dir']}/reads.trimed.seqstats", config=config),
        expand(f"{config['res_dir']}/reads.trimed.flagstat", config=config),
        expand(f"{config['res_dir']}/reads.trimed.sv_sniffles.vcf.stats", config=config),
        expand(f"{config['res_dir']}/reads.trimed.snp.vcf.stats", config=config)
