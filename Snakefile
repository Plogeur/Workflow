workdir: "./"
configfile: "config/config.yaml"
include: "rules/variants_calling.smk"

rule all:
	input:
		config['input_dir']