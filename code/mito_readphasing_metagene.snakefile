include: "common_setup.snakefile"

rule all:
    input:
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()),
        
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_stop_asite_metagene_overview.png", sample=samples.keys()),
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_stop_asite_metagene_profile.txt", sample=samples.keys()),
        
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_start_psite_metagene_overview.png", sample=samples.keys()),
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_start_psite_metagene_profile.txt", sample=samples.keys()),
        
        expand(config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.png", sample=samples.keys()),
        expand(config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.txt", sample=samples.keys())

include: "rules/genome.smk"
include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/metagene.smk"
include: "rules/phasing.smk"
