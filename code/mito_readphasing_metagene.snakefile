include: "common_setup.snakefile"

rule all:
    input:
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()),

        expand(config["results_dir"] + "/metagene/{sample}/{sample}_ND6stop_metagene_overview.png", sample=samples.keys()),
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_ND6stop_metagene_profile.txt", sample=samples.keys()),
        
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_ND4start_metagene_overview.png", sample=samples.keys()),
        expand(config["results_dir"] + "/metagene/{sample}/{sample}_ND4start_metagene_profile.txt", sample=samples.keys()),
        
        expand(config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.png", sample=samples.keys()),
        expand(config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.txt", sample=samples.keys())

include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/metagene.smk"
include: "rules/phasing.smk"