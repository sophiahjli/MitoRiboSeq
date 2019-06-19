# vi:syntax=snakemake
# Snakemake script to generate ribosomal occupancy counts using plastid

include: "common_setup.snakefile"

rule all:
   input:
          expand(config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_fw.wig",sample=samples.keys(), offset = config["offset"], mapping_function = config["mapping_function"]),
          expand(config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_rc.wig", sample=samples.keys(), offset = config["offset"], mapping_function = config["mapping_function"]),
          config["results_dir"] + "/codon_count/All_codoncount_table.txt"

include: "rules/wiggle.smk"
include: "rules/codon_count.smk"
