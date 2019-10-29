# vi:syntax=snakemake
# Snakemake script to generate ribosomal occupancy counts using plastid

include: "common_setup.snakefile"

rule all:
   input:
        expand(config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_fw.wig",
               sample=samples.keys(),
               offset=config["params"]["plastid"]["offset"],
               mapping_function=config["params"]["plastid"]["mapping_function"]),
        expand(config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_rc.wig",
               sample=samples.keys(),
               offset=config["params"]["plastid"]["offset"],
               mapping_function=config["params"]["plastid"]["mapping_function"]),
        config["results_dir"] + "/codon_count/All_codoncount_table.txt",
        mito_data_raw_coverage_bygene=config["results_dir"] + "/QC/mito_data_raw_coverage_bygene.csv",
        mito_data_raw_depth_bygene=config["results_dir"] + "/QC/mito_data_raw_depth_bygene.csv",
        figures=config["results_dir"] + "/figures",
        mito_occupancy_table=config["results_dir"] + "/tables/mito_occupancy_table.csv",
        mito_cumsum_table=config["results_dir"] + "/tables/mito_cumsum_table.csv"

rule all_qc:
    input:
        mito_data_raw_coverage_bygene=config["results_dir"] + "/QC/mito_data_raw_coverage_bygene.csv",
        mito_data_raw_depth_bygene=config["results_dir"] + "/QC/mito_data_raw_depth_bygene.csv",

rule all_figures:
    input:
        figures=config["results_dir"] + "/figures",

rule all_tables:
    input:
        mito_occupancy_table=config["results_dir"] + "/tables/mito_occupancy_table.csv",
        mito_cumsum_table=config["results_dir"] + "/tables/mito_cumsum_table.csv"


include: "rules/wiggle.smk"
include: "rules/codon_count.smk"
include: "rules/codon_occupancy.smk"
