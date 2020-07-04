# vi:syntax=snakemake
# Snakemake script to generate ribosomal occupancy counts using plastid

include: "common_setup.snakefile"

rule all:
   input:
        config["results_dir"] + "/qc/multiqc.html",
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

rule all_bedgraph:
    input:
        expand(config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}.wig",
               sample=samples.keys(),
               offset=config["params"]["plastid"]["offset"],
               mapping_function=config["params"]["plastid"]["mapping_function"],
               strand=['fw', 'rc']),

rule all_biotype_counts:
    input:
        expand(config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts.tsv",
               sample=samples.keys(),
               offset=config["params"]["plastid"]["offset"],
               mapping_function=config["params"]["plastid"]["mapping_function"]),
 

include: "rules/genome.smk"
include: "rules/wiggle.smk"
include: "rules/gene_counts.smk"
include: "rules/codon_count.smk"
include: "rules/codon_occupancy.smk"
include: "rules/multiqc.smk"
include: "rules/biotype_fastqc.smk"
