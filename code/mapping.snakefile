# The main entry point of your workflow for Tn5 RNA-Seq Pipeline
# After configuring, running snakemake -n in a clone of this repository should
# successfully execute a dry-run of the workflow.

include: "common_setup.snakefile"

rule all:
    input:
        #config["results_dir"] + "/multiqc.html",
        #config["results_dir"] + "/combined_gene_counts.tsv"
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()),
        expand(config["results_dir"] + "/mapped/{sample}.bai", sample=samples.keys())

include: "rules/genome.smk"
include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
