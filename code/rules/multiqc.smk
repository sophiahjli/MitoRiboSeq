
rule stage_multiqc_summary_counts:
    input:
        config["results_dir"] + \
            "/summarized_counts/{sample}_" + \
            config["params"]["plastid"]["mapping_function"] + \
            "_" + \
            str(config["params"]["plastid"]["offset"]) + \
            "_summary_counts.tsv",
    output:
        config["results_dir"] + "/summarized_counts/{sample}_summary_counts_mqc.tsv",
    shell:
        """
        cp {input:q} {output:q}
        """

def get_multiqc_files(wildcards):
    multiqc_files = (
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()) + 
        expand(config["working_dir"] + "/trimmed/{sample}.qc.txt", sample=samples.keys()) +
        expand(config["results_dir"] + "/samtools_stats/{sample}.txt", sample=samples.keys())
    )
    if config["biotype_summary_counts"]:
        multiqc_files.extend(expand(config["results_dir"] + "/summarized_counts/{sample}_summary_counts_mqc.tsv",
                       sample=samples.keys()))
    return(multiqc_files)

rule multiqc:
    input:
        get_multiqc_files
    output:
        config["results_dir"] + "/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        log_dir + "/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"

