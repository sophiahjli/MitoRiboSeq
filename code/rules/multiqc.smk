
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

rule multiqc:
    input:
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()),
        expand(config["working_dir"] + "/trimmed/{sample}.qc.txt", sample=samples.keys()),
        expand(config["results_dir"] + "/samtools_stats/{sample}.txt", sample=samples.keys()),
        # expand(config["results_dir"] + "/featurecounts/{sample}_counts.tsv.summary", sample=samples.keys()),
        expand(config["results_dir"] + "/summarized_counts/{sample}_summary_counts_mqc.tsv",
               sample=samples.keys()),
    output:
        config["results_dir"] + "/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        log_dir + "/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"

