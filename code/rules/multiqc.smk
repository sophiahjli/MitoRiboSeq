rule multiqc:
    input:
        expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip", sample=samples.keys()),
        expand(config["working_dir"] + "/trimmed/{sample}.qc.txt", sample=samples.keys()),
        expand(config["results_dir"] + "/samtools_stats/{sample}.txt", sample=samples.keys()),
        # expand(config["results_dir"] + "/featurecounts/{sample}_counts.tsv.summary", sample=samples.keys()),
        expand(config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_mqc.tsv",
               sample=samples.keys(),
               offset=config["params"]["plastid"]["offset"],
               mapping_function=config["params"]["plastid"]["mapping_function"]),
    output:
        config["results_dir"] + "/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        log_dir + "/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"

