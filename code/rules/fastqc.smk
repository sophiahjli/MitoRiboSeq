rule all_fastqc:
    input:
        html=expand(config["working_dir"] + "/fastqc/{sample}_fastqc.html", sample=samples.keys()),
        zip=expand(config["working_dir"] + "/fastqc/{sample}_fastqc.zip" , sample=samples.keys())


rule fastqc:
    input:
        get_fastq
    output:
        html=config["working_dir"] + "/fastqc/{sample}_fastqc.html",
        zip=config["working_dir"] + "/fastqc/{sample}_fastqc.zip"
    log:
        log_dir + "/fastqc/{sample}.log"
    params:
        config["params"]["fastqc"]["extra"]
    threads: 1
    wrapper:
        "0.63.0/bio/fastqc"

