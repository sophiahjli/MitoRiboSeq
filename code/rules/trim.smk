rule all_trim:
    input:
        fastq=expand(config["working_dir"] + "/trimmed/{sample}.fastq.gz", sample=samples.keys()),
        qc=expand(config["working_dir"] + "/trimmed/{sample}.qc.txt", sample=samples.keys())

rule cutadapt:
    input:
        get_fastq
    output:
        fastq=config["working_dir"] + "/trimmed/{sample}.fastq.gz",
        qc=config["working_dir"] + "/trimmed/{sample}.qc.txt",
    params:
        extra=config["params"]["cutadapt"]["params"]
    log:
        log_dir + "/cutadapt/{sample}.log"
    threads: 5
    wrapper:
        "0.63.0/bio/cutadapt/se"

