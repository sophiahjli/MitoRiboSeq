rule cutadapt:
    input:
        get_fastq
    output:
        fastq=config["working_dir"] + "/trimmed/{sample}.fastq.gz",
        qc=config["working_dir"] + "/trimmed/{sample}.qc.txt",
    params:
        extra=config["params"]["cutadapt"]["extra"]
    log:
        log_dir + "/cutadapt/{sample}.log"
    threads:
        6
    conda:
        "../envs/cutadapt.yml"
    shell:
      "cutadapt "
      "--adapter='CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG' "
      "--output={output.fastq:q} "
      "--cores={threads} "
      "{params.extra} "
      "{input:q} "
      "> {output.qc:q} 2> {log:q}"
