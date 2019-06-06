rule cutadapt:
    input:
        get_fastq
    output:
        fastq=config["working_dir"] + "/trimmed/{sample}.fastq.gz"
    params:
        extra=config["params"]["cutadapt"]["extra"]
    log:
        log_dir + "/cutadapt/{sample}.log"
    threads:
        32
    conda:
        "../envs/cutadapt.yml"
    shell:
      "cutadapt "
      "--adapter='CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG' "
      "--output={output:q} "
      "{params.extra} "
      "{input:q} "
      ">{log} 2>&1"
