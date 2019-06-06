rule bwa_index:
    input:
        "{genome}" + genome_ext
    output:
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt",
        "{genome}.pac",
        "{genome}.sa"
    log:
        "logs/bwa_index/{genome}.log"
    params:
        prefix="{genome}",
        algorithm="bwtsw"
    wrapper:
        "0.35.0/bio/bwa/index"

rule bwa_mem:
    input:
        reads=[config["working_dir"] + "/trimmed/{sample}.fastq.gz"],
        genome_amb=genome_dir + "/" + genome + ".amb",
        genome_ann=genome_dir + "/" + genome + ".ann",
        genome_bwt=genome_dir + "/" + genome + ".bwt",
        genome_pac=genome_dir + "/" + genome + ".pac",
        genome_sa=genome_dir + "/" + genome + ".sa"
    output:
        config["results_dir"] + "/mapped/{sample}.bam"
    log:
        log_dir + "bwa_mem/{sample}.log"
    params:
        index=genome_dir + "/" + genome,
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.35.0/bio/bwa/mem"
