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
        "0.51.3/bio/bwa/index"

rule bwa_aln:
    input:
        [config["working_dir"] + "/trimmed/{sample}.fastq.gz"],
        genome_amb=genome_dir + "/" + genome + ".amb",
        genome_ann=genome_dir + "/" + genome + ".ann",
        genome_bwt=genome_dir + "/" + genome + ".bwt",
        genome_pac=genome_dir + "/" + genome + ".pac",
        genome_sa=genome_dir + "/" + genome + ".sa"
    output:
        [config["working_dir"] + "/trimmed/{sample}.sai"]
    params:
        index=genome_dir + "/" + genome,
        extra=""
    log:
        log_dir + "bwa_aln/{sample}.log"
    threads: 8
    wrapper:
        "0.51.3/bio/bwa/aln"

rule bwa_samse:
    input:
        fastq=[config["working_dir"] + "/trimmed/{sample}.fastq.gz"],
        sai=[config["working_dir"] + "/trimmed/{sample}.sai"],
        genome_amb=genome_dir + "/" + genome + ".amb",
        genome_ann=genome_dir + "/" + genome + ".ann",
        genome_bwt=genome_dir + "/" + genome + ".bwt",
        genome_pac=genome_dir + "/" + genome + ".pac",
        genome_sa=genome_dir + "/" + genome + ".sa"
    output:
        config["results_dir"] + "/mapped/{sample}.bam"
    log:
        log_dir + "bwa_samse/{sample}.log"
    params:
        index=genome_dir + "/" + genome,
        extra=r"-r '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.51.3/bio/bwa/samse"

rule samtools_index:
    input: "{file}.bam"
    output: "{file}.bai"
    params:
        "" # optional params string
    wrapper:
        "0.51.3/bio/samtools/index"
