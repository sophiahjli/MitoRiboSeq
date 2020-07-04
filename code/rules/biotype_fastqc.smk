def get_biotypes_for_category(wildcards):
    biotypes=biotype_categories.loc[biotype_categories["category"]==wildcards.category]
    return expand(genome_dir + "/biotype_lists/{location}_{biotype}_list.txt",
                  biotype=biotypes.gene_biotype, **wildcards)


rule all_biotype_bed:
    input:
        expand(genome_dir + "/biotype_lists/{location}_{biotype.gene_biotype}_list.txt",
               biotype=biotype_categories.itertuples(),
               location=['cytosolic', 'mitochondrial']),
        expand(genome_dir + "/biotype_category_lists/{location}_{biotype.category}_list.bed",
               biotype=biotype_categories.itertuples(),
               location=['cytosolic', 'mitochondrial']),
        expand(config["working_dir"] + "/biotype_category_bam/{sample}_{location}_{biotype.category}.bam",
               sample=samples.keys(),
               biotype=biotype_categories.itertuples(),
               location=['cytosolic', 'mitochondrial']),
        expand(config["working_dir"] + "/biotype_category_bam/{sample}_unannotated.bam",
               sample=samples.keys()),
        expand(config["working_dir"] + "/biotype_category_bam/{sample}_unaligned.bam",
               sample=samples.keys()),
        expand(config["results_dir"] + "/qc/biotype_category_fastqc/{sample}_{location}_{biotype.category}_fastqc.zip",
               sample=samples.keys(),
               biotype=biotype_categories.itertuples(),
               location=['cytosolic', 'mitochondrial']),
        expand(config["results_dir"] + "/qc/biotype_category_fastqc/{sample}_unannotated_fastqc.zip",
               sample=samples.keys()),
        expand(config["results_dir"] + "/qc/biotype_category_fastqc/{sample}_unaligned_fastqc.zip",
               sample=samples.keys()),


rule cytosolic_biotype_list:
    input:
        gene_annotation_table
    output:
        genome_dir + "/biotype_lists/cytosolic_{biotype}_list.txt",
    shell:
        """
        grep {wildcards.biotype:q} {input:q} | \
        awk '$3 !~ /^MT.*/' > {output:q}
        """

rule mitochondrial_biotype_list:
    input:
        gene_annotation_table
    output:
        genome_dir + "/biotype_lists/mitochondrial_{biotype}_list.txt"
    shell:
        """
        grep {wildcards.biotype:q} {input:q} | \
        awk '$3 ~ /^MT.*/' > {output:q}
        """

rule biotype_category_bed:
    input:
        biotype_lists=get_biotypes_for_category,
        genome=genome_fasta_unzipped + ".fai"
    output:
        genome_dir + "/biotype_category_lists/{location}_{category}_list.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input.biotype_lists:q} | \
        awk -v OFS='\\t' '{{split($3, region, ":"); split(region[2], location, "-"); print region[1],location[1]-1,location[2],$1,$2,$5}}' | \
        bedtools sort -g {input.genome:q} > \
        {output:q}
        """

rule biotype_category_bam:
    input:
        bed=genome_dir + "/biotype_category_lists/{location}_{category}_list.bed",
        bam=config["results_dir"] + "/mapped/{sample}.bam",
        genome=genome_fasta_unzipped + ".fai"
    output:
        config["working_dir"] + "/biotype_category_bam/{sample}_{location}_{category}.bam"
    log:
        log_dir + "/biotype_category_bam/{sample}_{location}_{category}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.bam:q} -b {input.bed:q} -wa -sorted -g {input.genome:q} > {output:q} 2> {log:q}
        """

rule unannotated_bam:
    input:
        bed=nongenes_bed_file,
        bam=config["results_dir"] + "/mapped/{sample}.bam",
        genome=genome_fasta_unzipped + ".fai"
    output:
        config["working_dir"] + "/biotype_category_bam/{sample}_unannotated.bam"
    log:
        log_dir + "/biotype_category_bam/{sample}_unannotated.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.bam:q} -b {input.bed:q} -wa -sorted -g {input.genome:q} > {output:q} 2> {log:q}
        """


rule unaligned_bam:
    input:
        bam=config["results_dir"] + "/mapped/{sample}.bam",
    output:
        config["working_dir"] + "/biotype_category_bam/{sample}_unaligned.bam"
    log:
        log_dir + "/biotype_category_bam/{sample}_unaligned.log"
    params:
        "-b -f 4" # optional params string
    wrapper:
        "0.61.0/bio/samtools/view"

rule biotype_category_fastqc:
    input:
        config["working_dir"] + "/biotype_category_bam/{sample}_{category}.bam"
    output:
        html=config["results_dir"] + "/qc/biotype_category_fastqc/{sample}_{category}.html",
        zip=config["results_dir"] + "/qc/biotype_category_fastqc/{sample}_{category}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample}_{category}.log"
    threads: 1
    wrapper:
        "0.61.0/bio/fastqc"
