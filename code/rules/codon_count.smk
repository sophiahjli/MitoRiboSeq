# Generate all the codon counts to single nucleotide resolution based on the offset
rule codon_counts:
    input: gff=mito_gff_file,
           fasta=genome_fasta_unzipped,
           counts=config["results_dir"] + "/mapped/{sample}.bam",
           counts_bai=config["results_dir"] + "/mapped/{sample}.bai"
    output: config["results_dir"] + "/codon_count/{sample}_codon_count.txt"
    log: log_dir + "/codon_count/{sample}_codon_count.log"
    conda:
        "../envs/plastid.yml"
    shell: 'python code/BAM2codoncount.py \
            --bam "{input.counts}" \
            --gff "{input.gff}" \
            --fasta "{input.fasta}" \
            --min_length={config[params][plastid][min_length]} \
            --max_length={config[params][plastid][max_length]} \
            --offset {config[params][plastid][offset]} \
            --{config[params][plastid][mapping_function]} \
            -o "{output}" 2> "{log}"'

rule collapse_codon_counts:
    input: expand(config["results_dir"] + "/codon_count/{sample}_codon_count.txt", sample=samples.keys()),
    output: temp(config["results_dir"] + "/codon_count/All_codoncount_table.tmp")
    log: log_dir + "/codon_count/combine_codon_count.log"
    group: "all_codon_count"
    conda:
        "../envs/gawk.yml"
    shell:
        """
        head -1 {input[0]:q} | sed 's/$/\\tsample/' > {output:q}
        gawk '
            function basename(file, a, n) {{
                n = split(file, a, "/")
                return a[n]
            }}
        BEGIN {{ OFS="\\t" }}
        {{if ($1 !~ "transcript_id") {{ fn=basename(FILENAME); sub("_codon_count.txt", "", fn); print $0, fn }} }}' {input:q} >> {output:q} 2> {log:q}
        """

rule add_gene_names:
    input:
        gene_annotations=gene_annotation_table,
        codon_count_table=config["results_dir"] + "/codon_count/All_codoncount_table.tmp"
    output:
        config["results_dir"] + "/codon_count/All_codoncount_table.txt"
    log:
        log_dir + "/codon_count/add_gene_names.log"
    group: "all_codon_count"
    conda:
        "../envs/csvkit.yml"
    shell:
        """
        csvjoin --tabs --left --columns gene_id {input.codon_count_table:q} {input.gene_annotations:q} | 
        csvcut --not-columns Chromosome |
        csvformat --out-tabs > {output:q} 2> {log:q}
        """
