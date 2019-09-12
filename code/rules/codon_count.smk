# Generate all the codon counts to single nucleotide resolution based on the offset
rule codon_counts:
    input: gff=config["genome_annotation_gff3_file"],
           fasta=config["genome_fasta_file"],
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
            --add_three \
            --offset {config[params][plastid][offset]} \
            --{config[params][plastid][mapping_function]} \
            -o "{output}" 2> "{log}"'

rule collapse_codon_counts:
  input: expand(config["results_dir"] + "/codon_count/{sample}_codon_count.txt", sample=samples.keys()),
  output: config["results_dir"] + "/codon_count/All_codoncount_table.txt"
  log: log_dir + "/codon_count/combine_codon_count.log"
  shell:
    """
    awk '
        function basename(file, a, n) {{
            n = split(file, a, "/")
            return a[n]
        }}
    BEGIN {{ OFS="\\t" }}
    {{if (NR > 1) {{ fn=basename(FILENAME); sub("_codon_count.txt", "", fn); print $0, fn }} else {{ print $0, "sample" }} }}' {input:q} >> {output:q} 2> {log:q}
    """
