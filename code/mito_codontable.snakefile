# vi:syntax=snakemake
# Snakemake script to generate ribosomal occupancy counts using plastid

# Required inputs:
#   genome/chrM.gff3
#   genome/chrM.fa

#   data/mapping/SAMPLE.bam

# Software dependencies
#   Plastid
#   samtools

import glob

configfile:"code/mito_config.yml"
SAMPLES = config["samples"] 
output_dir = config["sample_dir"]


LENGTHS=range(config['min_length'], config['max_length'] + 1)
FASTA_FILE="{genome_dir}/{genome}.fa".format(genome_dir=config["genome_dir"], genome=config["genome"])
GFF_FILE="{genome_dir}/{annot_file}.gff3".format(genome_dir=config["genome_dir"], annot_file=config["annot_file"])


rule all:
   input:
          expand("{output_dir}/wiggle/{sample}_{mapping_function}_{offset}_map_fw.wig",output_dir = output_dir ,sample=SAMPLES, offset = config["offset"], mapping_function = config["mapping_function"]),
          expand("{output_dir}/wiggle/{sample}_{mapping_function}_{offset}_map_rc.wig",output_dir = output_dir, sample=SAMPLES, offset = config["offset"], mapping_function = config["mapping_function"]),

          expand("{output_dir}/codon_count/{sample}_codon_count.txt",output_dir = output_dir, sample=SAMPLES),
          expand("{output_dir}/codon_count/All_codoncount_table.txt",output_dir = output_dir)

rule wiggle_specify_mapping:
    input: counts="{output_dir}/mapping/{sample}.bam",
           counts_bai="{output_dir}/mapping/{sample}.bai",

    output: fw="{output_dir}/wiggle/{sample}_{mapping_function}_{offset}_map_fw.wig",
            rc="{output_dir}/wiggle/{sample}_{mapping_function}_{offset}_map_rc.wig",
            
    params: base="{output_dir}/wiggle/{sample}_{mapping_function}_{offset}_map"
    shell: 'make_wiggle -o "{params.base}" \
            --count_files "{input.counts}" \
            --min_length {config[min_length]} \
            --max_length {config[max_length]} \
            --{wildcards.mapping_function} \
            --offset {wildcards.offset} \
            --nibble {wildcards.offset}'

# Generate all the codon counts to single nucleotide resolution based on the offset
rule codon_counts:
    input: gff=GFF_FILE,
           fasta=FASTA_FILE,
           counts="{output_dir}/mapping/{sample}.bam",
           counts_bai="{output_dir}/mapping/{sample}.bai"
    output: "{output_dir}/codon_count/{sample}_codon_count.txt"
    log: "{output_dir}/codon_count/{sample}_codon_count.log"
    shell: 'python code/BAM2codoncount.py \
            --bam "{input.counts}" \
            --gff "{input.gff}" \
            --fasta "{input.fasta}" \
            --min_length={config[min_length]} \
            --max_length={config[max_length]} \
            --add_three \
            --offset {config[offset]} \
            --{config[mapping_function]} \
            -o "{output}" 2> "{log}"'

# Takes all the individual codon count table to one
######### Not yet done as I haven't modified the codes!######
rule collapse_codon_counts:
  input: counts_table = "{output_dir}/codon_count/{sample}_codon_count.txt"
  output: "{output_dir}/codon_count/All_codoncount_table.txt"
  log: "{output_dir}/codon_count/{sample}_codon_count.log"
  shell: 'python code/CollapseCodonCounts_MultipleFiles.py \
            --sample "{input.counts_table}" \
            -o "{output}" 2> "{log}"'
