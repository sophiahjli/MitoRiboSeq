

configfile:"mito_config.yml"
SAMPLES = config["samples"]
output_dir = config["sample_dir"]
GFF_FILE="{genome_dir}/{annot_file}.gff3".format(genome_dir=config["genome_dir"], annot_file=config["annot_file"])


rule all:
	input:
		expand("{output_dir}/metagene/{sample}/{sample}_metagene_overview.png",output_dir = output_dir, sample=SAMPLES),
        expand("{output_dir}/metagene/{sample}/{sample}_metagene_profile.txt",output_dir = output_dir, sample=SAMPLES),

        expand("{output_dir}/phasing_analysis/{sample}/{sample}_phasing.png",output_dir = output_dir, sample=SAMPLES),
        expand("{output_dir}/phasing_analysis/{sample}/{sample}_phasing.txt",output_dir = output_dir, sample=SAMPLES),

rule index_bam:
    input: "{file}.bam"
    output: "{file}.bai"
    shell: 'samtools index "{input}" "{output}"'

rule metagene_stop_ND6:
	input:
		counts = "{output_dir}/mapping/{sample}.bam",
		counts_bai = "{output_dir}/mapping/{sample}.bai"

	output: 
		profile = "{output_dir}/metagene/{sample}/{sample}_metagene_profile.txt",
		overview = "{output_dir}/metagene/{sample}/{sample}_metagene_overview.png"


	log:
		"{output_dir}/metagene/{sample}/metagene.log"


	shell:
		'metagene count "genome/chrM_orfs_stop_rois_ND6.txt" \
		"{wildcards.output_dir}/metagene/{wildcards.sample}/{wildcards.sample}" \
		--count_files "{input.counts}" --{config[mapping_function]} \
		--norm_region 30 200 --min_counts 10 --cmap Blues'


rule phasing_analysis:
    input: 
            counts = "{output_dir}/mapping/{sample}.bam",
            gff=GFF_FILE
    output: 
            profile = "{output_dir}/phasing_analysis/{sample}/{sample}_phasing.txt",
            overview = "{output_dir}/phasing_analysis/{sample}/{sample}_phasing.png"

    shell: 'phase_by_size {wildcards.output_dir}/phasing_analysis/{wildcards.sample}/{wildcards.sample} \
     --count_files {input.counts}  \
     --annotation_files {input.gff} \
     --annotation_format GFF3 \
     --{config[mapping_function]}  \
     --min_length {config[min_length]} \
     --max_length {config[max_length]} '