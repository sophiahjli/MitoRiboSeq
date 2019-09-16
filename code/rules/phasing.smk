rule phasing_analysis:
    input: 
        counts=config["results_dir"] + "/mapped/{sample}.bam",
        counts_index=config["results_dir"] + "/mapped/{sample}.bai",
        gff=config["genome_annotation_gff3_file"]
    output: 
        profile = config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.txt",
        overview = config["results_dir"] + "/phasing_analysis/{sample}/{sample}_phasing.png"
    conda:
        "../envs/plastid.yml"
    shell: 'phase_by_size {config[results_dir]}/phasing_analysis/{wildcards.sample}/{wildcards.sample} \
     --count_files {input.counts}  \
     --annotation_files {input.gff} \
     --annotation_format GFF3 \
     --{config[params][plastid][mapping_function]}  \
     --min_length {config[params][plastid][min_length]} \
     --max_length {config[params][plastid][max_length]} '
