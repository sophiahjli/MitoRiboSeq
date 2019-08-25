rule metagene_stop_ND6:
    input:
        counts = config["results_dir"] + "/mapped/{sample}.bam",
        counts_bai = config["results_dir"] + "/mapped/{sample}.bai",
        roi = config["metagene_roi_nd6_file"]
    output: 
        profile = config["results_dir"] + "/metagene/{sample}/{sample}_ND6stop_metagene_profile.txt",
        overview = config["results_dir"] + "/metagene/{sample}/{sample}_ND6stop_metagene_overview.png"
    log:
        log_dir + "/metagene/{sample}/metagene.log"
    conda:
        "../envs/plastid.yml"
    shell:
        'metagene count {input.roi:q} '
        '"{config[results_dir]}/metagene/{wildcards.sample}/{wildcards.sample}_ND6stop" '
        '--count_files "{input.counts}" --{config[mapping_function]} '
        '--norm_region 30 200 --min_counts 0 --cmap Blues'

rule metagene_start_ND4:
    input:
        counts = config["results_dir"] + "/mapped/{sample}.bam",
        counts_bai = config["results_dir"] + "/mapped/{sample}.bai",
        roi = config["metagene_roi_nd4_file"]
    output: 
        profile = config["results_dir"] + "/metagene/{sample}/{sample}_ND4start_metagene_profile.txt",
        overview = config["results_dir"] + "/metagene/{sample}/{sample}_ND4start_metagene_overview.png"
    log:
        log_dir + "/metagene/{sample}/metagene.log"
    conda:
        "../envs/plastid.yml"
    shell:
        'metagene count {input.roi:q} \
        "{config[results_dir]}/metagene/{wildcards.sample}/{wildcards.sample}_ND4start" \
        --count_files "{input.counts}" --{config[mapping_function]} \
        --norm_region 30 200 --min_counts 0 --cmap Blues'
