rule metagene_stop:
    input:
        counts = config["results_dir"] + "/mapped/{sample}.bam",
        counts_bai = config["results_dir"] + "/mapped/{sample}.bai",
        roi = stop_site_base + "_rois.txt"
    output: 
        profile = config["results_dir"] + "/metagene/{sample}/{sample}_stop_asite_metagene_profile.txt",
        overview = config["results_dir"] + "/metagene/{sample}/{sample}_stop_asite_metagene_overview.png"
    log:
        log_dir + "/metagene/{sample}/metagene_stop_asite.log"
    conda:
        "../envs/plastid.yml"
    shell:
        'metagene count {input.roi:q} '
        '"{config[results_dir]}/metagene/{wildcards.sample}/{wildcards.sample}_stop_asite" '
        '--count_files "{input.counts}" --{config[params][plastid][mapping_function]} '
        '--norm_region 30 200 --min_counts 0 --cmap Blues '
        '2> {log:q}'

rule metagene_start:
    input:
        counts = config["results_dir"] + "/mapped/{sample}.bam",
        counts_bai = config["results_dir"] + "/mapped/{sample}.bai",
        roi = start_site_base + "_rois.txt"
    output: 
        profile = config["results_dir"] + "/metagene/{sample}/{sample}_start_psite_metagene_profile.txt",
        overview = config["results_dir"] + "/metagene/{sample}/{sample}_start_psite_metagene_overview.png"
    log:
        log_dir + "/metagene/{sample}/metagene_start_psite.log"
    conda:
        "../envs/plastid.yml"
    shell:
        'metagene count {input.roi:q} '
        '"{config[results_dir]}/metagene/{wildcards.sample}/{wildcards.sample}_start_psite" '
        '--count_files "{input.counts}" --{config[params][plastid][mapping_function]} '
        '--norm_region 30 200 --min_counts 0 --cmap Blues '
        '2> {log:q}'
