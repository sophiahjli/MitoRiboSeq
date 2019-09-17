rule wiggle_specify_mapping:
    input: counts=config["results_dir"] + "/mapped/{sample}.bam",
           counts_bai=config["results_dir"] + "/mapped/{sample}.bai",
    output: fw=config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_fw.wig",
            rc=config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map_rc.wig",
    conda:
        "../envs/plastid.yml"            
    params: base=config["results_dir"] + "/wiggle/{sample}_{mapping_function}_{offset}_map"
    shell: 'make_wiggle -o "{params.base}" \
            --count_files "{input.counts}" \
            --min_length {config[params][plastid][min_length]} \
            --max_length {config[params][plastid][max_length]} \
            --{wildcards.mapping_function} \
            --offset {wildcards.offset} \
            --nibble {wildcards.offset}'
