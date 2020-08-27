rule bedgraph_specify_mapping:
    input: counts=config["results_dir"] + "/mapped_unique/{sample}.bam",
           counts_bai=config["results_dir"] + "/mapped_unique/{sample}.bai",
    output: fw=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_fw.wig",
            rc=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_rc.wig",
    conda:
        "../envs/plastid.yml"            
    params: base=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map"
    shell: 
       """
        make_wiggle \
            -o {params.base:q} \
            --count_files {input.counts:q} \
            --min_length {config[params][plastid][min_length]} \
            --max_length {config[params][plastid][max_length]} \
            --{wildcards.mapping_function} \
            --offset {wildcards.offset} \
            --output_format bedgraph
        """

rule sort_bedgraph:
    input:
        bedgraph=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}.wig",
        fai=genome_fasta_unzipped + ".fai"
    output:
        config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}_sorted.wig",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools sort \
            -faidx {input.fai:q} \
            -i {input.bedgraph:q} > \
            {output:q}
        """

