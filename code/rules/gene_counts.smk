rule gene_counts:
    input:
        genes_bed=genes_bed_file,
        bedgraph=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}.wig",
    output:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_{strand}_gene_counts.tsv",
    log:
        log_dir + "/gene_counts/{sample}/{sample}_{mapping_function}_{offset}_map_{strand}_gene_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools map -a {input.genes_bed:q} -b {input.bedgraph:q} -c 4 > {output:q} 2> {log:q}
        """
       
rule gene_counts_combined:
    input:
        fw=config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_fw_gene_counts.tsv",
        rc=config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_rc_gene_counts.tsv",
    output:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts.tsv",
    log:
        log_dir + "/gene_counts/{sample}/{sample}_{mapping_function}_{offset}_map_combined_gene_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        echo -e 'Chromosome\tStart\tEnd\tgene_id\tcount' > {output:q} &&
        cat {input.fw:q} {input.rc:q} | \
        bedtools sort | \
        sed 's/\.$/0/' | \
        bedtools groupby -g 1-4 -c 7 -o sum \
        >> {output:q} 2> {log:q} 
        """
 
rule gene_counts_annotated:
    input:
        gene_annotation_table=gene_annotation_table,
        gene_counts=config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts.tsv",
    output:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts_annotated.tsv",
    log:
        log_dir + "/gene_counts/{sample}/{sample}_{mapping_function}_{offset}_map_combined_gene_counts_annotated.log",
    conda:
        "../envs/csvkit.yml"
    shell:
        """
        csvjoin \
            --tabs \
            --no-inference \
            --columns gene_id \
            --outer \
            {input.gene_counts:q} \
            {input.gene_annotation_table} | \
        csvformat --out-tabs > \
        {output:q} 2> {log:q}
        """

rule biotype_counts:
    input:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts_annotated.tsv",
    output:
        config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_mqc.tsv",
    log:
        log_dir + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        (head -n 1 {input:q} && tail -n +2 {input:q} | sort -k9 ) | \
        bedtools groupby -g 9 -c 5 -o sum -inheader \
        > {output:q} 2> {log:q}
        """
