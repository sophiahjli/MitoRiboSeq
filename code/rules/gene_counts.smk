rule gene_counts:
    input:
        genes_bed=genes_bed_file,
        bedgraph=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}_sorted.wig",
        fai=genome_fasta_unzipped + ".fai"
    output:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_{strand}_gene_counts.tsv",
    log:
        log_dir + "/gene_counts/{sample}/{sample}_{mapping_function}_{offset}_map_{strand}_gene_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools map \
            -g {input.fai:q} \
            -a {input.genes_bed:q} \
            -b {input.bedgraph:q} \
            -c 4 \
        > {output:q} 2> {log:q}
        """

rule non_gene_counts:
    input:
        nongenes_bed=nongenes_bed_file,
        bedgraph=config["results_dir"] + "/bedgraph/{sample}_{mapping_function}_{offset}_map_{strand}_sorted.wig",
        fai=genome_fasta_unzipped + ".fai"
    output:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_{strand}_nongene_counts.tsv",
    log:
        log_dir + "/gene_counts/{sample}/{sample}_{mapping_function}_{offset}_map_{strand}_nongene_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools map \
            -g {input.fai:q} \
            -a {input.nongenes_bed:q} \
            -b {input.bedgraph:q} \
            -c 4 \
        > {output:q} 2> {log:q}
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

rule biotype_counts_mitochondrial:
    input:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts_annotated.tsv",
    output:
        config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_mitochondrial.tsv",
    log:
        log_dir + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        echo -e 'gene_biotype\\tcount' > {output:q} && \
        (head -n 1 {input:q} && tail -n +2 {input:q} | sort -k9 | grep "^{config[mito_chrom_id]}" ) \
        | bedtools groupby -g 9 -c 5 -o sum -inheader \
        >> {output:q} 2> {log:q}
        """

rule biotype_counts_cytosolic:
    input:
        config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_combined_gene_counts_annotated.tsv",
    output:
        config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_cytosolic.tsv",
    log:
        log_dir + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        echo -e 'gene_biotype\\tcount' > {output:q} && \
        (head -n 1 {input:q} && tail -n +2 {input:q} | sort -k9 | grep -v "^{config[mito_chrom_id]}" ) \
        | bedtools groupby -g 9 -c 5 -o sum -inheader \
        >> {output:q} 2> {log:q}
        """

rule biotype_counts_with_category:
    input:
        biotype_counts=config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_{location}.tsv",
        biotype_categories=config["biotype_categories"]
    output:
        biotype_counts=config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_with_category_{location}.tsv",
    conda:
        "../envs/csvkit.yml"
    shell:
        """
        csvjoin --tabs --columns 1 {input.biotype_counts:q} {input.biotype_categories:q} \
        | csvformat --out-tabs \
        | sort -k3,3 \
        > {output:q}
        """

rule category_counts:
    input:
        config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_biotype_counts_with_category_{location}.tsv",
    output:
        config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_category_counts_{location}.tsv",
    log:
        log_dir + "/biotype_counts/{sample}_{mapping_function}_{offset}_category_counts_{location}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools groupby -inheader -g 3 -c 2 -o sum -i {input:q} \
        | sed 's/^/{wildcards.location} /' \
        > {output:q} 2> {log:q}
        """

rule summarized_counts:
    input:
        category_counts_mitochondrial=config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_category_counts_mitochondrial.tsv",
        category_counts_cytosolic=config["results_dir"] + "/biotype_counts/{sample}_{mapping_function}_{offset}_category_counts_cytosolic.tsv",
        non_gene_counts_fw=config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_fw_nongene_counts.tsv",
        non_gene_counts_rc=config["results_dir"] + "/gene_counts/{sample}_{mapping_function}_{offset}_map_rc_nongene_counts.tsv",
        samtools_stats=config["results_dir"] + "/samtools_stats/{sample}.txt"
    output:
        config["results_dir"] + "/summarized_counts/{sample}_{mapping_function}_{offset}_summary_counts.tsv",
    conda:
        "../envs/gawk.yml"
    shell:
        """
        cat {input.category_counts_cytosolic:q} {input.category_counts_mitochondrial:q} > {output:q} && \
        cat {input.non_gene_counts_fw:q} {input.non_gene_counts_rc:q} \
            | gawk '{{s+=$4}} END {{printf "Unannotated\\t%.0f\\n", s}}' \
            >> {output:q} && \
        cat {input.samtools_stats:q} \
            | grep -P "^SN\\treads unmapped:" \
            | gawk 'BEGIN{{FS="\\t"}}{{printf "Unmapped\\t%s\\n", $3}}' \
            >> {output:q}
        """
