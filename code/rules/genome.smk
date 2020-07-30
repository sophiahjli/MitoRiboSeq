from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

ruleorder: get_fasta_ftp > gunzip

rule all_downloads:
    input:
        config["genome_fasta_file"],
        config["genome_annotation_gff3_file"],
        config["genome_annotation_gtf_file"],


rule get_fasta_ftp:
    input:
        FTP.remote(config["genome_fasta_ftp"], keep_local=True)
    output:
        config["genome_fasta_file"]
    shell:
        "mv {input:q} {output:q}"

rule index_fasta:
    input:
        genome_fasta_unzipped
    output:
        genome_fasta_unzipped + ".fai"
    params:
        ""
    wrapper:
        "0.60.1/bio/samtools/faidx"


rule get_gff_ftp:
    input:
        FTP.remote(config["genome_annotation_gff3_ftp"], keep_local=True)
    output:
        config["genome_annotation_gff3_file"]
    shell:
        "mv {input:q} {output:q}"

rule get_gtf_ftp:
    input:
        FTP.remote(config["genome_annotation_gtf_ftp"], keep_local=True)
    output:
        config["genome_annotation_gtf_file"]
    shell:
        "mv {input:q} {output:q}"

rule mito_gff_file:
    input:
        config["genome_annotation_gff3_file"]
    output:
        mito_gff_file
    params:
        mito_chrom_id=config["mito_chrom_id"]
    shell:
        """
        zcat -f {input:q} | \
        grep -E '(^##[^#]|^{params.mito_chrom_id})' >\
        {output:q}
        """

rule gunzip:
    input:
        "{file}.gz"
    output:
        # This crazy regex actually means the output can't end in `.gz`
        "{file,^((?!\.gz).)*$}"
    shell:
        "gunzip -c {input:q} > {output:q}"

# rule bgzip_gff:
#     input:
#         os.path.join(gff_dir, gff_file_ungz)
#     output:
#         os.path.join(gff_dir, gff_basename_noext + ".sorted.gff.gz")
#     conda:
#         "../envs/tabix.yml"
#     shell:
#         """
#         (grep ^"#" {input:q}; grep -v ^"#" {input:q} | sort -k1,1 -k4,4n) | bgzip > {output:q}
#         """

# rule tabix_gff:
#     input:
#         os.path.join(gff_dir, gff_basename_noext + ".sorted.gff.gz")
#     output:
#         os.path.join(gff_dir, gff_basename_noext + ".sorted.gff.gz.tbi")
#     conda:
#         "../envs/tabix.yml"
#     shell:
#         "tabix -p gff {input:q}"

rule mito_gff_file_utrs:
    input:
        mito_gff_file
    output:
        mito_gff_utr_file
    conda:
        "../envs/plastid.yml"
    script:
        "../scripts/simulate_utrs.py"

rule nd4_gff_file:
    input:
        mito_gff_utr_file
    output:
        nd4_base + ".gff",
    shell:
        """
        grep -E '(^##[^#]|{config[nd4_gene_id]})' {input:q} > {output:q}
        """

rule nd4_roi_file:
    input:
        nd4_base + ".gff",
    output:
        bed=nd4_base + "_rois.bed",
        txt=nd4_base + "_rois.txt"
    params:
        outbase=nd4_base
    conda:
        "../envs/plastid.yml"
    shell:
        """
        metagene generate \
                --landmark cds_start \
                --annotation_files {input:q} \
                --annotation_format GFF3 \
                {params.outbase:q}
        """

rule nd6_gff_file:
    input:
        mito_gff_utr_file
    output:
        nd6_base + ".gff",
    shell:
        """
        grep -E '(^##[^#]|{config[nd6_gene_id]})' {input:q} > {output:q}
        """

rule nd6_roi_file:
    input:
        nd6_base + ".gff",
    output:
        bed=nd6_base + "_rois.bed",
        txt=nd6_base + "_rois.txt"
    params:
        outbase=nd6_base
    conda:
        "../envs/plastid.yml"
    shell:
        """
        metagene generate \
                --landmark cds_stop \
                --annotation_files {input:q} \
                --annotation_format GFF3 \
                {params.outbase:q}
        """

# rule biotype_table:
#     input:
#         config["genome_annotation_gff3_file"]
#     output:
#         os.path.join(gff_dir, "{}.gene_biotype_table.txt".format(gff_basename_noext))
#     shell:
#         '''
#         echo "#Gene biotype table" > {output:q}
#         printf 'Geneid\\tgene_biotype\\n' >> {output:q}
#         zcat -f {input:q} | \
#             sed -rn 's/.*ID=gene:([^;]*);.*biotype=([^;]*).*/\\1\\t\\2/p' >> \
#             {output:q}
#         '''
#             # sed -rn 's/([^\\t]*).*ID=gene:([^;]*);.*biotype=([^;]*).*/\\2\\t\\3\\t\\1/p' >> \

rule gene_annotation_table:
    input:
        config["genome_annotation_gtf_file"]
    output:
        gene_annotation_table
    shell:
        """
        zcat -f {input:q} | \
            awk 'BEGIN{{FS="\\t"}}{{split($9,a,";"); if($3~"gene") print a[1]"\\t"a[3]"\\t"$1":"$4"-"$5"\\t"a[5]"\\t"$7}}' | \
            sed 's/gene_id "//' | \
            sed 's/gene_id "//' | \
            sed 's/gene_biotype "//' | \
            sed 's/gene_name "//' | \
            sed 's/gene_biotype "//' | \
            sed 's/"//g' | \
            sed 's/ //g' | \
            sed '1igene_id\\tGeneSymbol\\tChromosome\\tClass\\tStrand' > \
            {output:q}
        """

rule gene_bed_file:
    input:
        gtf=config["genome_annotation_gtf_file"],
        fai=genome_fasta_unzipped + ".fai"
    output:
        genes_bed_file
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        zcat -f {input.gtf:q} | \
            awk 'BEGIN{{FS="\\t"}}{{split($9,a,";"); if($3~"gene") print $1"\\t"$4-1"\\t"$5"\\t"a[1]"\\t.\\t"$7}}' | \
            sed 's/gene_id "//' | \
            sed 's/"//g' | \
            bedtools sort -faidx {input.fai:q} > \
            {output:q}
        """

rule non_gene_bed_file:
    input:
        genes=genes_bed_file,
        fai=genome_fasta_unzipped + ".fai"
    output:
        nongenes_bed_file
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools complement \
            -g {input.fai:q} \
            -i {input.genes:q} > \
            {output:q}
        """

