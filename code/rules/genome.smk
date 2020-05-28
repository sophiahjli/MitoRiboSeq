from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule get_fasta_ftp:
    input:
        FTP.remote(config["genome_fasta_ftp"], keep_local=True)
    output:
        "{genome}" + genome_ext
    shell:
        "mv {input:q} {output:q}"

rule get_gff_ftp:
    input:
        FTP.remote(config["genome_annotation_gff3_ftp"], keep_local=True)
    output:
        config["genome_annotation_gff3_file"]
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

