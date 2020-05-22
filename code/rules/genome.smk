from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule get_fasta_ftp:
    input:
        FTP.remote(config["genome_fasta_ftp"], keep_local=True)
    output:
        "{genome}" + genome_ext
    shell:
        "mv {input:q} {output:q}"

