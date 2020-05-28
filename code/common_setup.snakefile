#shell.executable("bash")

snakemake.utils.min_version("5.15.0")

def splitext_gz(path):
    if os.path.splitext(path)[1] == ".gz":
        return ([os.path.splitext(os.path.splitext(path)[0])[0],
                 os.path.splitext(os.path.splitext(path)[0])[1] + ".gz"])
    else:
        return os.path.splitext(path)


configfile: "code/mito_config.defaults.yml"
sample_files = snakemake.utils.listfiles(config["fastq_file_pattern"])
samples = dict((y[0], x) for x, y in sample_files)
assert len(samples) > 0, "ERROR: No fastq files were found using pattern '{}' (set in configfile)".format(config["fastq_file_pattern"])

log_dir = config["results_dir"] + "/logs"

genome_dir = os.path.dirname(config["genome_fasta_file"])
genome_fasta = os.path.basename(config["genome_fasta_file"])
genome_ext = splitext_gz(genome_fasta)[1]
genome = splitext_gz(genome_fasta)[0]

gff_dir = os.path.dirname(config["genome_annotation_gff3_file"])
gff_file = os.path.basename(config["genome_annotation_gff3_file"])
gff_ext = splitext_gz(gff_file)[1]
gff = splitext_gz(gff_file)[0]

mito_gff_file = os.path.join(gff_dir, "{}.mito.gff".format(gff))
mito_gff_utr_file = os.path.join(gff_dir, "{}.mito.added_utrs.gff".format(gff))
nd4_base = os.path.join(gff_dir, "{}.nd4".format(gff))
nd6_base = os.path.join(gff_dir, "{}.nd6".format(gff))

def get_fastq(wildcards):
    return samples[wildcards.sample]
