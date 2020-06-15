#shell.executable("bash")

snakemake.utils.min_version("5.15.0")

def splitext_gz(path):
    if os.path.splitext(path)[1] == ".gz":
        return ([os.path.splitext(os.path.splitext(path)[0])[0],
                 os.path.splitext(os.path.splitext(path)[0])[1],
                 ".gz"])
    else:
        return ([os.path.splitext(path)[0],
                 os.path.splitext(path)[1],
                 None])


configfile: "code/mito_config.defaults.yml"
sample_files = snakemake.utils.listfiles(config["fastq_file_pattern"])
samples = dict((y[0], x) for x, y in sample_files)
assert len(samples) > 0, "ERROR: No fastq files were found using pattern '{}' (set in configfile)".format(config["fastq_file_pattern"])

log_dir = config["results_dir"] + "/logs"

genome_dir = os.path.dirname(config["genome_fasta_file"])
genome_fasta = os.path.basename(config["genome_fasta_file"])
genome_ext = splitext_gz(genome_fasta)[1]
genome = splitext_gz(genome_fasta)[0]

if genome_fasta.endswith(".gz"):
    genome_fasta_unzipped = os.path.splitext(config["genome_fasta_file"])[0]
else:
    genome_fasta_unzipped = config["genome_fasta_file"]

gff_dir = os.path.dirname(config["genome_annotation_gff3_file"])
gff_basename = os.path.basename(config["genome_annotation_gff3_file"])
gff_file_gz = splitext_gz(gff_basename)[0] + splitext_gz(gff_basename)[1] + ".gz"
gff_file_ungz = splitext_gz(gff_basename)[0] + splitext_gz(gff_basename)[1]
gff_basename_noext = splitext_gz(gff_basename)[0]


gtf_dir = os.path.dirname(config["genome_annotation_gtf_file"])
gtf_basename = os.path.basename(config["genome_annotation_gtf_file"])
gtf_basename_noext = splitext_gz(gtf_basename)[0]
genes_bed_file = os.path.join(gtf_dir, "{}.genes.bed".format(gtf_basename_noext))
nongenes_bed_file = os.path.join(gtf_dir, "{}.nongenes.bed".format(gtf_basename_noext))
gene_annotation_table = os.path.join(gtf_dir, "{}.gene_annotation_table.txt".format(gtf_basename_noext))

mito_gff_file = os.path.join(gff_dir, "{}.mito.gff".format(gff_basename_noext))
mito_gff_utr_file = os.path.join(gff_dir, "{}.mito.added_utrs.gff".format(gff_basename_noext))
nd4_base = os.path.join(gff_dir, "{}.nd4".format(gff_basename_noext))
nd6_base = os.path.join(gff_dir, "{}.nd6".format(gff_basename_noext))

def get_fastq(wildcards):
    return samples[wildcards.sample]
