shell.executable("bash")

snakemake.utils.min_version("5.4.2")

configfile: "code/mito_config.defaults.yml"
sample_files = snakemake.utils.listfiles(config["fastq_file_pattern"])
samples = dict((y[0], x) for x, y in sample_files)
print(samples)
assert len(samples) > 0, "ERROR: No fastq files were found using pattern '{}' (set in configfile)".format(config["fastq_file_pattern"])

log_dir = config["results_dir"] + "/logs"

genome_dir = os.path.dirname(config["genome_fasta_file"])
genome_fasta = os.path.basename(config["genome_fasta_file"])
genome_ext = os.path.splitext(genome_fasta)[1]
genome = os.path.splitext(genome_fasta)[0]


def get_fastq(wildcards):
    return samples[wildcards.sample]
