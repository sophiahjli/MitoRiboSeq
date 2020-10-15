# MitoRiboSeq

This workflow is used to generate codon- and gene-based analysis from mitochondrial ribosome profiling data generated using the protocols described in Monitoring mitochondrial translation with ribosome profiling in Nature Protocol by Li et al. 2020.

Dependencies are installed using [Bioconda](https://bioconda.github.io/).
The workflow is written using [Snakemake](https://snakemake.readthedocs.io/).


## Overview

This workflow is designed to take FASTQ files from the Illumina sequencer and map each mitoribosome footprint to its occupied A-site. This mapping will allow detailed monitoring of mitoribosome translation dynamics. 

Starting with FASTQ files, the workflow is divided into three main parts: QC and metagene analysis, A-site assignment and codon count table generation, and downstream codon occupancy analysis. In the first part, the workflow will run QC on the input FASTQ files, then align all of the reads to the genome (nuclear and micorhondrial genome together), and finally perform a metagene analysis. The metagene analysis helps the user to determine the offset from the 3' end of the read to the ribosomal A-site. Once the offset is determined, all of the reads will be assigned to the A-site at the nucleotide level and the counts for each codon in each gene are arranged in a table to facilitate downstream analysis. Once the codon count table is generated, two common analyses, which are presented in Figure 7 of the manuscript, It includes codon occupancy analysis and cumulative mitoribosome footprint along the transcripts, are provided to visaulize mitoribosome distribution on mitochondria-encoded genes.

### Inputs

*   FASTQ files for each sample
*   Configuration file(s) in YAML format
*   Reference sequence in FASTA format (by default this is obtained from Ensembl)
*   Gene model in GTF format (by default this is obtained from Ensembl)

### Outputs

*   `mapped` - BAM alignment files 
*   `metagene` - [plastid `metagene count`] output used to determine the A-site (using stop codon) and P-site (using start codon) offsets 
*   `codon_count` - codon count table of all samples, separated and combined as one file
*   `phasing_analysis` - [plastid `phaze_by_size`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.phase_by_size.html#module-plastid.bin.phase_by_size)
    output to estimate sub-codon phasing, stratified by read length.
*   `bedgraph` - [plastid `make_wiggle`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.make_wiggle.html#module-plastid.bin.make_wiggle) output. Genome browser tracks from read alignments, using mapping rules to extract ribosomal A-sites from the alignments
*   `qc` - the quality control analysis, including read depth and coverage
*   `figures` - all the figures
*   `tables` - cumulative codon counts for each gene and codon occupancy analysis results
*   `logs` - logs from each of the workflow steps, used in troubleshooting

#### Intermediate outputs

*    `trimmed` - Fastq files that have been trimmed of adapter and low quality sequences

### Workflow

1.  **FASTQ summary and QC metrics** - Use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to determine some basic QC metrics from the raw FASTQ files
2.  **Trim reads** - Trim adapter sequences and low quality bases from fastq files using [cutadapt](https://cutadapt.readthedocs.io/en/stable/).
3.  **Align reads** - Use [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) to align reads to the genome
4.  **Metagene analysis** - Use the [plastid](https://plastid.readthedocs.io/en/latest/) scripts 
    [`metagene`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.metagene.html#module-plastid.bin.metagene) to provide a profile of counts relative to start and stop codons to determine the A-site and P-site offset for the experiment
5.  **Read Phasing Analysis** - Use the [plastid](https://plastid.readthedocs.io/en/latest/) script
    [`phase_by_size`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.phase_by_size.html#module-plastid.bin.phase_by_size) 
    to estimate [sub-codon phasing](https://plastid.readthedocs.io/en/latest/glossary.html#term-sub-codon-phasing), stratified by read length
6.  **Generate codon count table** - Use [plastid](https://plastid.readthedocs.io/en/latest/) to generate codon count table with the offset determined in the Metagene analysis (assign each reads to the corresponding ribosomal A-site)
7.  **Visualize** - Visualize the ribosomal A-site coverage by generating genome browser tracks (wiggle format)
8.  **Codon occupancy Analysis** - Written in R to generate three figures. The first one uses hierarchical clustering and heatmap to compare samples by their codon occupancy. The second one shows the codon occupancy result for each sample individually. The third one identifies ribosome stalling site along the transcripts by plotting the cumulative ribosome counts.


## Install prerequisites

1. Install conda

    *   If you have [Anaconda](https://www.anaconda.com/distribution/) installed, you already have it.
    *   Otherwise, install the [Miniconda](https://conda.io/en/latest/miniconda.html) package.

2.  Enable the [Bioconda](https://bioconda.github.io/#using-bioconda) channel
    (requires  64-bit Linux or Mac OS, Windows is not supported)

    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

## Setup environment and run workflow

1.  Clone workflow into working directory

    ```bash
    git clone https://github.com/sophiahjli/MitoRiboSeq.git
    cd MitoRiboSeq
    ```

2.  Input data

    *FASTQ* files - the FASTQ data from the sequencer should
    be stored in `data/fastq` in `fastq.gz` format, one file
    per sample.

3.  Edit configuration files as needed

    ```bash
    cp code/mito_config.defaults.yml code/mito_config.yml
    nano code/mito_config.yml
    
    # Only if running on a cluster
    cp cluster_config.yml mycluster_config.yml
    nano mycluster_config.yml
    ```

4.  Install dependencies into an isolated environment

    ```bash
    conda env create --file code/mitoriboseq_environment.yml
    ```

    Note: If you are updating the workflow, you may need to update the conda environment
    ```bash
    conda env update --file code/mitoriboseq_environment.yml
    ```

5.  Activate the environment

    ```bash
    source activate MitoRiboSeq
    ```

6.  Execute the trimming, mapping, read phasing, and metagene workflow

    ```bash
    snakemake \
         -s code/mito_readphasing_metagene.snakefile \
         --configfile "code/mito_config.yml" \
         --use-conda \
         --cores 4
    ```
    
    If `--configfile` is not specified, the defaults are used.

7.  Execute the codon occupancy workflow

    ```bash
    snakemake \
         -s code/mito_codontable.snakefile \
         --configfile "code/mito_config.yml" \
         --use-conda \
         --cores 4
    ```

## Common snakemake options


*   `--cores` -  Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the number of available CPU cores. **Required**
*   `--configfile "myconfig.yml"` - Override defaults using the configuration found in `myconfig.yml`
*   `--dryrun` - Do not execute anything, and display what would be done. If you have a very large workflow, use `--dryrun --quiet` to just print a summary of the DAG of jobs.
*   `--use-conda` - Use [conda](http://conda.io) to create an environment for each rule, installing and using the exact version of the software required (*recommended*)
*   `--cluster` - Execute snakemake rules with the given submit command, *e.g.* `qsub`. Snakemake compiles jobs into scripts that are submitted to the cluster with the given command, once all input files for a particular job are present. The submit command can be decorated to make it aware of certain job properties (input, output, params, wildcards, log, threads and dependencies (see the argument below)), *e.g.*: `$ snakemake –cluster ‘qsub -pe threaded {threads}’`.
*   `--cluster-config` - A JSON or YAML file that defines the wildcards used in `cluster` for specific rules, instead of having them specified in the Snakefile. For example, for rule `job` you may define: `{ ‘job’ : { ‘time’ : ‘24:00:00’ } }` to specify the time for rule `job`. You can specify more than one file. The configuration files are merged with later values overriding earlier ones.
*   `--set-threads [RULE=THREADS [RULE=THREADS ...]]` -  Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule. (default: None)
*   `--drmaa` - Execute snakemake on a cluster accessed via [DRMAA](https://en.wikipedia.org/wiki/DRMAA). Snakemake compiles jobs into scripts that are submitted to the cluster with the given command, once all input files for a particular job are present. `ARGS` can be used to specify options of the underlying cluster system, thereby using the job properties input, output, params, wildcards, log, threads and dependencies, *e.g.*: `--drmaa ‘ -pe threaded {threads}’`. Note that `ARGS` must be given in quotes and with a leading whitespace.

See the [Snakemake documentation for a list of all options](https://snakemake.readthedocs.io/en/stable/executable.html#all-options).


## Snakemake Examples 

### Running workflow on a single computer with 4 threads

```bash
snakemake --configfile "myconfig.yml" --use-conda --cores 4
``` 

### Running workflow using SLURM

```bash
snakemake \
    --configfile "myconfig.yml" \
    --cluster-config "mycluster_config.yml" \
    --cluster "sbatch --cpus-per-task={cluster.n} --mem={cluster.memory} --time={cluster.time}" \
    --use-conda \
    --cores 100
``` 

### Running workflow on Princeton LSI cluster using [DRMAA](https://en.wikipedia.org/wiki/DRMAA)

```bash
snakemake \
    --configfile "myconfig.yml" \
    --cluster-config "cetus_cluster.yml" \
    --drmaa " --cpus-per-task={cluster.n} --mem={cluster.memory} --qos={cluster.qos} --time={cluster.time}" \
    --use-conda \
    --cores 1000 \
    --output-wait 60
```
