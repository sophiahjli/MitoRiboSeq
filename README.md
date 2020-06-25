# MitoRiboSeq
This workflow is used to generate codon- and gene-based analysis from mitochondrial ribosome profiling data generated using the protocols described in Monitoring mitochondrial translation with ribosome profiling in Nature Protocol by Li et al. 2019 [ref link].

Dependencies are installed using [Bioconda](https://bioconda.github.io/).
The workflow is written using [Snakemake](https://snakemake.readthedocs.io/).


## Overview

This workflow is designed to take FASTQ files from the Illumina sequencer to map each mitoribosome footprint to its occupied A-site. This mapping will allow detailed monitor of mitoribosome translation dynamics. 

Starting with FASTQ files, the workflow is divided grossly into three parts: QC and metagene analysis, A-site assignment and codon count table generation, and downstream codon occupancy and cumulative ribosome occupancy analysis. In the first part, the workflow will run a QC, align all the reads to the genome, and do metagene gene analysis to produce a file where the user can determine the best offset to assign the ribosomal A-site. Once the offset is determined, all the reads will be assigned to the A-site at nucleotide level and a codon count table where the counts for each codon in each gene are arranged in a table to facilitate downstream analysis. Two common analysis we routinely used such as codon occupancy analysis and cumulative mitoribosome footprint along the transcript will be provided to visaulize mitoribosome distribution on mitochondria-encoded genes.

### Inputs

*   Reference sequence in FASTA format
*   FASTQ files for each sample
*   ROI files `metagene_roi_nd4_file` and `metagene_roi_nd6_file` created using the 
    [plastid `metagene generate`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.metagene.html#module-plastid.bin.metagene) command
*   Gene model in GTF format
*   Configuration file(s) in YAML format

### Outputs

*   `mapping` - BAM alignment files 
*   `metagene` - [plastid `metagene count`] output for both the `metagene_roi_nd4_file` and the `metagene_roi_nd6_file` specified in the config
*   `phasing_analysis` - [plastid `phaze_by_size`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.phase_by_size.html#module-plastid.bin.phase_by_size)
    output to estimate sub-codon phasing, stratified by read length.
*   `wiggle` - [plastid `make_wiggle`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.make_wiggle.html#module-plastid.bin.make_wiggle) output. Genome browser tracks from read alignments, using mapping rules to extract ribosomal P-sites from the alignments.
*   `codon_count` - It contains the codon count table of all samples, separated and combined as one file.
*   `metagene` - It contains the metagene analysis result from the plastid metagene program.
*   `QC` - It contains all the QC of analysis, including read depth and coverage.
*   `figures` - It contains all the figures.
*   `tables` - It contains the cumulative codon counts for each gene and codon occupancy analysis results.

#### Intermediate outputs

*    `trimmed` - Fastq files that have been trimmed of adapter and low quality sequences

### Workflow

1.  **FASTQ summary and QC metrics** - Use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to determine some basic QC metrics from the raw FASTQ files
2.  **Trim reads** - Trim adapter sequences and low quality bases from fastq files using [cutadapt](https://cutadapt.readthedocs.io/en/stable/).
3.  **Align reads** - Use [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) to align reads to the genome
4.  **Metagene analysis** - Use the [plastid](https://plastid.readthedocs.io/en/latest/) scripts 
    [`metagene`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.metagene.html#module-plastid.bin.metagene) to provide a profile of counts over stop codons
5.  **Read Phasing Analysis** - Use the [plastid](https://plastid.readthedocs.io/en/latest/) script
    [`phase_by_size`](https://plastid.readthedocs.io/en/latest/generated/plastid.bin.phase_by_size.html#module-plastid.bin.phase_by_size) 
    to estimate [sub-codon phasing](https://plastid.readthedocs.io/en/latest/glossary.html#term-sub-codon-phasing), stratified by read length
6.  **Visualize** - visalize the ribosomal P-site coverage by generating genome browser tracks (wiggle format)
7.  **Codon Analysis** - Use [plastid](https://plastid.readthedocs.io/en/latest/) to generate ribosomal occupancy counts per codon
8.  **Read Length Distribution** - *Lorem ipsum dolor sit amet, consectetur adipiscing elit. Proin auctor.*
9.  **Codon Occupancy** - Written in R, *lorem ipsum dolor sit amet, consectetur adipiscing elit. Mauris eleifend.*


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

    1.  *FASTQ* files - the FASTQ data from the sequencer should
        be stored in `data/fastq` in `fastq.gz` format, one file
        per sample.
    2.  *Genome* files
        1. The genome reference sequence in fasta format
        2. The gene annotations in GFF v3 format
        3. Mitochondrial gene info file  # **TODO** Describe format
        4. Mitochondrial amino acid codon table  # **TODO** Describe format
        5. Mitochondrial orfs start rois ND4  # **TODO** Describe format
        6. Mitochondrial orfs stop rois ND6  # **TODO** Describe format

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
    snakemake --configfile "code/mito_config.yml" --use-conda -s code/mito_readphasing_metagene.snakefile
    ```
    
    If `--configfile` is not specified, the defaults are used.

7.  Execute the codon occupancy workflow

    ```bash
    snakemake --configfile "code/mito_config.yml" --use-conda -s code/mito_codontable.snakefile
    ```

## Common snakemake options

*   `--configfile "myconfig.yml"` - Override defaults using the configuration found in `myconfig.yml`
*   `--use-conda` - Use [conda]() to create an environment for each rule, installing and using the exact version of the software required (recommended)
*   `--cores` - Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores.
*   `--cluster` - Execute snakemake rules with the given submit command, *e.g.* `qsub`. Snakemake compiles jobs into scripts that are submitted to the cluster with the given command, once all input files for a particular job are present. The submit command can be decorated to make it aware of certain job properties (input, output, params, wildcards, log, threads and dependencies (see the argument below)), *e.g.*: `$ snakemake –cluster ‘qsub -pe threaded {threads}’`.
*   `--drmaa` - Execute snakemake on a cluster accessed via [DRMAA](https://en.wikipedia.org/wiki/DRMAA). Snakemake compiles jobs into scripts that are submitted to the cluster with the given command, once all input files for a particular job are present. `ARGS` can be used to specify options of the underlying cluster system, thereby using the job properties input, output, params, wildcards, log, threads and dependencies, *e.g.*: `--drmaa ‘ -pe threaded {threads}’`. Note that `ARGS` must be given in quotes and with a leading whitespace.
*   `--cluster-config` - A JSON or YAML file that defines the wildcards used in `cluster` for specific rules, instead of having them specified in the Snakefile. For example, for rule `job` you may define: `{ ‘job’ : { ‘time’ : ‘24:00:00’ } }` to specify the time for rule `job`. You can specify more than one file. The configuration files are merged with later values overriding earlier ones.
*   `--dryrun` - Do not execute anything, and display what would be done. If you have a very large workflow, use `--dryrun --quiet` to just print a summary of the DAG of jobs.

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
