# MitoRiboSeq
Codes for processing mitochondrial ribosome profiling data

The workflow for the data processing from sequenced data is:
1. Alignment and extract mitochondria-specific data
2. Metagene analysis to assign the mitoribosomal A-site position from the mapped reads
3. Generate codon table at nucleotide resolution
4. Codon occupancy analysis and gene-specific cumulative count visualization

### 1. Alignment and extract mitochondria-specific data
This is done using several tools and the workflow is available via GALAXY

### 2. Metagene analysis to assign the mitoribosomal A-site position from the mapped reads

## Install python in Anaconda/Miniconda
(Lance please add something here)

## Install R and Rstudio
Some subsequent analysis after codon count table generation is done in the R environment including read depth and coverage analysis, and plots for codon occupancy and gene-specific cumulative occupancy along the transcript.

## Create an environment for sequencing analysis
Here we use conda to install an environment named riboseq that has been defined in the `riboseq_environment.yml` file.
```
conda env creat -f riboseq_environment.yml
```
Now enter the environment by

```
conda activate riboseq
```

## Storing BAM files for analysis
BAM files are recommended to be stored in `data/mapping`. After each steps, corresponding data will be stored in subfolders under `data` for clean organization.

## Specifying specific parameters in the configuration file
`mito_config.yml` contains specific parameters that are experiment/analysis specific.  

## Analyzing the sequencing data
The analysis begins with BAM files. We use the python package `plastid` to facilitate our analysis. There are two snakefiles already made to iteratively call different functions from plastid. The user will only need to specify the setting in the `mito_config.yml` file. According to the 

## Metagene and readphasing analysis
This part wraps around the `plastid` functions including `phase_by_size` and `metagene count`. Simply execute the pipeline by typing:
```
snakemake -f mito_readphasing_metagene.snakefile
```
Note you can specify the samples and the range of RNA fragments to be analyzed in `mito_config.yml`.

## 3. Generate codon table at nucleotide resolution

Once the offset is determined, 


