#!/usr/bin/env Rscript --vanilla

library(tidyverse, quietly=TRUE)
library(ggrepel, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("Generate codon occupancy tables and figures")

p <- add_argument(p, "codon_counts", help="Codon count data table")
p <- add_argument(p, "gene_annotations", help="Mapping of gene IDs to gene names")
p <- add_argument(p, "mitochondria_codon_table", help="Mitochondria codon table")

# Add an optional arguments
p <- add_argument(p, "--raw_coverage_bygene",
                  help="Output filename for table of raw codon coverage by gene",
                  default="mito_data_raw_coverage_bygene.csv")
p <- add_argument(p, "--raw_depth_bygene",
                  help="Output filename for table of raw codon depth by gene",
                  default="mito_data_raw_depth_bygene.csv")
p <- add_argument(p, "--occupancy_table",
                  help="Output filename for table of mitochondrial codon occupancy",
                  default="mito_occupancy_table.csv")
p <- add_argument(p, "--cumsum_table",
                  help="Output filename for table of cumulative sum of mitochondrial codon occupancy",
                  default="mito_cumsum_table.csv")
p <- add_argument(p, "--figures_directory",
                  help="Directory where figures will be saved",
                  default="figures")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    snakemake@input[["codon_counts"]],
    snakemake@input[["gene_annotations"]],
    snakemake@input[["mitochondria_codon_table"]],
    "--raw_coverage_bygene", snakemake@output[["raw_coverage_bygene"]],
    "--raw_depth_bygene", snakemake@output[["raw_depth_bygene"]],
    "--occupancy_table", snakemake@output[["occupancy_table"]],
    "--cumsum_table", snakemake@output[["cumsum_table"]],
    "--figures_directory", snakemake@params[["figures_directory"]]
  ))
} else if (interactive()) {
    # Arguments supplied inline (for debug/testing when running interactively)
    print("Running interactively...")
    codon_counts_file <- "results/codon_count/All_codoncount_table.txt"
    gene_annotations_file <- "genome/Homo_sapiens.GRCh38.100.gene_annotation_table.txt"
    mitochondria_codon_table <- "genome/AA_Codon_HumMito.csv"
    argv <- parse_args(p, c(codon_counts_file,
                            gene_annotations_file,
                            mitochondria_codon_table))
    print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)  
}


############
#Load data #
############
all_raw_data <- read_tsv(argv$codon_counts)

# Load the mapping where we can turn IDs into names that are easily identifiable
mito_info <-read_tsv(argv$gene_annotations)

raw_data <- all_raw_data %>%
  inner_join(., mito_info, by="gene_id")

# Load the mitochondria codon table
mito_aminoacid_codon <- read_csv(argv$mitochondria_codon_table)


####################################
# Read depth and coverage analysis #
####################################

# use raw data to analyze the quality of riboseq data
# coverage is defined as the percentage of codons that have at least one mapped read per gene
# depth is defined as the average number of mapped reads per gene

mito_data_raw_coverage <- raw_data %>% group_by(sample,GeneSymbol) %>%
  summarise(coverage = sum(codon_count_sum> 0)/n(), depth = sum(codon_count_sum)/n()) 

mito_data_raw_coverage_bygene <- mito_data_raw_coverage %>%
  select(sample, GeneSymbol,coverage) %>%  spread(GeneSymbol,coverage)
write_csv(mito_data_raw_coverage_bygene, path = argv$raw_coverage_bygene)

mito_data_raw_depth_bygene <- mito_data_raw_coverage %>%
  select(sample, GeneSymbol,depth) %>% spread(GeneSymbol,depth)
write_csv(mito_data_raw_depth_bygene, path = argv$raw_depth_bygene)

##########################################
# Cumulative sum of codon count analysis #
##########################################

# Normalize counts in each sample to RPM
mito_data <- raw_data %>% select(-contains("position")) %>% group_by(sample) %>%
  mutate(RPM = codon_count_sum/sum(codon_count_sum, na.rm = TRUE) * 10^6 + 1) %>%
  group_by(sample,GeneSymbol) %>%
  mutate(RPM_normgene = RPM/sum(RPM, na.rm = TRUE),
         RPM_cumsum = cumsum(RPM),
         RPM_cumsum_normgene = RPM_cumsum/sum(RPM, na.rm = TRUE)) %>% mutate(codon_seq = toupper(codon_seq))

mito_data_cumsum_plot <- mito_data %>% ggplot(aes(x=codon_index,y=RPM_cumsum_normgene,col=sample)) +
  geom_line() +
  scale_y_continuous(labels = function(n) format(n,digits=2,scientific=T)) +
  theme(strip.background = element_blank(), aspect.ratio = 0.8) +
  facet_wrap(~GeneSymbol, scales = "free") +
  facet_wrap(~GeneSymbol, scales = "free_x") +
  ylab("Normalized cumulative ribosome occupancy")

if (! dir.exists(argv$figures_directory)) {
  dir.create(argv$figures_directory)
}
save_plot(filename = file.path(argv$figures_directory, "norm_cumsum_plot.pdf"),
          mito_data_cumsum_plot, base_height = 8, base_width = 12)
write_csv(x = mito_data, path = argv$cumsum_table)

#############################
# Codon occupancy analysis #
#############################


calc_codon_occupancy <- function(data){
  # Two methods are used to calculate codon occupancy in this function
  # Occupancy is defined as the ratio between actual counts to expected counts, which is dependent on the codon frequency
  
  # This one takes the occupancy for codons in each gene first
  occupancy_bygene <- data %>% group_by(GeneSymbol) %>% mutate(gene_length = n()) %>% 
    group_by(GeneSymbol,codon_seq,gene_length) %>% summarise(codon_num = n(),RPM_normgene = sum(RPM_normgene)) %>% 
    mutate(codon_freq = codon_num/gene_length,occupancy = RPM_normgene/codon_freq)
  # Then use the mean of occupancy from all genes to calculate the final occupancy
  occupancy_bygene_summary <- occupancy_bygene %>% group_by(codon_seq) %>% summarise(occupancy_bygene = mean(occupancy))
  
  # This one directly takes all the codons from all genes and calculate the occupancy
  occupancy_bygenome <- data %>% ungroup() %>% mutate(genome_size = n()) %>% 
    group_by(codon_seq,genome_size) %>% summarise(codon_num = n(),RPM = sum(RPM)) %>% 
    mutate(codon_freq = codon_num/genome_size,occupancy_bygenome = RPM/codon_freq/10^6)
  occupancy_data <- occupancy_bygenome %>% select(codon_seq,occupancy_bygenome) %>% inner_join(.,occupancy_bygene_summary)
  return(occupancy_data)
}  

plot_occupancy <- function(sample,occupancy){
  # This function defines how to plot the occupancy
  Met <- c("ATG","ATA")
  Stop <- c("TAA","TAG","AGA","AGG")
  # Label four different styles of plotting depending on what type of codon it is
  occupancy <- occupancy %>% mutate(sig = case_when(occupancy_bygenome >= mean(.$occupancy_bygenome) + 1.5*sd(.$occupancy_bygenome)~ 1,TRUE ~ 0)) %>%
    mutate(is_Met = codon_seq %in% Met, is_stop = codon_seq %in% Stop,
           plot_col = case_when(is_Met == T ~ "2",is_stop ==T ~ "3",sig ==1 ~ "4",T~"1")) %>%
    mutate(codon_type = case_when(codon_seq %in% Met ~ "Start",
                                  codon_seq %in% Stop ~ "Stop",
                                  sig == 1 ~ "Sig",
                                  T ~ "Codons"))
  
  plot_style <- list(geom_point(),
                     geom_hline(yintercept = 1,linetype = 2),
                     theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                           legend.position = "right",legend.title = element_blank(),aspect.ratio = 0.442),
                     xlab(""),ggtitle(sample),
                     scale_colour_manual(values = c("Codons" = "#000000",
                                                    "Start" = "#0021fa",
                                                    "Stop" = "#009e25",
                                                    "Sig" = "#fa0000")))
  
  output_plot <- ggplot(data = occupancy, aes(x=reorder(codon_seq,occupancy_bygenome),y = occupancy_bygenome, col = factor(codon_type))) +
    plot_style +  geom_text_repel(data = occupancy %>% filter(sig>0),aes(label=codon_seq),show.legend = FALSE) +
    ylab("Codon occupancy") 
  
  return(output_plot)
}
# Group and nest the data to calculate occupancy
mito_occupancy <- mito_data %>% group_by(sample) %>% nest() %>% mutate(occupancy = map(data,~calc_codon_occupancy(.))) %>% 
  select(-data) 
# Take the occupancy to plots
mito_occupancy_plot <- mito_occupancy %>% mutate(occu_plot = map2(sample,occupancy,~plot_occupancy(.x,.y)))

# Write the occupancy plot of each sample to PDF
paths <- stringr::str_c(argv$figures_directory, "/", mito_occupancy$sample, "_occupancy_plot.pdf")
pwalk(list(paths,mito_occupancy_plot$occu_plot), ggsave)

# Save the occupancy table to a CSV file
mito_occupancy_tab <- mito_occupancy %>% unnest(cols = c(occupancy))
write_csv(x = mito_occupancy_tab, path = argv$occupancy_table)
