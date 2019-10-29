library(tidyverse)
library(ggrepel)
library(cowplot)

############
#Load data #
############
all_raw_data <- read_tsv(snakemake@input[["all_codon_count"]])

# Load the mapping where we can turn IDs into names that are easily identifiable
mito_info <-read_tsv(snakemake@input[["mito_info"]])

raw_data <- all_raw_data %>%
  inner_join(., mito_info, by=c("gene_id"="UCSC_id")) %>%
  select(gene_name=AssociatedGeneName, GeneType, sample, codon_seq, codon_index, codon_count_sum, Complex)

# Load the mitochondria codon table
mito_aminoacid_codon <- read_csv(snakemake@input[["mito_aminoacid_codon"]])


####################################
# Read depth and coverage analysis #
####################################

# use raw data to analyze the quality of riboseq data
# coverage is defined as the percentage of codons that have at least one mapped read per gene
# depth is defined as the average number of mapped reads per gene

mito_data_raw_coverage <- raw_data %>% group_by(sample,gene_name) %>%
  summarise(coverage = sum(codon_count_sum> 0)/n(), depth = sum(codon_count_sum)/n()) 

mito_data_raw_coverage_bygene <- mito_data_raw_coverage %>%
  select(sample, gene_name,coverage) %>%  spread(gene_name,coverage)
write_csv(mito_data_raw_coverage_bygene, path = snakemake@output[["mito_data_raw_coverage_bygene"]])

mito_data_raw_depth_bygene <- mito_data_raw_coverage %>%
  select(sample, gene_name,depth) %>% spread(gene_name,depth)
write_csv(mito_data_raw_depth_bygene, path = snakemake@output[["mito_data_raw_depth_bygene"]])

##########################################
# Cumulative sum of codon count analysis #
##########################################

# Normalize counts in each sample to RPM
mito_data <- raw_data %>% select(-contains("position")) %>% group_by(sample) %>%
  mutate(RPM = codon_count_sum/sum(codon_count_sum, na.rm = TRUE) * 10^6 + 1) %>%
  group_by(sample,gene_name) %>%
  mutate(RPM_normgene = RPM/sum(RPM, na.rm = TRUE),
         RPM_cumsum = cumsum(RPM),
         RPM_cumsum_normgene = RPM_cumsum/sum(RPM, na.rm = TRUE)) %>% mutate(codon_seq = toupper(codon_seq))

mito_data_cumsum_plot <- mito_data %>% ggplot(aes(x=codon_index,y=RPM_cumsum_normgene,col=sample)) +
  geom_line() +
  scale_y_continuous(labels = function(n) format(n,digits=2,scientific=T)) +
  theme(strip.background = element_blank(), aspect.ratio = 0.8) +
  facet_wrap(~gene_name, scales = "free") +
  facet_wrap(~gene_name, scales = "free_x") +
  ylab("Normalized cumulative ribosome occupancy")

dir.create(snakemake@output[["figures"]])
save_plot(filename = file.path(snakemake@output[["figures"]], "norm_cumsum_plot.pdf"),
          mito_data_cumsum_plot, base_height = 8, base_width = 12)
write_csv(x = mito_data, path = snakemake@output[["mito_cumsum_table"]])

#############################
# Codon occupancy analysis #
#############################


calc_codon_occupancy <- function(data){
  # Two methods are used to calculate codon occupancy in this function
  # Occupancy is defined as the ratio between actual counts to expected counts, which is dependent on the codon frequency
  
  # This one takes the occupancy for codons in each gene first
  occupancy_bygene <- data %>% group_by(gene_name) %>% mutate(gene_length = n()) %>% 
    group_by(gene_name,codon_seq,gene_length) %>% summarise(codon_num = n(),RPM_normgene = sum(RPM_normgene)) %>% 
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
           plot_col = case_when(is_Met == T ~ "2",is_stop ==T ~ "3",sig ==1 ~ "4",T~"1"))
  
  plot_style <- list(geom_point(),
                     geom_hline(yintercept = 1,linetype = 2),
                     theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                           legend.position = "right",legend.title = element_blank(),aspect.ratio = 0.442),
                     xlab(""),ggtitle(sample),
                     scale_colour_manual(values = c("#000000","#0021fa","#009e25","#fa0000"),labels =c("Codons","Start","Stop","Sig")))
  
  output_plot <- ggplot(data = occupancy, aes(x=reorder(codon_seq,occupancy_bygenome),y = occupancy_bygenome,col = plot_col)) +
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
paths <- stringr::str_c(snakemake@output[["figures"]], "/", mito_occupancy$sample, "_occupancy_plot.pdf")
pwalk(list(paths,mito_occupancy_plot$occu_plot), ggsave)

# Save the occupancy table to a CSV file
mito_occupancy_tab <- mito_occupancy %>% unnest()
write_csv(x = mito_occupancy_tab, path = snakemake@output[["mito_occupancy_table"]])