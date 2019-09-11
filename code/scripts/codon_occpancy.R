library(tidyverse)
library(ggrepel)
library(cowplot)

all_raw_data <- read_tsv(snakemake@input[["all_codon_count"]])

# Load the mapping where we can turn IDs into names that are easily identifiable
mito_info <-read_tsv(snakemake@input[["mito_info"]])

raw_data <- all_raw_data %>%
    left_join(., mito_info, by=c("gene_id"="UCSC_id")) %>%
    select(gene_name=AssociatedGeneName, GeneType, sample, codon_seq, codon_index, codon_count_sum, Complex)

# Load the mitochondria codon table
mito_aminoacid_codon <- read_csv(snakemake@input[["mito_aminoacid_codon"]])

condition_base <- "WT"

# use raw data to analyze the quality of riboseq data
# coverage is defined as the percentage of codons that have at least one mapped read per gene
# depth is defined as the average number of mapped reads per gene

mito_data_raw_coverage <- raw_data %>% group_by(sample,gene_name) %>%
  summarise(coverage = sum(codon_count_sum> 0)/n(), depth = sum(codon_count_sum)/n()) %>%
  separate(sample, into= c("condition","replicates"),remove = FALSE)

mito_data_raw_coverage_bygene <- mito_data_raw_coverage %>%
  select(sample, gene_name,coverage) %>%  spread(gene_name,coverage)
write_csv(mito_data_raw_coverage_bygene, path = snakemake@output[["mito_data_raw_coverage_bygene"]])

mito_data_raw_depth_bygene <- mito_data_raw_coverage %>%
  select(sample, gene_name,depth) %>% spread(gene_name,depth)
write_csv(mito_data_raw_depth_bygene, path = snakemake@output[["mito_data_raw_depth_bygene"]])

mito_data_raw_coverage_summary_bycondition <- mito_data_raw_coverage %>%
  group_by(condition) %>% summarise_at(vars(coverage),funs(mean,sd,min,max,median), na.rm = TRUE)
write_csv(mito_data_raw_coverage_summary_bycondition, path = snakemake@output[["mito_data_raw_coverage_summary_bycondition"]])

mito_data_raw_depth_summary_bycondition <- mito_data_raw_coverage %>%
  group_by(condition) %>% summarise_at(vars(depth),funs(mean,sd,min,max,median), na.rm = TRUE)
write_csv(mito_data_raw_depth_summary_bycondition, path = snakemake@output[["mito_data_raw_depth_summary_bycondition"]])



# Normalize counts in each sample to RPM
mito_data <- raw_data %>% select(-contains("position")) %>% group_by(sample) %>%
  mutate(RPM = codon_count_sum/sum(codon_count_sum, na.rm = TRUE) * 10^6 + 1) %>%
  group_by(sample,gene_name) %>%
  mutate(RPM_normgene = RPM/sum(RPM, na.rm = TRUE),
         RPM_cumsum = cumsum(RPM),
         RPM_cumsum_normgene = RPM_cumsum/sum(RPM, na.rm = TRUE)) %>%
  separate(sample,c("condition","replicate"),sep="-",remove = FALSE) %>%
  mutate_at(vars(codon_seq),funs(toupper))

# Summarise the replicates and generate cumsum.
mito_data_samplemean <- mito_data %>% group_by(condition,codon_index,gene_name,codon_seq) %>%
    summarise(RPM_samplemean = mean(RPM, na.rm = TRUE),
              RPM_samplemean_se = sd(RPM, na.rm = TRUE)/sqrt(n()),
              RPM_normgene_samplemean = mean(RPM_normgene, na.rm = TRUE),
              RPM_normgene_samplemean_se = sd(RPM_normgene, na.rm = TRUE)/sqrt(n())) %>%
    group_by(condition,gene_name) %>%
    mutate(RPM_cumsum_samplemean = cumsum(RPM_samplemean),
           RPM_cumsum_samplemean_se = cumsum(RPM_samplemean_se),
           RPM_cumsum_samplemean_normgene = cumsum(RPM_samplemean)/sum(RPM_samplemean, na.rm = TRUE),
           RPM_cumsum_samplemean_normgene_se = cumsum(RPM_samplemean_se)/sum(RPM_samplemean, na.rm = TRUE))


mito_data_cumsum <- mito_data  %>%
    group_by(condition,gene_name,codon_index) %>%
    summarise(avg_cumsum = mean(RPM_cumsum_normgene))

mito_data_cumsum_plot <- mito_data_cumsum %>% ggplot(aes(x=codon_index,y=avg_cumsum,col=condition)) +
  geom_line() +
  scale_y_continuous(labels = function(n) format(n,digits=2,scientific=T)) +
  theme(strip.background = element_blank(), aspect.ratio = 0.8) +
  facet_wrap(~gene_name, scales = "free") +
  facet_wrap(~gene_name, scales = "free_x") +
  ylab("Normalized cumulative ribosome occupancy")

dir.create(snakemake@output[["figures"]])
save_plot(filename = file.path(snakemake@output[["figures"]], "norm_cumsum_plot.pdf"),
          mito_data_cumsum_plot, base_height = 8, base_width = 12)

# This function takes the average of codon occupancy from each gene

calc_codon_occupancy <- function(data){
  # Calculate the codon frequency from all the genes taken in calculation
    codon_freq <- data %>% ungroup() %>% select(gene_name,codon_seq) %>%
        group_by(codon_seq) %>%
        summarise(codon_num = n()) %>% ungroup() %>%
        mutate(codon_freq = codon_num/sum(codon_num))

    occupancy_data <- data %>% group_by(gene_name,sample) %>% mutate(codon_num = n()) %>%
        group_by(gene_name,sample,codon_seq) %>%
       # The following is a critical step calculating codon occupancy per codon per gene. This value compares the theoretical ribosome density (which should be dependent only on codon frequency) to the measured density.
        summarise(codon_occupancy_pergene = sum(RPM_normgene)/(n()/mean(codon_num)),codon_count = n()) %>%
        group_by(sample,codon_seq) %>%
        summarise(codon_occupancy = mean(codon_occupancy_pergene,na.rm=TRUE))  %>%
        left_join(.,codon_freq,by="codon_seq") %>%
        left_join(.,mito_aminoacid_codon,by="codon_seq")
    return(occupancy_data)
}

summarise_codonoccupancy <- function(data){
    data_summary <- data %>% separate(sample,c("condition","replicate")) %>%
        group_by(condition,codon_seq,aminoacid) %>%
        summarise(avg_codon_occupancy = mean(codon_occupancy),
                  sd_codon_occupancy = sd(codon_occupancy),
                  se_codon_occupancy = sd_codon_occupancy/sqrt(n())) %>%
      group_by(codon_seq,aminoacid) %>%
      mutate(occupancy_ratio = avg_codon_occupancy/avg_codon_occupancy[condition == condition_base],
             occupancy_ratio_se = occupancy_ratio *
               sqrt((se_codon_occupancy/avg_codon_occupancy)^2 +
                      (se_codon_occupancy[condition ==condition_base]/avg_codon_occupancy[condition ==condition_base])^2))
    return(data_summary)
}

plot_occupancy <- function(occupancy_condition,occupancy_data,plot_condition,plot_ratio,add_Met){
  NNG <- c("AAA","AAG","TTA","TTG","CAA","CAG","GAA","GAG","TGA","TGG")
  Met <- c("ATG","ATA")
  occupancy_data <- occupancy_data %>%
      mutate(plot_type = case_when(!codon_seq %in% c(NNG,Met) ~"A", codon_seq %in% NNG ~ "B", codon_seq %in% Met ~ "C")) %>%
      unite(aa_codon,c("aminoacid","codon_seq"),sep = "-",remove = FALSE)
  plot_style <- list(geom_point(),
                     geom_hline(yintercept = 1),
                     theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
                           legend.position = "none",aspect.ratio = 0.442),
                     xlab(""),ggtitle(occupancy_condition))

  # if not plot ratio then we plot cell lines individually with the given input of plot_condition
  if(missing(plot_ratio)){
    output_plot <- occupancy_data  %>%
      ggplot(aes(x=reorder(codon_seq,avg_codon_occupancy),y = avg_codon_occupancy,col = plot_type)) +
      plot_style  + geom_text_repel(data = occupancy_data %>% filter(plot_type== "B"),aes(label=aa_codon)) +
      ylab("Codon occupancy")
    return(output_plot)
  }
  if(plot_ratio){
    output_plot <- occupancy_data  %>%
      ggplot(aes(x=reorder(codon_seq,occupancy_ratio),y = occupancy_ratio,col = plot_type)) +
      plot_style +
      geom_text_repel(data = occupancy_data %>% filter(plot_type== "B"),aes(label=aa_codon),col = "black") +
      ylab("Relative codon occupancy") +
      scale_color_manual(values=c("black","#FF0000","black"))
    if(add_Met){
      output_plot <- output_plot + scale_color_manual(values=c("black","#FF0000","blue")) +
        geom_text_repel(data = occupancy_data %>% filter(plot_type=="C"),aes(label=aa_codon),col = "black")
    }
    return(output_plot)
  }
}

mito_occupancy <- mito_data %>% calc_codon_occupancy(.)
mito_occupancy_summary <- mito_occupancy %>% summarise_codonoccupancy(.)

occupancy_ratio_summary_plot <- mito_occupancy_summary %>%
    group_by(condition) %>% nest() %>%
    mutate(occu_plot = map2(condition,data,~plot_occupancy(.x,.y,plot_ratio = T, add_Met = T) +
                            scale_y_continuous(limits = c(0,10))))
paths <- stringr::str_c(snakemake@output[["figures"]], "/", occupancy_ratio_summary_plot$condition, "_occupancy_plot.pdf")
pwalk(list(paths,occupancy_ratio_summary_plot$occu_plot), ggsave)

mito_occupancy_table <- mito_occupancy %>% select(sample,codon_seq,codon_occupancy) %>% spread(sample,codon_occupancy)
write_csv(x = mito_occupancy_table, path = snakemake@output[["mito_occupancy_table"]])
