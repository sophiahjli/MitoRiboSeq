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
  occupancy_data <- occupancy_data %>% mutate(plot_type = case_when(!codon_seq %in% c(NNG,Met) ~"A", codon_seq %in% NNG ~ "B", codon_seq %in% Met ~ "C")) %>% unite(aa_codon,c("aminoacid","codon_seq"),sep = "-",remove = FALSE)
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

