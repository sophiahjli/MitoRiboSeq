#!/usr/bin/env Rscript --vanilla

library(tidyverse, quietly=TRUE)
library(RColorBrewer, quietly=TRUE)
library(pheatmap, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(Rfast, quietly=TRUE)
library(argparser, quietly=TRUE)

p <- arg_parser("Create heatmaps from codon count data")
p <- add_argument(p, "codon_counts", help="Codon count data table")

# Add an optional arguments
p <- add_argument(p, "--heatmap_pdf",
                  help="Heatmap PDF filename",
                  default="Heatmap_codoncount_ordercodonfreq.pdf")
p <- add_argument(p, "--heatmap_scaled_pdf",
                  help="Scaled heatmap PDF filename",
                  default="Heatmap_codoncount_ordercodonfreq_scaled.pdf")
p <- add_argument(p, "--heatmap_scaled_compressed_pdf",
                  help="Scaled and compressed PDF filename",
                  default="Heatmap_codoncount_ordercodonfreq_scaled_compressed.pdf")

# Parse arguments (interactive, snakemake, or command line)
if (exists("snakemake")) {
  # Arguments via Snakemake
  argv <- parse_args(p, c(
    snakemake@input[["codon_counts"]],
    "--heatmap_pdf", snakemake@output[["heatmap_pdf"]],
    "--heatmap_scaled_pdf", snakemake@output[["heatmap_scaled_pdf"]],
    "--heatmap_scaled_compressed_pdf", snakemake@output[["heatmap_scaled_compressed_pdf"]]
  ))
} else if (interactive()) {
  # Arguments supplied inline (for debug/testing when running interactively)
  print("Running interactively...")
  input_file <- "results/codon_count/All_codoncount_table.txt"
  argv <- parse_args(p, c(input_file))
  print(argv)
} else {
  # Arguments from command line
  argv <- parse_args(p)  
}

# Read input data
mito_data_merge <- read_tsv(argv$codon_counts)

make_data_matrix <- function(data,summarize_var,column_to_make){
  col <- enquo(column_to_make)
  sum_var <- enquo(summarize_var)
  data_t <- data %>% group_by(sample,!!col) %>% summarise_at(vars(!!sum_var),funs(sum))
  matrix_t <- data_t %>% pivot_wider(names_from = !!col, values_from = !!sum_var) %>% as.data.frame()
  rownames(matrix_t) <- matrix_t$sample
  matrix_t <- matrix_t[,-1]
  return(matrix_t)
}

# Add codon annotation for start and stop codons 
mito_data_merge <- mito_data_merge %>% mutate(codon_seq = case_when(
  codon_seq %in% c("AUG","AUA") ~ paste0(codon_seq,"_Met"),
  codon_seq %in% c("UAA","UAG","AGA","AGG") ~ paste0(codon_seq,"_Stop"),
  TRUE ~ codon_seq ))

# Calculate total codon number and calculate frequency
codon_frequency <- mito_data_merge %>% mutate(ones = 1)  %>% make_data_matrix(.,ones,codon_seq)
codon_frequency <- codon_frequency/rowSums(codon_frequency) * 100
codon_frequency_order <- order(codon_frequency[1,],decreasing = T)

gene_length <- mito_data_merge %>% mutate(ones = 1) %>% make_data_matrix(.,ones,gene_id)
gene_frac <- gene_length/rowSums(gene_length) * 100
gene_frac_order <- order(gene_frac[1,],decreasing = T)

# Sum up all the raw counts first for each gene
mito_data_count_bycodon <- mito_data_merge %>% make_data_matrix(.,codon_count_sum,codon_seq) + 1
mito_data_merge_sample_name <- rownames(mito_data_count_bycodon)

# Then normalize by sample
mito_data_count_bycodon_norm <- (mito_data_count_bycodon)/rowSums(mito_data_count_bycodon)*10^2


## Make a heatmap of codon occupancy that is ordered by codon frequency
mito_data_count_bycodon_norm_order <- mito_data_count_bycodon_norm[,codon_frequency_order]

mito_data_count_bycodon_norm_order_log10 <- log10(mito_data_count_bycodon_norm_order)
mito_data_count_bycodon_norm_order_log2 <- log2(mito_data_count_bycodon_norm_order)


plot_heatmap_fraccodoncount_bycodonfrequency <- function(data_mat,...){
  pheatmap(t(data_mat),cluster_rows = F,legend=T,cellwidth = 10,cellheight = 10,...)}

codon_count_heatmap <- mito_data_count_bycodon_norm_order %>% plot_heatmap_fraccodoncount_bycodonfrequency(.,breaks = seq(0,8,length.out = 100),main = "Percentage of counts ordered by codon frequency")

mito_data_count_bycodon_norm_order_log2_rowCentered <- mito_data_count_bycodon_norm_order_log2 %>% scale(scale=F)
codon_count_scaled_heatmap <- mito_data_count_bycodon_norm_order_log2_rowCentered %>% plot_heatmap_fraccodoncount_bycodonfrequency(.,)

# We want to put a limit on the relative values to highlight the data by compressing the scale
abs_value <- 1.6
compressed_matrix <- mito_data_count_bycodon_norm_order_log2_rowCentered
compressed_matrix[compressed_matrix > abs_value] <- abs_value
compressed_matrix[compressed_matrix < -abs_value] <- -abs_value

max_inrange <- Rfast::nth(unique(as.vector(compressed_matrix)), 2, descending = T)
min_inrange <-Rfast::nth(unique(as.vector(compressed_matrix)), 2, descending = F)

range <- max(abs(c(max_inrange,min_inrange)))

paletteLength <- 100

myColor <- c("#4559b4",colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(paletteLength-1),"#d76227")
myBreaks <- c(-abs_value,seq(-range, range, length.out = paletteLength),abs_value)

codon_count_scaled_heatmap_compressed <- compressed_matrix %>%
  plot_heatmap_fraccodoncount_bycodonfrequency(.,color = myColor,breaks = myBreaks)


save_plot(
  filename = argv$heatmap_pdf,
  plot = codon_count_heatmap,
  base_height = 14,
  base_width = 14
)
save_plot(
  filename = argv$heatmap_scaled_pdf,
  plot = codon_count_scaled_heatmap,
  base_height = 14,
  base_width = 14
)
save_plot(
  filename = argv$heatmap_scaled_compressed_pdf,
  plot = codon_count_scaled_heatmap_compressed,
  base_height = 14,
  base_width = 14
)
