#Correcting the duplicated input file to only contain the top isoform for each gene
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
path <- args[2]

setwd(wd)

dup_counts <- read.delim("./mageck_output_dup/CRISPRa_screen.dup_count.txt") 
dup_norm_counts <- read.delim("./mageck_output_dup/CRISPRa_screen.dup_count_norm.txt")
Input_dup <- read.delim(paste0(path, "/files/Library_file_dup.txt"))

#Reduced the mle table to include only the top gene per gene 
mle_red <- read.delim("./mageck_output_dup/CRISPRa_screen.gene_summary.txt") %>%
  select(Gene, Top.beta, Top.wald.fdr) %>%
  separate(Gene, c("Transcript", "Gene"), sep = "_") %>%
  group_by(Gene) %>%
  slice_min(order_by = Top.wald.fdr, n = 1, with_ties = FALSE) %>%
  unite("Gene", c(Transcript, Gene), sep = "_")


index_red_Till <- Input_dup %>% 
  select(sequence) %>% 
  unique() %>% 
  mutate(sgRNA = sequence, gene = sequence)  %>%
  write_delim("./Library_file_reduced.txt", delim = "\t")

red_counts <- dup_counts %>% 
  filter(Gene %in% c(mle_red$Gene, "NA_negative_control")) %>%
  write_delim("./mageck_output_final/CRISPRa_screen.final_count.txt", delim = "\t")

red_norm_counts <- dup_norm_counts %>% 
  filter(Gene %in% c(mle_red$Gene, "NA_negative_control")) %>%
  write_delim("./mageck_output_final/CRISPRa_screen.final_count_norm.txt", delim = "\t")