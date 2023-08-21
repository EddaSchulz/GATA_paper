#Duplicates sgRNAs targeting multiple isoforms in count table
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
path <- args[2]

setwd(wd)

norm_counts <- read.delim("./mageck_output_raw/CRISPRa_screen.count_normalized.txt")
counts <- read.delim("./mageck_output_raw/CRISPRa_screen.count.txt") 
Input_dup <- read.delim(paste0(path, "/files/Library_file_dup.txt"))

#Duplicate raw counts  
seq_counts <- inner_join(counts, Input_dup) %>%
  select(-sgRNA, -Gene)

dup_counts <- left_join(Input_dup, seq_counts) %>%
  unique() %>%
  select(-sequence) %>%
  write_delim("./mageck_output_dup/CRISPRa_screen.dup_count.txt", delim = "\t")


#Duplicate normalized counts
seq_norm_counts <- inner_join(norm_counts, Input_dup) %>%
  select(-sgRNA, -Gene)

dup_norm_counts <- left_join(Input_dup, seq_norm_counts) %>%
  unique() %>%
  select(-sequence) %>%
  write_delim("./mageck_output_dup/CRISPRa_screen.dup_count_norm.txt", delim = "\t")
