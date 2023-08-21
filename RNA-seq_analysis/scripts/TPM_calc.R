#Calculates TPM for TX Timecourse data
library(tidyverse)
library(Rsubread)
library(Rgb)

args <- commandArgs(trailingOnly = TRUE)
bam_dir <- args[1]
wd <- args[2]


gencode_xert <- paste0(wd, "files/GENCODE_vM25_plus_Xert.gtf")

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

setwd(bam_dir)

temp = list.files(pattern="*.bam$") #Files should be named after the following pattern: [Sex]_[Day]_[Replicate]_sorted.bam

feature_counts <- featureCounts(temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, isPairedEnd = TRUE,
                                GTF.featureType = "exon", strandSpecific = 2, allowMultiOverlap = TRUE)
counts <- data.frame(feature_counts$counts, feature_counts$annotation) %>%
  select(-Chr, -Start, -End, -Strand)
names(counts) <- gsub(x = names(counts), pattern = "\\.", replacement = "_")
names(counts) <- gsub(x = names(counts), pattern = "_sorted_bam", replacement = "")

TPM <- counts_rna %>%
  mutate(Length = Length / 1000) %>%
  mutate_at(vars(-GeneID, -Length), funs( (./Length) / (sum(./Length) / 1000000))) %>%
  select(-Length) %>% 
  inner_join(gene_names) %>% 
  select(gene, gene_id = GeneID, everything())


setwd(wd)
write_delim(TPM, "RNAseq_TPM.txt", delim = "\t")
write_delim(counts_rna, "RNAseq_counts.txt", delim = "\t")


  
