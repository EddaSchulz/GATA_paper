#Prepares peak files for use with AME
library(Rsubread)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
input_dir <- args[2]


input_dir <- 
peak_dir <- "/project/agsgpa/Till/CutnTag/Liat_Gata/ame/peaks/"
bam_dir <- "/project/agsgpa/Till/CutnTag/Liat_Gata/output_final/bam/"
output_dir <- "/project/agsgpa/Till/CutnTag/Liat_Gata/ame/input_beds/"

peak_files <- list.files(paste0(input_dir, "merged_peaks/"), pattern = ".*Gata[0-9].*.bed",
                         full.names = TRUE)
chr_vec <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")

sample_name <- function(bed_path) { str_extract(bed_path, "[A-Z0-9_a-z]*(?=.bed)")
}
center_peaks <- function(bed_path) {
  col_names <- c("chr", "start", "end")
  peak_file <- read.delim(bed_path,
                          header = FALSE,  col.names = col_names) %>% 
    mutate(center = round((start + end)/2)) %>% 
    mutate(start = center - 250, end = center + 250)  %>% 
    transmute(GeneID = paste0(chr, ":", start, "-", end),chr ,start, end, strand = ".") %>% 
    filter(chr %in% chr_vec) 
}
count_reads <- function(peak_file, sample_id) {
  setwd(paste0(input_dir, "final_bam/"))
  bam_list = list.files(pattern = paste0(sample_id, ".*.bam$"), full.names = TRUE)
  feature_counts <- featureCounts(bam_list, nthreads = 3,
                                      annot.ext = peak_file, isPairedEnd = TRUE)


}
count_rpm <- function(feature_counts) {
  total_reads <- feature_counts$stat %>% 
    rename(r1 = 2, r2 = 3) %>% 
    summarize(r1 = sum(r1), r2 = sum(r2)) %>% 
    pivot_longer(everything(), names_to = "rep", values_to = "total_reads")
  rpm_df <- data.frame(feature_counts$counts) %>% 
    rownames_to_column("peaks") %>% 
    dplyr::rename(r1 = 2, r2 = 3) %>% 
    pivot_longer(c(-1), names_to = "rep", values_to = "reads") %>% 
    left_join(total_reads) %>% 
    mutate(rpm = reads * 1000000 / total_reads) %>% 
    select(-reads, -total_reads) %>% 
    pivot_wider(names_from = rep, values_from = rpm) %>% 
    mutate(rpm = (r1 + r2)/ 2) %>% 
    select(-r1, -r2) 

}
write_output_files <- function(rpm_df, sample_id) {
  output_file <- rpm_df %>% 
    arrange(desc(rpm)) %>% 
    separate(peaks, c("chr", "location"), sep = ":") %>% 
    separate(location, c("start", "end"), sep = "-")
  
  write_delim(output_file, paste0(sample_id, "_ame_input.bed"), delim = "\t", col_names = FALSE)
}
  


sample_list <- sample_name(peak_files)
peaks_centered <- lapply(peak_files, center_peaks)
feature_counts_list <- mapply(count_reads, peaks_centered, sample_list, SIMPLIFY = FALSE)
rpm_df_list <- lapply(feature_counts_list, count_rpm)

setwd(wd)

mapply(write_output_files, rpm_df_list, sample_list, SIMPLIFY = FALSE)
