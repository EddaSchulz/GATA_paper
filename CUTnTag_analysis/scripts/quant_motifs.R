#This script calculates RPM within CUT&Tag peaks and plots motif vs no-motif
library(Rsubread)
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 6),
                  plot.title = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(size = 6)))

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
input_dir <- args[2]

#Counts reads in the different GATA CUT&Tag samples in peaks with/without motif
gata_vec <- c("Gata2", "Gata3", "Gata4", "Gata6")

fimo_function <- function(a) {
  setwd(wd)
  col_names <- c("chr", "start", "end", "trash", "trash2", "trash3", "trash4", "trash5", "color")
  
  peak_file = list.files(pattern=paste0(".*", a, ".*_rgb_motif.bed$"), full.names = TRUE)
  
  fimo_gata <- read.delim(peak_file,
                           header = FALSE,  col.names = col_names) %>% 
    mutate(gata_motif = ifelse(color == "150,150,150", FALSE, TRUE))
  
  fimo_gata_true <- fimo_gata %>% 
    filter(gata_motif == TRUE) %>% 
    transmute(GeneID = paste0(chr, ":", start, "-", end),chr ,start, end, strand = ".")
  
  fimo_gata_false <- fimo_gata %>% 
    filter(gata_motif == FALSE) %>% 
    transmute(GeneID = paste0(chr, ":", start, "-", end),chr ,start, end, strand = ".")
  
  setwd(paste0(input_dir, "final_bam/"))
  temp_gata = list.files(pattern=paste0(".*", a, ".*.bam$"), full.names = TRUE) 
  
  feature_counts_gata_true <- featureCounts(temp_gata, 
                                             annot.ext = fimo_gata_true, isPairedEnd = TRUE)
  
  feature_counts_gata_false <- featureCounts(temp_gata,
                                              annot.ext = fimo_gata_false, isPairedEnd = TRUE)
  
  total_reads_sums <- colSums(feature_counts_gata_true$stat[,2:3])
  
  total_reads_gata <- data.frame(rep = c("r1", "r2"), total_reads = total_reads_sums)
  
  
  counts_gata_true <- data.frame(feature_counts_gata_true$counts) %>% 
    mutate(gata_motif = TRUE)
  
  counts_gata_false <- data.frame(feature_counts_gata_false$counts) %>% 
    mutate(gata_motif = FALSE)
  
  counts_gata <- rbind(counts_gata_true, counts_gata_false) 
  names(counts_gata) <- gsub(x = names(counts_gata), pattern = "\\.", replacement = "_")
  names(counts_gata) <- gsub(x = names(counts_gata), pattern = "_bam", replacement = "")
  
  sum_counts_gata <- counts_gata %>% 
    rownames_to_column("peak") %>% 
    pivot_longer(c(2:3), names_to = "sample", values_to = "reads") %>% 
    separate(sample, c("line", "sex", "factor", "rep"), sep = "_") %>% 
    left_join(total_reads_gata)
  
  return(sum_counts_gata)
}

total_counts_gata <- bind_rows(lapply(gata_vec, fimo_function))

setwd(wd)

total_rpm_gata <- total_counts_gata %>% 
  mutate(rpm = reads / (total_reads / 1000000))

plot_gata <- total_rpm_gata %>% 
  select(-reads, -total_reads) %>% 
  pivot_wider(names_from = rep, values_from = rpm)  %>% 
  mutate(rpm = (r1 + r2)/ 2) %>% 
  ggplot() +
  facet_wrap(~factor) +
  geom_density(aes(x = log2(rpm), color = factor(gata_motif, levels = c(TRUE, FALSE)))) +
  labs(x = "RPM (log2)", y = "Peak density (a.u.)") +
  scale_color_manual(values = c("#1D71B8", "grey"), labels = c("GATA-motif", "No GATA-motif")) +
  scale_x_continuous(limits = c(-1, 9), breaks = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4))


fix <- set_panel_size(plot_gata, width = unit(3, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE5c_GATA_signal_enrichment.pdf", fix, dpi = 300, useDingbats = FALSE)

gata_txt <- total_counts_gata %>% 
  group_by(gata_motif, line, sex, factor) %>% 
  summarize(mean_reads = mean(reads)) %>% 
  pivot_wider(names_from = gata_motif, values_from = mean_reads) %>% 
  rename(RPM_motif = `TRUE`, RPM_no_motif = `FALSE`) %>% 
  mutate(LFC_motif = round(log2(RPM_motif/RPM_no_motif), 2))

write_delim(gata_txt, "GATA_motif_expression.txt")

peaks_xic <- total_rpm_gata %>% 
  separate(peak, c("chr", "location"), sep = ":") %>% 
  separate(location, c("start", "end"), sep = "-") %>% 
  filter(chr == "chrX" & end >= 103182701 & start <= 103955531) %>% 
  select(-reads, -total_reads) %>% 
  pivot_wider(names_from = rep, values_from = rpm) %>% 
  mutate(rpm = (r1 + r2)/ 2) %>% 
  select(-r1, -r2) %>% 
  unite("Condition", line, sex, factor) %>% 
  select(Condition, chr, start, end, rpm, gata_motif)

write_delim(peaks_xic, "GATA_XIC_peaks.txt", delim = "\t")
