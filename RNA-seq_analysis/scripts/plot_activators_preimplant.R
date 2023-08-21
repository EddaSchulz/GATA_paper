#Plots heatmap for activators in Deng 2014 preimplantation RNAseq
library(tidyverse)
library(egg)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
wd <- args[2]

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6),
                  axis.title = element_blank(), panel.border = element_rect(color = "black", size = 0.5, fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  legend.text = element_text(size = 6), strip.background = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5,  color = "black"), axis.ticks = element_blank(),
                  strip.text = element_text(size = 6), legend.key.size = unit(0.5,"line")))

tpm_table <- read.delim(paste0(path, "/files/embryo_rnaseq_TPM.txt"))

setwd(wd)

tpm_transposed <- tpm_table %>% 
  pivot_longer(c(-1), names_to = "sample", values_to = "tpm") %>% 
  separate(sample, c("trash1", "trash2", "sample", "rep"), sep = "_") %>% 
  group_by(gene, sample) %>% 
  summarize(tpm = mean(tpm))

target_genes <- c("Gata1", "Xist", "Esx1", "Jpx", "Cdx4", "Nup62cl", "Tslrn1", "Ripply1", "Rlim", 
                  "Nr0b1", "Ftx", "Mid2", "Bcorl1", "Nexmif", "Map3k15", "Bcor", "Gm9785", "Eras", "Stag2")

stages <- c("Zygote", "2.cell", "4.cell", "8.cell", "16.cell", "E3.5.ICM", "E4.0.ICM", "E3.5.TE")

stage_labels <- c("Zygote", "2-cell", "4-cell", "8-cell", "16-cell", "E3.5-ICM", "E4.0-ICM", "E3.5-TE")

tpm_filt <- tpm_transposed %>% 
  filter(gene %in% target_genes)

heatmap_act <- tpm_filt %>%
  filter(sample %in% stages) %>% 
  ggplot(aes(x = factor(sample, levels = stages), y = factor(gene, levels = rev(target_genes)), fill = log2(tpm + 1))) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "white", high = "red") +
  scale_x_discrete(expand = c(0,0), labels = stage_labels) +
  scale_y_discrete(expand = c(0,0)) 


fix <- set_panel_size(heatmap_act, width = unit(1.6, "cm"), height = unit(3.8, "cm"))
print(grid.arrange(fix))
ggsave("FigE2I_Activators_heatmap_embryo.pdf", useDingbats=FALSE, onefile = FALSE)
