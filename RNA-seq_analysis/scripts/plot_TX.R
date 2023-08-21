#Plots heatmap for RNAseq timecourse
#Plots lineplot of Gatafactors
library(tidyverse)
library(egg)
library(gridExtra)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
wd <- args[2]

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), panel.border = element_rect(color = "black", size = 0.5, fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  legend.text = element_text(size = 6), strip.background = element_blank(),
                  strip.text = element_text(size = 6), legend.key.size = unit(0.5,"line")))

TPM <- read.delim(paste0(path, "/files/RNAseq_TPM.txt"))

setwd(wd)

TPM_long <- TPM %>%
  select(-GeneID) %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
  separate(sample, c("rep", "timepoint", "sex"), sep = "_")


#Plots lineplot of Xist and Gata factors
genes <- c("Gata6", "Gata5", "Gata4", "Gata3", "Gata2", "Gata1", "Xist")


line_tpm <- TPM_long %>%
  filter(gene %in% genes) %>% 
  ggplot(aes(x = timepoint, y = tpm, color = gene)) +
  facet_wrap(~factor(sex, levels = c("XX", "XO"))) +
  stat_summary(aes(group = gene), geom = "line", fun = "mean") +
  scale_color_manual(values = c(brewer.pal(6, "Blues"), "deeppink3")) +
  scale_x_discrete(labels = c(0, 1, 2, 3, 4)) +
  labs(x = "Day", y = "TPM")


fix <- set_panel_size(line_tpm, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE5E_Gata_lineplot.pdf", fix, dpi = 300, useDingbats=FALSE)

#Plotting the hits as heatmap...
auto_conts <- c("Nanog", "Zfp42", "Esrrb", "Yy1", "Sox2", "Myc", "Klf4", "Pou5f1", "Prdm14", "Eed", "Ctcf", "Chd8",
                "Kat8", "Msl1", "Msl2", "Kansl1", "Kansl3", "Dnmt1", "Mcrs1", "Spen", "Lbr", "Hnrnpu", "Hnrnpk")
x_conts <- c("Xist", "Rlim", "Ftx", "Jpx")


mle_final <- read.delim(paste0(path, "/files/CRISPRa_screen.gene_summary.txt")) %>%
  separate(Gene, c("Transcript", "Gene"), sep = "_") %>%
  mutate(color = ifelse(Gene %in% auto_conts, "Auto_cont", 
                        ifelse(Gene %in% x_conts, "X_cont", "X-linked")))

col_df <- mle_final %>% 
  filter(Top.wald.fdr <= 0.05 & Top.beta >= 0) %>%
  mutate(color = ifelse(color == "Auto_cont", "#CC0000",
                        ifelse(color == "X-linked", "#009999", "#FF9900"))) %>% 
  arrange(Top.beta) 

col_vec <- col_df$color
gene_vec <- col_df$Gene

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), strip.background = element_blank(),
                  panel.border = element_rect(size = 0.5, color = "black", fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.ticks = element_blank(), axis.title.y = element_blank(),
                  legend.key.size = unit(0.5,"line"), axis.title.x = element_text(size = 6),
                  axis.text.y = element_text(color = col_vec), strip.text = element_text(size = 6),
                  legend.text = element_text(size = 6, color = "black")))

combi_df <- col_df %>% 
  select(gene = Gene, Top.beta, Top.wald.fdr) %>% 
  inner_join(TPM_long)

heatmap_act <- combi_df %>%
  filter(Top.beta >= 0) %>% 
  ggplot(aes(x = timepoint, y = factor(gene, levels = gene_vec), fill = log2(tpm + 1))) +
  facet_wrap(~factor(sex, levels = c("XX", "XO"))) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red") +
  theme() +
  scale_x_discrete(expand = c(0,0), labels = c(0, 1, 2, 3, 4)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "Day")



fix <- set_panel_size(heatmap_act, height = unit(3.8, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("FigE2E_Xist_activator_heatmap.pdf", fix, dpi = 300, useDingbats=FALSE)