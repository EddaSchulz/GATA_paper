#Plots log2 foldchanges of the screen hits (including guides)
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 6)))


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
path <- args[2]

setwd(path)

auto_conts <- c("Nanog", "Zfp42", "Sox2", "Myc", "Klf4", "Esrrb", "Pou5f1", 
                "Prdm14", "Ctcf", "Yy1", "Chd8", "Myst1", "Msl1", "Msl2",
                "Kansl3", "Kansl1", "Mcrs1", "Dnmt1", "Spen", "Lbr", "Hnrnpu",
                "Hnrnpk", "Kat8")
x_conts <- c("Jpx", "Ftx", "Rlim", "Xist")


mle_final <- read.delim("./files/CRISPRa_screen.gene_summary.txt") %>%
  separate(Gene, c("Transcript", "Gene"), sep = "_") %>%
  mutate(color = ifelse(Gene %in% auto_conts, "Auto_cont", 
                        ifelse(Gene %in% x_conts, "X_cont", "X-linked")))

col_df <- mle_final %>% 
  filter(Top.wald.fdr <= 0.05 & Top.beta >= 0) %>%
  mutate(color = ifelse(color == "Auto_cont", "#CC0000",
                        ifelse(color == "X-linked", "#009999", "#FF9900"))) %>% 
  arrange(Top.beta) %>% 
  select(Gene, color)

col_vec <- col_df$color
gene_vec <- col_df$Gene

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#In this part, I am plotting the fold change of the individual guides
red_norm_counts <- read.delim("./files/CRISPRa_screen.final_counts_norm.txt")
lfc_guides <- red_norm_counts %>%
  transmute(sgRNA, Gene, R1_Top = R1_Top/R1_Unsorted, R2_Top = R2_Top/R2_Unsorted, R3_Top = R3_Top/R3_Unsorted) %>%
  filter(Gene != "NA_negative_control") %>% 
  separate(Gene, c("Transcript", "Gene"), sep = "_") %>%
  transmute(sgRNA, Gene, Transcript, LFC = log2((R1_Top + R2_Top + R3_Top) / 3)) %>%
  mutate(Gene = ifelse(Gene == "1700019B21Rik", "Tslrn1", 
                       ifelse(Gene == "C77370", "Nexmif", Gene))) %>% 
  inner_join(col_df)

theme_set(theme_classic() +
            theme(axis.title.y = element_text(size = 6, color = "black"), 
                  axis.title.x = element_blank(),
                  panel.border = element_rect(size = 0.5, color = "black", fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.text.y = element_text(color = "black"), legend.key.size = unit(0.5,"line"),
                  axis.text.x = element_text(color = rev(col_vec), angle = 90, vjust = 0.5, hjust=1),
                  legend.position = "none"))

setwd(wd)

guide_plot <- lfc_guides %>% 
  ggplot(aes(x = factor(Gene, levels = rev(gene_vec)), y = LFC, color = color)) +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", alpha = 0.2, aes(fill = color), width = 0.75) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#009999", "#FF9900")) +
  scale_fill_manual(values = c("#009999", "#FF9900"))

fix <- set_panel_size(guide_plot, width = unit(3.8, "cm"), height = unit(3, "cm"))
grid.arrange(fix)
ggsave("Fig1B_activators_lfc.pdf", fix, dpi = 300, useDingbats=FALSE)

#In this part, I am plotting the fold change of the individual guides for the repressors
rep_col_df <- mle_final %>% 
  filter(Top.wald.fdr <= 0.05 & Top.beta <= 0) %>%
  mutate(-sgRNA, color = ifelse(color == "Auto_cont", "#CC0000",
                                ifelse(color == "X-linked", "#009999", "#FF9900"))) %>% 
  arrange(Top.beta) %>%
  select(Gene, Top.beta, color)

rep_col_vec <- rep_col_df$color
rep_gene_vec <- rep_col_df$Gene

rep_lfc_guides <- red_norm_counts %>% 
  transmute(sgRNA, Gene, R1_Top = R1_Top/R1_Unsorted, R2_Top = R2_Top/R2_Unsorted, R3_Top = R3_Top/R3_Unsorted) %>%
  filter(Gene != "NA_negative_control") %>% 
  separate(Gene, c("Transcript", "Gene"), sep = "_") %>%
  transmute(sgRNA, Gene, Transcript, LFC = log2((R1_Top + R2_Top + R3_Top) / 3)) %>%
  inner_join(rep_col_df)

theme_set(theme_classic() +
            theme(axis.title.y = element_text(size = 6, color = "black"), 
                  axis.title.x = element_blank(),
                  panel.border = element_rect(size = 0.5, color = "black", fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.text.y = element_text(color = "black"), legend.key.size = unit(0.5,"line"),
                  axis.text.x = element_text(color = rep_col_vec, angle = 90, vjust = 0.5, hjust=1),
                  legend.position = "none"))

rep_guide_plot <- rep_lfc_guides %>% 
  ggplot(aes(x = factor(Gene, levels = rep_gene_vec), y = LFC, color = color)) +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot", alpha = 0.2, aes(fill = color), width = 0.75) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#009999", "#CC0000")) +
  scale_fill_manual(values = c("#009999", "#CC0000"))

fix <- set_panel_size(rep_guide_plot, width = unit(8.6, "cm"), height = unit(3, "cm"))
grid.arrange(fix)
ggsave("Fig1C_repressors_lfc.pdf", fix, dpi = 300, useDingbats=FALSE)
