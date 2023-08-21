#Plots a volcano plot of the mageck mle results
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

setwd(wd)

volcano <- mle_final %>%
  ggplot(aes(x = Top.beta, y = -log10(Top.wald.fdr))) +
  geom_point(aes(color = color), size = 0.8) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dotted") +
  geom_text(data = mle_final[mle_final$Top.wald.fdr < 0.05,], aes(label = Gene, color = color), 
            size = 6 / 2.8) +
  xlim(-1.05, 1.05)  + 
  ylim(0, NA) +
  labs(y = "-log10 FDR Top/Unsorted", x = "Beta score Top/Unsorted", aes(color = color, text = Gene)) +
  scale_color_manual(values = c("#CC0000", "#009999", "#FF9900"))

fix <- set_panel_size(volcano, width = unit(5, "cm"), height = unit(5, "cm"))
grid.arrange(fix)
ggsave("FigE1G_CRISPRa_screen_volcano.pdf", fix, dpi = 300, useDingbats = FALSE)