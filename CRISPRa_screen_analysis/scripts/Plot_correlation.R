#Plots correlation between samples for CRISPRa screen
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic()+
            theme(axis.title = element_text(size = 6), 
                  panel.border = element_rect(color = "black", size = 0.5, fill = NA), 
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  strip.background = element_blank(), strip.text = element_text(size = 6)))


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
path <- args[2]

ggsave("FigE1E_R1R3_Correlation.pdf", fix, dpi = 300, useDingbats = FALSE)


cor_plot_r2r3 <- counts_long %>%
  ggplot(aes(x = R2, y = R3)) +
  facet_wrap(~factor(fraction, levels = c("Unsorted", "Top")), nrow = 1) +
  geom_point(size = 0.2, alpha = 0.1) +
  labs(x = "Norm. counts (R2)", y = "Norm. counts (R3)") +
  scale_y_continuous(breaks = c(0, 2000, 4000), limits = c(0, 4500)) +
  scale_x_continuous(breaks = c(0, 2000, 4000), limits = c(0, 4500))


fix <- set_panel_size(cor_plot_r2r3, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE1E_R2R3_Correlation.pdf", fix, dpi = 300, useDingbats = FALSE)
