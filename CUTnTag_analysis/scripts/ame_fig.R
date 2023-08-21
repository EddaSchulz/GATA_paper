library(tidyverse)
library(egg)
library(gridExtra)
library(ggrepel)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 6),
                  plot.title = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(size = 6)))

args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]


setwd(paste0(wd, "ame_data/")) 

dir_list <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]

read_ame <- function(a) {
  gata = str_extract(a, "Gata[0-9]") 
  type = ifelse(str_detect(a, "TS"), "TS", "XEN")
  
  ame <- read.delim(paste0(a, "/ame.tsv")) %>% 
    select(motif_ID, motif_alt_ID, E.value) %>% 
    mutate(motif_fam = ifelse(str_detect(motif_alt_ID, "GATA"), "GATA-motif", "Non-GATA-motif"), 
           factor = gata, line = type) %>% 
    na.omit() %>% 
    filter(E.value <= 10) %>% 
    rowid_to_column("ID")
}

gata_df <- bind_rows(lapply(dir_list, read_ame))

ame_plot <- gata_df %>% 
  ggplot(aes(x = ID, y = -log10(E.value), color = motif_fam)) +
  facet_wrap(line~factor, scales = "free_x") +
  geom_point(size = 0.5) +
  geom_text_repel(data = gata_df[gata_df$ID<=3,], aes(label = motif_alt_ID), size = 2.142857) +
  labs(x = "Rank", y = "E-Value (-log10)") +
  scale_color_manual(values = c("#1D71B8", "grey"))

fix <- set_panel_size(plot_gata, width = unit(3, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE5d_GATA_ame.pdf", fix, dpi = 300, useDingbats = FALSE)

