#This script plots the correlation matrix for CutNTag
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), axis.title = element_blank(),
                  panel.border = element_rect(size = 0.5, color = "black", fill = NA), legend.key.size = unit(0.5,"line"),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  legend.text = element_text(size = 6), axis.ticks = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5)))


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]

setwd(wd)

input_matrix <- "Gata_CnT_pearson.txt" 

pearson_df <- read.delim(input_matrix)

pearson_df_long <- pearson_df %>%
  rename(Sample = X) %>% 
  pivot_longer(-Sample, names_to = "Sample2", values_to = "pearson")

hclust <- pearson_df %>%
  select(-X)

data <- scale(t(hclust))
ord <- hclust(dist(data, method = "euclidean"), method = "ward.D")$order

sample_levels <- unique(pearson_df$X)[ord]

pearson_plot <- pearson_df_long %>%
  ggplot(aes(x = factor(Sample, levels = sample_levels), y = factor(Sample2, levels = sample_levels), fill = pearson)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1,1), midpoint = 0) +
  labs(fill = "Pearson correlation")


fix <- set_panel_size(pearson_plot, width = unit(3.5, "cm"), height = unit(3.5, "cm"))
grid.arrange(fix)
ggsave("FigE5b_GATA_CUTnTag_correlation.pdf", fix, dpi = 300, useDingbats = FALSE)

