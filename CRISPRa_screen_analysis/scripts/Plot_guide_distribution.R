#Plots quality control metrics for CRISPRa screen
library(tidyverse)
library(egg)
library(gridExtra)


theme_set(theme_classic() +
    theme(legend.title = element_blank(), legend.text = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
          strip.background = element_blank(), axis.title = element_text(size = 6),  
          axis.text = element_text(size = 6), strip.text = element_text(size = 6)))


args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]
path <- args[2]

setwd(path)


norm_counts <- read.delim("./files/CRISPRa_screen.count_normalized.txt")
raw_counts <- read.delim("./files/CRISPRa_screen.count.txt")

setwd(wd)

#Prepares count table for plotting
counts_long <- norm_counts %>%
  pivot_longer(-c(sgRNA, Gene), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("replicate", "fraction"))

unsorted_counts <- counts_long %>%
  filter(fraction == "Unsorted") %>%
  pivot_wider(names_from = fraction, values_from = counts)

#Calculates log2 distribution width of the guides
width_df <- norm_counts %>%
  filter(Gene != "NA_negative_control") %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)

sample_levels = rev(c("R3_Top", "R3_Unsorted", "R2_Top", "R2_Unsorted", 
                      "R1_Top", "R1_Unsorted", "PlasmidLib"))

width_plot <- width_df %>%
  filter(sample %in% sample_levels) %>% 
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point(size = 2) +
  labs(x = "Distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1.5, 3)

fix <- set_panel_size(width_plot, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE1F_Distribution_width.pdf", fix, dpi = 300, useDingbats = FALSE)


nt_df <- norm_counts %>%
  filter(Gene == "NA_negative_control") %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)

ntc_plot <- nt_df %>%
  filter(sample %in% sample_levels) %>% 
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point(size = 2) +
  labs(x = "NTC distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1.5, 3)


fix <- set_panel_size(ntc_plot, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE1F_NTC_Distribution_width.pdf", fix, dpi = 300, useDingbats = FALSE)


#Plots cumulative distributions for the individual fractions
width_df_2 <- width_df %>%
  filter(sample != "PlasmidLib") %>% 
  separate(sample, c("replicate", "fraction"), sep = "_") %>% 
  filter(fraction == "Top")

cum_plot_df <- counts_long %>%
  filter(fraction != "Unsorted") %>%
  left_join(unsorted_counts) %>%
  rename(Sorted = counts) %>%
  pivot_longer(c(Sorted, Unsorted), names_to = "curve", values_to = "counts") %>% 
  filter(fraction == "Top") %>% 
  left_join(width_df_2)

cumulative_plot <- cum_plot_df %>%
  ggplot(aes(x = log2(counts), color = curve)) +
  facet_wrap(~ replicate, ncol = 3) +
  stat_ecdf(geom = "line", size = 0.3) +
  scale_color_manual(values = c("red", "black")) +
  xlim(5, 13) +
  geom_vline(aes(xintercept = width_low), linetype = "dashed", size = 0.3) +
  geom_vline(aes(xintercept = width_high), linetype = "dashed", size = 0.3) +
  labs(x = "Normalized sgRNA counts (log2)", y = "Cumulative sgRNA frequency")

fix <- set_panel_size(cumulative_plot, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE1D_Cumulative_distribution.pdf", fix, dpi = 300, useDingbats = FALSE)

width_plasmid <- raw_counts %>%
  summarize(width_low = log2(quantile(PlasmidLib, 0.1)), width_high = log2(quantile(PlasmidLib, 0.9))) %>%
  mutate(width = width_high - width_low)

cumulative_plasmid <- raw_counts %>% 
  ggplot(aes(x = log2(PlasmidLib))) +
  stat_ecdf(geom = "line", size = 0.3, color = "black") +
  xlim(5, 11) +
  geom_vline(aes(xintercept = width_plasmid$width_low), linetype = "dashed", size = 0.3) +
  geom_vline(aes(xintercept = width_plasmid$width_high), linetype = "dashed", size = 0.3) +
  labs(x = "Raw sgRNA counts (log2)", y = "Cumulative sgRNA frequency")

fix <- set_panel_size(cumulative_plasmid, width = unit(2, "cm"), height = unit(2, "cm"))
grid.arrange(fix)
ggsave("FigE1C_Cumulative_distribution_plasmid.pdf", fix, dpi = 300, useDingbats = FALSE)