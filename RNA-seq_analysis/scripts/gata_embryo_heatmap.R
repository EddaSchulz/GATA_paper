library(EnvStats)
library(Rgb)
library(tidyverse)
library(egg)
library(gridExtra)
library(SingleCellExperiment)
library(scater)


args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
wd <- args[2]


theme_set(theme_classic() +
            theme(axis.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6, family = "Arial"),
                  axis.ticks = element_blank(), 
                  plot.title = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(size = 6, family = "Arial"), axis.text.x = element_text(angle = 90, vjust = 0.5)))
setwd(wd)

TPM <- read.delim(paste0(path, "/files/Deng_rnaseq_TPM.txt")) %>% 
  rename(X16.cell = "16-cell", X8.cell = "8-cell", X4.cell = "4-cell", X2.cell = "2-cell")

heatmap_deng <- TPM %>%
  filter(gene %in% c("Gata1", "Gata2", "Gata3", "Gata4", "Gata5", "Gata6")) %>% 
  pivot_longer(-1, names_to = "stage", values_to = "TPM") %>% 
  filter(stage %in% c("Zygote", "2-cell", "4-cell", "8-cell", "16-cell")) %>% 
  ggplot(aes(x = factor(stage, levels = c("Zygote", "2-cell", "4-cell", "8-cell", "16-cell")), 
                        y = factor(gene, levels = c("Gata6", "Gata5", "Gata4", "Gata3", "Gata2", "Gata1"))
             , fill = TPM)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#B23E8F") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) 


fix <- set_panel_size(heatmap_deng, width = unit(1.578947, "cm"), height = unit(1.894737, "cm"))
print(grid.arrange(fix))

ggsave("Fig5a_preimplant_GATA_heatmap.pdf", useDingbats=FALSE, onefile = FALSE)



#Plots argelaguet 2019 data from sce object
sce <- readRDS(paste0(path, "/files/SingleCellExperiment.rds"))


meta.col = colData(sce)
lin = meta.col$lineage
stage = meta.col$stage
sample = meta.col$sample
qc = meta.col$pass_rnaQC

meta.row = rowData(sce)
gene = meta.row$symbol

# Extract normalized log-transformed expression data
data.matrix = exprs(sce)

anno.cell = data.frame(id = colnames(data.matrix), stage, lin, sample, qc)
anno.gene = data.frame(ens = rownames(data.matrix), gene)

temp = data.frame(data.matrix,ens=rownames(data.matrix))

data.long = gather(temp,'id','counts',-ens) %>%
  left_join(anno.cell) %>% inner_join(anno.gene)

data.long.gata = data.long %>% 
  filter(gene %in% c('Gata1','Gata2','Gata3','Gata4','Gata5','Gata6'), qc==T)

lin2 = data.long.gata$lin
lin2[data.long.gata$lin=='Visceral_endoderm'] = 'Primitive_endoderm'

data.long.gata = data.long.gata %>% mutate(lin2=lin2)

data.mean = data.long.gata %>% group_by(lin,gene,stage,lin2) %>% summarize(mean.counts=mean(2^counts-1))

argelaguet_plot = data.mean %>% 
  filter((stage %in% c('E4.5','E5.5','E6.5')),lin2 %in% c('Epiblast','Primitive_endoderm')) %>%
  ggplot(aes(x = stage, y = gene, fill = mean.counts)) +
  facet_wrap(~lin2) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#B23E8F") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) 


fix <- set_panel_size(argelaguet_plot, width = unit(0.9473685, "cm"), height = unit(1.894737, "cm"))
print(grid.arrange(fix))

ggsave("Fig5b_argelaguet_GATA_heatmap.pdf", useDingbats=FALSE, onefile = FALSE)


