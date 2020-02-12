#load data saved by "load_data_from_source.R
load("gut_load_data_7.25.19.Rda")

#load packages
library(tidyverse)
library(phyloseq)
theme_set(theme_minimal())

#calculate compositional PCA
atch_pcoa <- ordinate(pseq, method="PCoA", distance=dist(otu_clr))

#attach number of reads
atch_pcoa_tbl <- atch_pcoa$vectors[,1:2] %>%
  as_tibble() %>%
  mutate(n_reads=colSums(otu), 
         n_reads_cut = cut(colSums(otu),breaks=c(0,500,1000,5000,1e4, Inf),
                     labels = c("<500","500 - <1,000", "1,000 - <5,000",
                                "5,000 - <10,000","10,000+"))) %>%
  group_by(n_reads_cut) %>%
  mutate(n=n(), min=min(n_reads)) %>%
  ungroup() %>%
  arrange(min) %>%
  mutate(n_reads_cut = paste0(n_reads_cut, " (", n, ")"),
         n_reads_cut = factor(n_reads_cut, levels=unique(n_reads_cut)))

#plot
ggplot(atch_pcoa_tbl, aes(x=Axis.1, y=Axis.2, color=n_reads_cut)) +
  geom_point() +
  scale_color_discrete(name="# reads in sample (N)") +
  labs(x=paste0("Axis 1 (", round(100*atch_pcoa$values$Relative_eig[1], 1), "%)"),
       y=paste0("Axis 2 (", round(100*atch_pcoa$values$Relative_eig[2], 1), "%)"))
