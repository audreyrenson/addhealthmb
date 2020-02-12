load("all_estimates.Rda")
load("gut_load_data_7.25.19.Rda")

library(tidyverse)
library(ggrepel)
library(scales)
theme_set(theme_minimal())

#Retreive prevalence information
prevalence <- list(Species=otu_species,
                   Genus=otu_genus,
                   Family=otu_family,
                   Phylum=otu_phylum) %>%
  map(~rowMeans(.>0)) %>%
  map(~data.frame(tax=names(.), Prevalence=unname(.), stringsAsFactors = FALSE)) %>%
  bind_rows(.id="tax_level") %>%
  filter(tax %in% all_estimates$tax)


#pull variables with no absolute logFC at least .8, and no log10 svalue at least 2
#These are labeled "other gene expression markers" in the figure
lump <- p1_dat %>% 
  group_by(var) %>% 
  summarize(l=max(abs(lfc)),x=max(x)) %>% 
  filter(l<.8 | x<2) %>% 
  pull(var) %>% 
  as.character

#color palette
col <- rev(RColorBrewer::brewer.pal(6, "Dark2")) %>% c("black")

#Figure 3
figure3 <- p1_dat %>% 
  left_join(prevalence) %>%
  mutate(lab=ifelse(abs(lfc)>.8 & x>2,gsub("\\+|-", "", lab),NA),
         var=ifelse(var %in% lump, "Other gene expression markers", as.character(var)),
         grp=ifelse(var %in% var_trans[1:4], "Metabolic", "Immune")) %>%
  #reverse code HDL
  mutate(lfc=ifelse(var=="HDL", -lfc, lfc),
         var=recode(var, HDL="HDL*")) %>%
  arrange(var=="Other gene expression markers") %>%
  mutate(var=factor(var, levels=unique(var))) %>%
  split(.$grp) %>%
  map(~ggplot(., aes(x=lfc, y=x, color=var)) + 
        geom_point(aes(size=Prevalence, alpha=is.na(lab))) +
        geom_label_repel(aes(label=lab), alpha=0.8, force=10, size=2.5, 
                         show.legend = FALSE, min.segment.length = 0)  +
        scale_color_manual(values=col, name="Biomarker") +
        scale_size_continuous(name="Taxon prevalence", breaks=c(0,0.1,0.3,0.6)) +
        scale_alpha_manual(values=c(1, 0.25), guide="none") + 
        theme(strip.text = element_blank()) +
        labs(x="log2-Fold Change", y="-log10(s-value)"))

svg("svg/Figure3.svg", width=13, height=4.5); v$Immune | v$Metabolic; dev.off()