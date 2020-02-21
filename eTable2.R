load("gut_load_data_7.25.19.Rda")
load("deseq2_results.Rda")
load("all_estimates.Rda")

library(tidyverse)

des_all <- list(otu=des_otu,
                Species=des_species,
                Genus=des_genus,
                Family=des_family,
                Phylum=des_phylum)
otu_all <- list(otu=otu_otu,
                Species=otu_species,
                Genus=otu_genus,
                Family=otu_family,
                Phylum=otu_phylum)

#get mean relative abundance for all taxa
relabund_all <- map(otu_all, ~sweep(., 2, colSums(.), "/"))

#calculate HPDIs for all estimates
posteriors <- des_all %>%
  map_depth(attr,.depth = 2, "ash") %>%
  map_depth(ashr::get_post_sample, .depth = 2, 1e4) 

hpdi <- map_depth(posteriors, .depth=2, apply, 2, HDInterval::hdi, credMass=.9) %>% 
  map_depth(2, t) %>%
  map_depth(2,as.data.frame) %>%
  imap(~map(.x, mutate, tax=rownames(des_all[[.y]][[1]]))) %>%
  map(bind_rows, .id="var") %>%
  bind_rows(.id="tax_level")

#assemble
eTable2 <- all_estimates %>% 
  left_join(hpdi) %>%
  mutate(text_lfc_hpid = paste0(formatC(log2FoldChange, digits=1, format="f"),
                                " [", formatC(lower, digits=1, format="f"), "; ",
                                formatC(upper, digits=1, format="f"),"]")) %>%
  mutate(prevalence = map2_chr(tax, tax_level, ~try(mean(otu_all[[.y]][.x,]>0)))) %>%
  mutate(mean_abundance = map2_chr(tax, tax_level, ~try(relabund_all[[.y]][.x,] %>% {mean(.[.>0])}))) %>%
  mutate_at(vars(prevalence,mean_abundance), as.numeric) 


write.csv(eTable2, file="eTable2.csv", quote=FALSE, row.names = FALSE)  