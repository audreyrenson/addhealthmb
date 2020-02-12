library(tidyverse)
amigo_gene_list <- read.table("amigo_gene_list.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE) %>% as_tibble()
load("gut_load_data_7.25.19.Rda")

rna_ah <- rna[, rna[["AID"]] %in% sample_names(pseq_subset)]

eTable1 <- amigo_gene_list %>%
  filter(!duplicated(Gene)) %>%
  mutate(number_positive = rowSums(exprs(rna_ah)!=0)[Gene],
         number_positive = ifelse(is.na(number_positive), 0, number_positive),
         detected_in_sample=number_positive != 0) %>%
  select(Gene, Set, Subset, Direction, number_positive)

write.csv(eTable1, file="j_gerontology_microbiome/eTable1.csv", quote = FALSE, row.names = FALSE)
