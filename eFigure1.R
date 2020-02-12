load("deseq2_results.Rda")

#labeller for biomarker variables
var_trans <- c(ldl="LDL",hdl="HDL",
               glu_dec="Glucose",hba1c="HbA1c",
               log_crp="CRP",infl="Hallmark Inflammation",
               n_vs_e_cd8="Naive vs. Effector CD8+",
               cd4_cd8="CD4+ vs. CD8+",
               th1_vs_n="Th1 vs. Naive CD4+",
               treg_vs_n="Treg vs. Naive CD4+",
               klrg1_vs_n="KLRG1 high vs. Naive CD8+",
               cd57_nk="CD57+ vs. - NK",
               hm_ros="Hallmark ROS",
               hm_dna="Hallmark DNA Repair",
               rna_age="Transcriptomic Age")

#eFigure 1: pairwise correlations
des_all <- list(Species=des_species,
                Genus=des_genus,
                Family=des_family,
                Phylum=des_phylum)
des_cor <- des_all %>% 
  map_depth(2,as.data.frame) %>% 
  map_depth(2, tibble::rownames_to_column, "tax") %>%
  map(bind_rows, .id="var") %>% bind_rows(.id="tax_level") %>%
  select(tax, var, log2FoldChange) %>%
  spread(key=var, value=log2FoldChange) %>%
  select(-tax) %>%
  as.matrix() %>%
  set_colnames(var_trans[colnames(.)]) %>%
  cor()
var_cor <- ah[,vars_test] %>%
  as.matrix() %>%
  set_colnames(var_trans[colnames(.)]) %>% 
  cor(use="pairwise.complete")

#this plot shows that because of the shrinkage in the hierarchical model,
#correlations between microbiome signals do not simply reflect correlations between
#variables themselves.
inner_join(reshape2::melt(var_cor, value.name="var_cor"),
           reshape2::melt(des_cor, value.name="des_cor")) %>%
  filter(var_cor < 1) %>%
  ggplot(aes(x=var_cor, y=des_cor)) + geom_point() + geom_abline(slope=1,intercept=0)

#eFigure 1
des_cor %>% gplots::heatmap.2(trace="none",col=gplots::bluered, margins=c(13,13))
