load("deseq2_results.Rda")

library(tidyverse)
library(ggrepel)
theme_set(theme_minimal())

#function to sample from a prior distribution
sample_from_g <- function(g, R=1e5) {
  which.dist <- sample(1:length(g$pi),size = R, 
                       replace = TRUE, prob=g$pi)
  obs <- sapply(which.dist, function(d) rnorm(1, g$mean[d], g$sd[d]))
}

#draw 10,000 samples from priors for each taxonomic level
priors <- list(phylum=des_phylum,
               family=des_family,
               genus=des_genus,
               species=des_species) %>%
  map(~map(., attr, "g")) %>%
  map(~map(., sample_from_g, R=10000)) %>%
  map(~map(., ~data.frame(value=.))) %>%
  map(bind_rows, .id="Phenotype") %>%
  bind_rows(.id="tax_level") %>%
  mutate(tax_level=factor(tax_level, levels=unique(tax_level))) %>%
  mutate(Phenotype=var_trans[Phenotype] %>% factor(levels=var_trans)) %>%
  as_tibble()

#find modal density, for positioning of labels
modes <- priors %>% filter(tax_level=="genus",value>-3e-2,value<3e-2) %>%
  group_by(Phenotype) %>% nest %>% 
  mutate(dens=map(data, ~density(.$value,)), 
         modal_dens=map_dbl(dens, ~max(.$y)),
         mode=map_dbl(dens, function(d) d$x[which.max(d$y)]))

#Figure 1.
plot_priors <- priors %>%
  filter(tax_level=="genus") %>%
  ggplot(aes(x=value, color=Phenotype)) +
  geom_density() +
  geom_label_repel(data=modes, aes(x=mode, y=modal_dens, label=Phenotype),
                   box.padding = 1,force=10, ylim=c(0,150), size=2.5) +
  labs(x=expression(paste("Prior ", beta, " values")),y="Prior density") +
  xlim(c(-3e-2,3e-2)) +
  ylim(c(0,150)) +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(color=guide_legend(label.position = "left", title.hjust = 1))

svg(filename = "svg/Figure1.svg", height=4, width=6); plot_priors; dev.off()


#get MLE estimates, for discussion
all_models <- list(family=des_family, genus=des_genus, species=des_species) %>%
  map(~map(., attr, "ash")) %>%
  map(~map(., "result")) %>%
  map(bind_rows, .id="pheno") %>%
  bind_rows(.id="tax_level") %>%
  mutate(Phenotype=var_trans[pheno]) %>%
  as_tibble()

#density of MLE estimates, shown to be wider than priors, indicating strong shrinkage
ggplot(all_models, aes(x=betahat, color=Phenotype)) +
  geom_density() +
  facet_wrap(~tax_level) +
  coord_cartesian(xlim=c(-2,2))


#showing that priors contribute significantly more information than data
mle<-all_models %>%
  group_by(tax_level, Phenotype) %>%
  summarize(data_info=1/var(betahat)) 

priors %>% 
  group_by(tax_level, Phenotype) %>%
  summarize(prior_info=1/var(value)) %>%
  inner_join(mle) %>%
  ggplot(aes(x=log(data_info), y=log(prior_info))) + 
  geom_point()