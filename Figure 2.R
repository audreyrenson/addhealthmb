load("gut_load_data_7.25.19.Rda")
load("deseq2_results.Rda")

library(tidyverse)
library(ggtree)
library(ggrepel)
theme_set(theme_minimal())

###########################################
# 1. Assemble the data for manhattan plot #
###########################################

otu_tax <- tax_table(pseq)@.Data

all_estimates <- list(Phylum=des_phylum, Family=des_family, Genus=des_genus, Species=des_species) %>%
  map(map, tibble::rownames_to_column, "otu_id") %>%
  map(bind_rows, .id="pheno") %>%
  bind_rows(.id="tax_level") %>%
  select(tax_level, pheno, otu_id, log2FoldChange, svalue, pvalue, lwr, upr) %>%
  mutate(otu_id_higher = c(str_match(otu_id, "_[A-Z]$")),
         otu_id_higher = gsub("_","", otu_id_higher),
         is_ambig = !is.na(otu_id_higher),
         otu_id = gsub("_[A-Z]$", "", otu_id)) %>%
  #disambiguate genera that exist in multiple families
  mutate(families = map(otu_id, 
                        function(o) unique(otu_tax[otu_tax[,"Genus"]==o, "Family"])),
         families = map(families, ~gsub("f__","",.)),
         pattern = paste0("^", otu_id_higher),
         family = map2(pattern, families, ~grep(.x,.y, value=TRUE)[1]),
         family=unlist(ifelse(!is.na(otu_id_higher), family, NA))) %>%
  ### This is to start the process of randomly picking an OTU to
  #represent larger cateogires 
  #start by grabbing the right section of the taxa table
  mutate(phylo_group = map2(tax_level, otu_id, 
                            function(.x,.y) otu_tax[otu_tax[,.x]==.y,,drop=FALSE]))

#disambiguate where necessary (Genus)
a<-all_estimates %>% filter(!is.na(otu_id_higher))
p <- a$phylo_group
f <- a$family 
disambig_phylo <- function(p, f) p[p[,"Family"]==paste0("f__", f),,drop=FALSE]
a$phylo_group <-map2(p,f, disambig_phylo)
all_estimates[all_estimates$is_ambig, ] <- a

#disambiguate where necessary (Species)
s <- all_estimates %>% filter(tax_level == "Species") %>%
  mutate(split=strsplit(otu_id, " "),
         species=paste0("s__", map_chr(split, 2)), 
         genus = paste0("g__", map_chr(split, 1)),
         phylo_group = map2(genus, species, 
                            function(g,s) otu_tax[otu_tax[,"Genus"]==g &
                                                    otu_tax[,"Species"]==s,,drop=FALSE])) %>%
  select(-split,-species,-genus)
all_estimates[all_estimates$tax_level=="Species",] <- s

#now select the random otu (make it the same random otu for a given taxon)
set.seed(10)
all_estimates %<>%
  mutate(random_otu = map_chr(phylo_group, ~try(sample(rownames(.),1)))) %>%
  group_by(otu_id) %>%
  mutate(random_otu = dplyr::first(random_otu)) %>%
  ungroup() %>%
  select(id=random_otu, tax=otu_id, tax_level, var=pheno, log2FoldChange, svalue, lwr,upr)

all_estimates <- bind_rows(list(all_estimates) 


#group the tree according phylum, for tree colors
tr <- phy_tree(pseq)
tax <- tax_table(pseq)
tr <- groupOTU(.data=tr,
               .node=split(tr$tip.label, tax[tr$tip.label, "Phylum"]),
               group_name = "Phylum")
attr(tr, "Phylum") <- as.character(attr(tr, "Phylum"))
attr(tr, "Phylum")[attr(tr, "Phylum")=="0"] <- "unclassified"
attr(tr, "Phylum") <- factor(attr(tr, "Phylum"))
ordered_tips <- rev(na.omit(tr$tip.label[tr$edge[,2]])) #this ensures ordering of dots is the same as the tree
tr <- ape::rotateConstr(tr, ordered_tips)


#setup the data for dots
p1_dat <- all_estimates %>%
  mutate(x=-log10(svalue), 
         id=factor(id, levels=ordered_tips),
         var=factor(var_trans[var], levels=var_trans),
         sign=ifelse(log2FoldChange<0, "-", "+"),
         lab=paste(sign,gsub("[a-z]__","",tax)),
         tax_level=factor(tax_level, levels=c("Phylum","Family","Genus","Species",
                                              "OTU"))) %>%
  filter(!tax_level %in% c("OTU"))%>%
  select(id, tax,var, x, sign, lab,tax_level, lfc=log2FoldChange)

####################################
# 2. Create tree portion of plot   #
####################################

#This will be saved and pasted together with the manhattan plot in inkscape

p2 <- ggtree(tr, aes(color=Phylum), ladderize = FALSE) +
  theme_tree2() +
  #scale_color_manual(values=tr.pallete) +
  theme(legend.position = "right")

# # check that the dots line up
# p2<- facet_plot(p2, panel="dot", data=p1_dat, 
#                 geom = geom_point, mapping=aes(x=x)) + 
#   theme(strip.text = element_blank(), legend.position="right") + 
#   ggtitle(seed)

#save this
svg(filename="svg/phylo_allpheno.svg", width = 3.5,hei=7); p2; dev.off()

####################################
# 2. Create dot portion of plot   #
####################################

#lump measures together what have no log10 svalues greater than 3
lump <- p1_dat %>% 
  group_by(var) %>% 
  summarize(m=max(x)) %>% 
  filter(m<3) %>% 
  pull(var) %>% 
  as.character

#dot plot
p1b <- p1_dat %>%
  mutate(lab=ifelse(x>3,lab,NA),
         lab=gsub("\\[|\\]","", lab),
         var=ifelse(var %in% lump, "Other gene expression markers", as.character(var))) %>%
  ggplot(aes(x=id, y=x, color=var, 
             group=var, shape=tax_level)) + 
  geom_point(alpha=0.8) +
  geom_label_repel(aes(label=lab,x=id,y=x, color=var), 
                   inherit.aes = FALSE, force=20, size=2.5, 
                   show.legend = FALSE, min.segment.length = 0) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "right",
        panel.grid.major.x = element_blank()) +
  labs(x=NULL, y="-log10(S-value)") +
  scale_color_brewer(palette = "Dark2",direction = -1) +
  guides(color=guide_legend(title="Biomarker"),
         shape=guide_legend(title="Taxonomic rank"))

svg(filename = "svg/Figure2.svg", width=11, height=5); p1b; dev.off()

#save data assembly
save(all_estimates, p1_dat, file="all_estimates.Rda")
