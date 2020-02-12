load("gut_load_data_7.25.19.Rda")

##################################################
# 1. set up variables, filtering parameters, and #
#    functions for the whole pipeline            #
##################################################

library(DESeq2)
library(phyloseq)
library(tidyverse)
library(furrr)
plan(multiprocess)

#OTU data for all differential abundance analyses
otu  <- as(otu_table(pseq), "matrix")
#metadata for all differential abundance analyses
metadata <- data.frame(sample_data(pseq))
#remove hba1c samples greater than 12 (n=4)
metadata$hba1c[metadata$hba1c > 12] <- NA

#load a function to perform GMPR normalization
source("https://raw.githubusercontent.com/jchen1981/GMPR/master/GMPR.R")

#filter for all taxonomic levels: keep taxa with >3 counts present in >25 samples
otu_filter <- function(otu, min_samp=25, min_ct=3) rowSums(otu > min_ct) > min_samp



#variables to check associations with in all differential abundance analyses
vars_test <- c("ldl","hdl","glu_dec","hba1c", 
               "log_crp",names(df_rna)[-1]
)  %>% set_names(.,.)
#covariates for all DA analyses
confounders <- c(list(ldl=c("by_ctr","sex","lip_med","edu_w4"),
                      hdl=c("by_ctr","sex","lip_med","edu_w4"),
                      glu_dec=c("by_ctr","sex","glu_med","edu_w4"),
                      hba1c=c("by_ctr","sex","glu_med","edu_w4"),
                      log_crp=c("by_ctr","sex","glu_med","lip_med","crp_med","edu_w4")),
                 map(names(df_rna[,-1]) %>% set_names(.,.), 
                     ~c("by_ctr","sex","edu_w4",grep("^w5", names(ah),val=T))))[vars_test]

#setting up the pipeline                    
deseq2_pipeline <- function(exposure, 
                            covariates, 
                            otu_tab,
                            norm_factors, #user-supplied normalization factors
                            metadata,
                            standardize=is.numeric(metadata[[exposure]]),
                            svalue=TRUE,
                            ci=TRUE, 
                            ci_level=0.89,
                            keep=NULL) {
  
  #function to add user-supplied normalization factors to DESeq2
  set_size_factors <- function(dds,sf) {sizeFactors(dds) <- sf; dds}
  
  #drop cases with missing values
  if(is.null(keep))  keep <- complete.cases(metadata[,c(exposure, covariates)]) & 
      !is.na(norm_factors)
  
  #scale and center biomarker values
  if(standardize) metadata[[exposure]] <- scale(metadata[[exposure]])
  
  #generate regression formula
  form <- as.formula(paste0("~", paste(c(covariates, exposure),
                                       collapse="+")))
  #produce standard DESeq2 results using user-supplied normalization factors
  result<-DESeqDataSetFromMatrix(countData = otu_tab[,keep],
                                 colData = metadata[keep,],
                                 design= form) %>%
    set_size_factors(norm_factors[keep]) %>%
    estimateDispersions 
  
  dispersion <- as.data.frame(result@rowRanges@elementMetadata@listData)
  attr(dispersion, "dispersionFunction") <- dispersionFunction(result)
  
  result <- result %>%
    nbinomWaldTest %>%
    results(cooksCutoff = FALSE) %>%
    as.data.frame
  
  #get ASHR posterior distribution 
  dds_ash <- ashr::ash(betahat = result$log2FoldChange, #MLE estimates are used to construct the prior
                       sebetahat = result$lfcSE, "normal",
                       pointmass=FALSE, #this option sets the prior to be a Gaussian mixture 
                       outputlevel=3) 
  result$log2FoldChange <- ashr::get_pm(dds_ash) #set point estimates to be the posterior mode
  result$lfcSE <- ashr::get_psd(dds_ash) #set SE to be the posterior SD
  attr(result, "g") <- ashr::get_fitted_g(dds_ash) #return the adaptive prior
  
  if(svalue) result$svalue <- ashr::get_svalue(dds_ash) #get s-values
  if(ci) { #get highest density intervals
    ci <- ashr::ashci(dds_ash, level = ci_level)
    result$lwr <- ci[,1]
    result$upr <- ci[,2]
  }
  
  attr(result, "dispersion") <- dispersion
  attr(result, "ash") <- dds_ash
  attr(result, "call") <- match.call()
  attr(result, "n") <- sum(keep)
  result
}

##########################################
# 2. Run the pipeline at each taxa level #
##########################################

### OTU level (not reported in paper, only in eTable 1)
otu_otu <- otu[otu_filter(otu), ]
norm_factors_otu <- GMPR(otu_otu) 
des_otu <- future_map2(.x=vars_test, 
                       .y=confounders,
                       .f=~try(deseq2_pipeline(.x,.y,
                                               otu_tab=otu_otu,
                                               norm_factors=norm_factors_otu,
                                               metadata=metadata)))
save(des_otu, file="deseq2_results.Rda")

### Species level
pseq_species <- tax_glom(pseq, "Species")
pseq_species <- prune_taxa(tax_table(pseq_species)@.Data[,"Species"]
                           !="s__", pseq_species)
spec_names <- apply(tax_table(pseq_species)@.Data[,c("Genus","Species")], 1, function(x)
  paste(gsub("[a-z]__","", x), collapse=" "))
taxa_names(pseq_species) <- spec_names

otu_species <- as(otu_table(pseq_species),"matrix") 
otu_keep_species <- otu_filter(otu_species)
norm_factors_species <- GMPR(otu_species)

des_species <-  future_map2(.x=vars_test, 
                            .y=confounders,
                            .f=deseq2_pipeline,
                            .progress = TRUE,
                            otu_tab=otu_species[otu_keep_species,],
                            norm_factors = norm_factors_species,
                            metadata=metadata)
save(des_otu, des_species, file="deseq2_results.Rda")


### Genus level
pseq_genus <- tax_glom(pseq,"Genus")
pseq_genus <- prune_taxa(tax_table(pseq_genus)@.Data[,"Genus"]
                         !="g__", pseq_genus)
#manage duplicates
dups<-duplicated(tax_table(pseq_genus)@.Data[,"Genus"]) | 
  duplicated(tax_table(pseq_genus)@.Data[,"Genus"], fromLast=TRUE)
f_dups <- ifelse(dups, tax_table(pseq_genus)@.Data[,"Family"] %>% substr(4,4) %>% paste0("_",.), "")
taxa_names(pseq_genus) <- paste0(tax_table(pseq_genus)@.Data[,"Genus"], f_dups)

otu_genus <- as(otu_table(pseq_genus),"matrix") 
otu_keep_genus <- otu_filter(otu_genus)
norm_factors_genus <- GMPR(otu_genus)

des_genus <-  future_map2(.x=vars_test, 
                          .y=confounders,
                          .f=deseq2_pipeline,
                          .progress = TRUE,
                          otu_tab=otu_genus[otu_keep_genus,],
                          norm_factors = norm_factors_genus,
                          metadata=metadata)
save(des_otu, des_species, des_genus, file="deseq2_results.Rda")

### Family level
pseq_family <- tax_glom(pseq,"Family")
pseq_family <- prune_taxa(tax_table(pseq_family)@.Data[,"Family"]
                          !="f__", pseq_family)
taxa_names(pseq_family) <- tax_table(pseq_family)@.Data[,"Family"]
otu_family <- as(otu_table(pseq_family),"matrix") 
otu_keep_family <- otu_filter(otu_family)
norm_factors_family <- GMPR(otu_family,intersect.no =  5)
des_family <- future_map2(.x=vars_test, 
                          .y=confounders,
                          .f=deseq2_pipeline,
                          .progress = TRUE,
                          otu_tab=otu_family[otu_keep_family,],
                          norm_factors=norm_factors_family,
                          metadata=metadata)
save(des_otu, des_species, des_genus, des_family, file="deseq2_results.Rda")


### Phylum level
pseq_phylum <- tax_glom(pseq,"Phylum")
pseq_phylum <- prune_taxa(tax_table(pseq_phylum)@.Data[,"Phylum"]
                          !="p__", pseq_phylum)
taxa_names(pseq_phylum) <- tax_table(pseq_phylum)@.Data[,"Phylum"]
otu_phylum <- as(otu_table(pseq_phylum),"matrix") 
otu_keep_phylum <- otu_filter(otu_phylum)
norm_factors_phylum <- GMPR(otu_phylum,intersect.no =  2)
des_phylum <- future_map2(.x=vars_test, 
                          .y=confounders,
                          .f=deseq2_pipeline,
                          .progress = TRUE,
                          otu_tab=otu_phylum[otu_keep_phylum,],
                          norm_factors=norm_factors_phylum,
                          metadata=metadata)


save(des_otu, des_species, des_genus, des_family, des_phylum,
     pseq_species, pseq_genus, pseq_family, pseq_phylum, 
     otu, otu_otu, otu_species, otu_genus, otu_family, otu_phylum,
     otu_keep_species, otu_keep_genus, otu_keep_family, otu_keep_phylum,
     metadata, confounders, vars_test,
     file="deseq2_results.Rda")