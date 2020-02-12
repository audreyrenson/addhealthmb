######################
# 1. Gene expression #
######################
library(Biobase)
library(tidyverse)
library(magrittr)
library(phyloseq)

### 1a. Load data

#load gene sets downloaded from amiGO2 search
amigo_gene_list <- read.table("amigo_gene_list.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)

#load Add Health gene expression data
rna <- readRDS("../../../RNA/share/zurich_subjects/dt.rds")

#load gene sets downloaded from Molecular Signature Database 
gene_sets <- c(n_vs_e_cd8="GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
               cd4_cd8="GSE8835_CD4_VS_CD8_TCELL_UP",
               th1_vs_n="GSE14308_TH1_VS_NAIVE_CD4_TCELL_UP",
               treg_vs_n="GSE20366_TREG_VS_NAIVE_CD4_TCELL_UP",
               klrg1_vs_n="GSE10239_NAIVE_VS_KLRG1HIGH_EFF_CD8_TCELL_UP",
               cd57_nk="GSE23695_CD57_POS_VS_NEG_NK_CELL_UP",
               infl="HALLMARK_INFLAMMATORY_RESPONSE",
               hm_ros="HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
               hm_dna="HALLMARK_DNA_REPAIR") %>%
  map(~read.table(paste0("gene_sets/",.,".txt"), 
                  stringsAsFactors = FALSE,
                  header=FALSE)$V1) 

### 1b. Filtering and normalization

#removing these red blood cell traits because they make up a huge majority and don't vary
rna_norb <- rna[!featureNames(rna) %in% c("HBA2","HBA1","HBB"),] 

#perform centered log-ratio normalization on gene expression data, to be used in most subsequent analyses
rna_clr <- rna_norb %>% `exprs<-`(as(t(compositions::clr(t(exprs(rna_norb)))),"matrix"))

#create a separate object containing only genes in amigo_gene_list
rna_infl <- rna_norb[featureNames(rna_norb) %in% amigo_gene_list$Gene,]


### 1c. Create Molecular Signature Database gene set PCs

#calculate PCA for each of the gene sets
rna_pca <- map(gene_sets, ~prcomp(exprs(rna_clr[featureNames(rna_clr) %in% .,])))
#obtain percent variance explained by each PC for each gene set
rna_pca %>% map(function(x) ((x$sdev^2)/sum((x$sdev^2)))[1])
#create a data frame containing the participant ID and the first PC of each gene set
df_rna <- map(rna_pca, function(x) unname(x$rotation[,1])) %>% 
{c(list(aid=rna_clr[["AID"]]), .)} %>%
  bind_cols()

### 1d. Transcriptomic age

#scale and center the gene expression data, still containing the red blood traits "HBA2","HBA1","HBB"
rna_std <- rna
exprs(rna_std) <- scale(exprs(rna_std)) 

#download the regression model coefficients for transcriptomic age
library(httr)
url1<-"https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9570/MediaObjects/41467_2015_BFncomms9570_MOESM440_ESM.xlsx"
GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))
tr_xls <- readxl::read_excel(tf, sheet = 2, skip = 3)
#filter keeping genes present in the Add Health sample
tr_B <- tr_xls[tr_xls$GeneID %in% featureNames(rna_std),] 
tr_exprs <- exprs(rna_std)[tr_B$GeneID,]

#calculate transcriptomic age by multiplying the regression coefficients with
#scaled RNA-seq values in Add Health
rna_age <- c(t(tr_exprs) %*% tr_B$`PREDICTOR-for-NEW-COHORTS` )

#retrieve chronological ages of participants to standardize the distribution of RNA age
w5_ages <- 116-Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave1/allwave1.xpt")$h1gi1y
w5_ages <- w5_ages[w5_ages>25]
rna_age <- mean(w5_ages) + (rna_age - mean(rna_age)) * (sd(w5_ages)/sd(rna_age))

#add transcriptomic age to the data frame of gene expression variables
df_rna <- df_rna %>% mutate(rna_age=rna_age)


######################
# 2. Microbiome      #
######################

### 2a. OTU table, basically in good form ##
otu <- "fecal.041719" %>%
  read.table(header = TRUE, row.names = 1, check.names = FALSE) %>%
  purrr::discard(~all(is.na(.x))) %>%
  as.matrix()

### 2b. Taxonomy table, needs some processing ##
# first get the taxon and ASV sequence columns
tax <- read.table("fecal-tables-3-04-19/taxonomy.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE)$Taxon
tax_seq <- read.table("fecal-tables-3-04-19/taxonomy.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE)$Feature.ID

#getting ready to split the taxon string into columns,
#first adding a blank of the form [kpcofgs]__ for taxa without one of those
for(i in c("k__","p__","c__","o__","f__","g__","s__")) {
  tax[!grepl(i, tax)] <- paste0(tax[!grepl(i, tax)],"; ", i) 
}
#execute the split
tax <- do.call(rbind, strsplit(tax, "; ")) %>%
  set_rownames(tax_seq) %>%
  set_colnames(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

### 2c. Phylogenetic tree 
ph_tree <- ape::read.tree("fecal-tables-3-04-19/root-tree.nwk")

### 2d. Assemble PHYLOSEQ object
pseq <- phyloseq(otu_table(as(otu,"matrix"), taxa_are_rows = TRUE), 
                 tax_table(tax),
                 phy_tree(ph_tree))

### 2e. Do some further processing of the taxonomy table
#drop empty clade classifications
tt <- tax_table(pseq)@.Data
tt[grep("__$",tt)] <- NA

#create a meaningful list of taxa names
tn <- apply(tt, 1, function(row) rev(na.omit(row))[1])
for(i in unique(unname(tn))) tn[tn==i] <- paste0(tn[tn==i], seq_len(length(tn[tn==i])))
taxa_names(pseq) <- unname(tn[taxa_names(pseq)])

### 2f. Filter read depth
#check compositional PCA location according to read depth
otu_clr <- as(compositions::clr(t(otu)), "matrix")

atch_pcoa <- ordinate(pseq, method="PCoA", distance=dist(otu_clr))

atch_pcoa_tbl <- atch_pcoa$vectors[,1:2] %>%
  as_tibble() %>%
  mutate(n_reads=cut(colSums(otu),breaks=c(0,500,1000,5000,1e4, Inf)))
ggplot(atch_pcoa_tbl, aes(x=Axis.1, y=Axis.2, color=n_reads)) +
  geom_point()

#remove samples with less than 1,000 reads
pseq <- prune_samples(colSums(otu) > 999, pseq)


################################
# 3. Add Health Interview Data #
################################

sch <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/inschool/inschool.xpt")
wave1 <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave1/allwave1.xpt")
wave4 <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/wave4.xpt")
w4meds <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/w4meds.xpt")
w4vars <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/w4vars.xpt")
wave5 <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave5/wave5.xpt")
crp <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/crp_ebv.xpt")
glu <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/glu_a1c.xpt")
lip <- Hmisc::sasxport.get("/ifs/sec/cpc/addhealth/addhealthdata/wave4/lipids.xpt")
data_dict <- c("aid"="unique_id_addhealth",
               "X.SampleID"="unique_id_microbiome",
               "shannon.diversity"="shannon",
               "faith_pd"="faith_pd",
               "evenness.pielou_e" = "evenness.pielou_e",
               "pc49c.1"="allergy_w1", 
               "pc49d.1"="asthma_w1",
               "h4id9f"="allergy_w4",
               "h4id5f"="asthma_w4_everdiag",
               "crp"="crp",
               "log_crp"="log_crp",
               "crp.flag"="crp_lowquality") #high between-duplicate differences and extreme average values

# code medication variables
w5_med_vars <- grep("(?=.*id)(?=.*m$)", names(wave5), value=TRUE, perl = TRUE)
names(w5_med_vars) <- sapply(wave5[,w5_med_vars], attr, "label")

w5_meds <- wave5[,c("aid", w5_med_vars)] %>% 
  set_names(c("aid", names(w5_med_vars))) %>%
  filter(aid %in% sample_names(pseq)) %>%
  mutate(w5_lip_med =  as.numeric(
    `S5Q6BM MED 4 HIGH CHOLEST/LIPID/TRIG-W5`==1 |
      `S5Q6CM MED FOR HIGH BLOOD PRESSURE-W5` == 1 |
      `S5Q6NM MED FOR STROKE-W5` == 1 |
      `S5Q6OM MED FOR HEART FAILURE-W5` == 1 |
      `S5Q6PM MED FOR AFIB-W5` == 1 |
      `S5Q6QM MED FOR AORTIC ANEURYSM-W5` == 1 ),
    w5_psych_med = as.numeric(
      `S5Q6GM MED FOR DEPRESSION-W5`==1 |
        `S5Q6HM MED FOR PST TRAU STR DIS/PTSD-W5` == 1 |
        `S5Q6IM MED FOR PANIC/ANXIETY DISORDER-W5` == 1 ),
    w5_asthma_med = as.numeric(`S5Q6FM MED 4 ASTHMA/CHRN BRONC/EMPHY-W5`==1),
    w5_glu_med = as.numeric(`S5Q6DM MED FOR WITH DIABETES-W5`==1)) %>%
  select(aid, starts_with("w5_"))
#merge data
wave1 <- data.frame(aid=as.character(wave1$aid),
                    birthyear=as.numeric(wave1$h1gi1y),
                    sex=as.numeric(wave1$bio.sex),
                    stringsAsFactors = FALSE)
wave4 <- data.frame(aid=as.character(wave4$aid),
                    edu_w4=case_when(wave4$h4ed2 %in% 1:3 ~ "le_hs",
                                     wave4$h4ed2 %in% 4:6 ~ "some_college",
                                     wave4$h4ed2 %in% 7:8 ~ "bach",
                                     wave4$h4ed2 %in% 9:13 ~ "ge_masters"),
                    stringsAsFactors = FALSE)
crp <- data.frame(aid=as.character(crp$aid), crp=as.numeric(crp$crp), 
                  crp_med = as.numeric(crp$crp.med8),
                  crp.flag=as.numeric(crp$crp.flag), stringsAsFactors = FALSE)
glu <- data.frame(aid=as.character(glu$aid), glucose=as.numeric(glu$glucose),
                  glu_dec=as.numeric(cut(glu$glucose, quantile(glu$glucose, probs=seq(0,1,0.1), na.rm=TRUE))),
                  hba1c=as.numeric(glu$hba1c), glu_med=glu$c.joint, 
                  stringsAsFactors = FALSE)
lip <- data.frame(aid=as.character(lip$aid), tg=lip$tg, tc=lip$tc, hdl=lip$hdl, ldl=lip$ldl,
                  non.hdl=lip$non.hdl, tc.hdl=lip$tc.hdl, lip_med=lip$c.joint2,
                  stringsAsFactors = FALSE)
ah <- wave1 %>%
  left_join(wave1) %>%
  left_join(wave4) %>%
  left_join(crp) %>%
  left_join(glu) %>%
  left_join(lip)  %>%
  left_join(w5_meds) %>%
  filter(aid %in% sample_names(pseq))


#### ---- RECODING ---- #####
ah$crp[ah$crp>997] <- NA # drop implausible values of CRP
ah$log_crp <- log(ah$crp) # take the log
ah$glucose[ah$glucose==999] <- NA #recode missing code to NA
ah$hba1c[ah$hba1c==99] <- NA#recode missing code to NA
for(i in c("tg","tc","hdl","ldl","non.hdl","tc.hdl")) ah[[i]][ah[[i]]==99] <- NA #recode missing code to NA
ah$by_ctr <- ah$birthyear - median(ah$birthyear) #centering birthyear at the median

########################
# 4. Assemble          #
########################

#add gene expression data
ah %<>% left_join(df_rna) 

rownames(ah) <- ah$aid

sample_data(pseq) <- ah[sample_names(pseq),]
sampleNames(rna_infl) <- rna_infl[["AID"]]
sampleNames(rna) <- sampleNames(rna_clr) <- rna[["AID"]]
rna_infl  <- rna_infl[, rna_infl[["AID"]] %in% sample_names(pseq)]

#create a separate phyloseq object containing only samples with both microbiome and gene expression data
pseq_subset <- prune_samples(sample_names(pseq) %in% rna_infl[["AID"]],
                             pseq)

#########################
# 5. Save and cleanup   #
#########################
save(rna, rna_infl, rna_clr, 
     gene_sets, df_rna, 
     pseq, pseq_subset, otu, otu_clr,
     ah, data_dict, 
     file = "gut_load_data_7.25.19.Rda")

rm(list=ls())