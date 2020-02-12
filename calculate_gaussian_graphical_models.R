load("gut_load_data_7.25.19.Rda")
load("deseq2_results.Rda") #just loading this for the variable 'confounders'

library(tidyverse)
library(huge)

#functions------------------
edge.weights <- function(community, network, weight.within=100, weight.between=1) {
  bridges=igraph::crossing(communities = community, graph=network)
  weights=ifelse(test=bridges, yes=weight.between, no=weight.within)
}
get_degree_dists <- function(huge_obj) {
  list_degree <- lapply(huge_obj$path, function(x) 
    data.frame(degree=igraph::degree(graph_from_adjacency_matrix(x))))
  names(list_degree) <- round(huge_obj$lambda,4)
  list_degree
}
plot_degree_dist <- function(huge_obj) {
  get_degree_dists(huge_obj) %>%
    bind_rows(.id="lambda") %>%
    ggplot(aes(x=degree)) + geom_histogram() + facet_wrap(~lambda)
}
plot_scalefree_fit <- function(huge_obj) {
  get_degree_dists(huge_obj) %>% 
    map(unlist) %>%
    map(~try(WGCNA::scaleFreeFitIndex(.))) %>%
    bind_rows(.id="lambda") %>%
    ggplot(aes(x=lambda,y=Rsquared.SFT, group=1)) +
    geom_line() + 
    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1))
}

boot.huge <- function(dat, R=10, ...) {
  replicate(R, huge(dat[sample(1:nrow(dat), nrow(dat), TRUE)],  ...))
}


#################
# 1. DATA SETUP #
#################


#### 1a. SETUP GROUPS OF VARIABLES ####

genes_keep <- read.table("gene_sets/amigo_gene_list_unique.txt", #this is just the unique entries in amigo_gene_list.tsv
                         header = FALSE, stringsAsFactors = FALSE) %>% 
  pull %>% 
  unique %>% 
  `[`(. %in% featureNames(rna_clr))

#g is the genes
#CLR transforming the expression data, with all genes in the denominator
g_geom_mean = apply(exprs(rna),2, function(x) exp(mean(log1p(x))))
g<- t(log1p(sweep(exprs(rna_infl), 2, g_geom_mean[sampleNames(rna_infl)], "/")))

#m is the microbiome data
pseq_network <- filter_taxa(pseq_subset, function(x) sum(x>3) > 7, TRUE)
m <- compositions::clr(t(otu_table(pseq_network)@.Data))

#Z is a vector of confounders to residualize over
z_vars <- unique(unlist(confounders))
z <- data.frame(sample_data(pseq_subset))[z_vars] %>% 
  mutate(edu_w4=as.numeric(factor(edu_w4, levels=c("le_hs","some_college","bach","ge_masters")))) %>%
  mutate_all(scale) %>%
  mutate_all(as.vector)


#### 1b. PERFORM SINGLE IMPUTATION ON MISSING DATA ####

set.seed(10)
m1 <- mice::complete(mice::mice(cbind(z),m = 1))
z <- m1[z_vars]


#### 1c. ASSEMBLE DATA SET ####

#filter and combine genes and taxa
net_dat <- cbind(m[,colSums(m)>0],g[,colSums(g)>0])

#keep features expressed in greater than 10 individuals
net_dat <- net_dat[,colSums(abs(net_dat)>0)>10]

#residualize
net_dat_RESID <- apply(net_dat, 2, function(i)
  resid(lm(i~.,data=z)))

#perform nonparanormal transform
set.seed(41)
net_dat_RESID_npn <- huge.npn(net_dat_RESID)


#### 1d. SPLIT DATA ####

set.seed(51)
in_train <- sample(1:nrow(net_dat_RESID_npn), size = 177)

net_dat_tr <- net_dat_RESID_npn[in_train,]
net_dat_te <- net_dat_RESID_npn[-in_train,]


########################
# 2. FIT NETWORK MODEL #
########################

#### 2a. FIT NETWORK MODEL TO TRAINING DATA

set.seed(152)
huge_fit_tr <- huge(net_dat_tr, lambda=seq(0.25,0.39,length.out = 60), method = "mb")
#check fit to a scale free distribution and adjust lambda above to include only the range with sufficiently large R^2
plot_scalefree_fit(huge_fit_tr)
#then perform data-based model selection
huge_select_tr <- huge.select(huge_fit_tr, criterion = "stars", 
                              stars.thresh = 0.01, rep.num = 100)
#save, because this can take a little while
save(huge_select_tr, file= "huge_select_tr_R100.Rda")


#### 2b. FIT TO TEST DATA USING SAME LAMBDA

load("~/huge_select_tr_R100.Rda")
huge_fit_te <- huge(net_dat_te, lambda=huge_select_tr$opt.lambda)
huge_select_te <- huge_fit_te$path[[1]]


#### 2c. GET PROPORTION OF SHARED EDGES between training and test, by type

v_types <- ifelse(grepl("^[A-Z]", colnames(net_dat)), "gene",
                  ifelse(grepl("^[a-z]__", colnames(net_dat)), "otu", "phen"))
tr_edges <- as_edgelist(SpiecEasi::adj2igraph(huge_select_tr$refit))
te_edges <- as_edgelist(SpiecEasi::adj2igraph(huge_select_te))

all_edges <- list(train = tr_edges,
                  test = te_edges) %>%
  map(as.data.frame) %>%
  map(set_names, c("v1","v2")) %>%
  map(tibble::rownames_to_column,"eid") %>%
  bind_rows(.id="run") %>%
  mutate(edge_symbol = paste(v1, "-", v2)) %>%
  mutate(v1_type=v_types[v1], v2_type=v_types[v2], 
         edge_type=paste(v1_type, "-", v2_type))

unique_edges <- unique(all_edges[,c("edge_symbol", "edge_type")]) %>%
  mutate(in_train = edge_symbol %in% all_edges[all_edges$run=="train","edge_symbol"],
         in_test = edge_symbol %in% all_edges[all_edges$run=="test","edge_symbol"])

unique_edges %>% group_by(edge_type) %>% 
  summarise(
    n_replicated = sum(in_train & in_test),
    overlap=mean(in_train & in_test),
    replication=overlap / mean(in_train))

all_edges %>%
  filter(run=="train") %>%
  mutate(replicated=edge_symbol %in% all_edges[all_edges$run=="test","edge_symbol"])


#### 2d. FINAL NETWORK WITH ONLY REPLICATED EDGES ####

all_edges %>%
  filter(run=="train") %>%
  mutate(replicated=edge_symbol %in% all_edges[all_edges$run=="test","edge_symbol"]) %>%
  group_by(edge_type) %>% 
  summarise(rep_rate=mean(replicated))

edge_tbl <- all_edges %>%
  filter(run=="train") %>%
  mutate(replicated=edge_symbol %in% all_edges[all_edges$run=="test","edge_symbol"]) %>%
  filter(replicated==TRUE) 

eid_final <- edge_tbl$eid

#add node type (taxon vs. gene) and name to nodes
final_graph <- SpiecEasi::adj2igraph(huge_select_tr$refit)
vertex_attr(final_graph, "name") <- colnames(net_dat_tr)
final_graph <- igraph::delete.edges(final_graph,setdiff( E(final_graph), eid_final))
final_graph <- igraph::delete.vertices(final_graph, igraph::degree(final_graph)==0)
vertex_attr(final_graph, "shape") <- ifelse(grepl("^[a-z]__", vertex_attr(final_graph, "name")), 
                                            "circle", #taxa are circles, genes are squares
                                            ifelse(grepl("^[A-Z]", vertex_attr(final_graph, "name")),
                                                   "square",
                                                   "crectangle"))
#estimate covariances
set.seed(25)
huge_cov <- huge(net_dat_RESID_npn, lambda = huge_select_tr$opt.lambda-.1,
                 method="glasso", cov.output = TRUE)$cov[[1]] %>% 
  set_rownames(colnames(net_dat_RESID_npn)) %>%
  set_colnames(colnames(net_dat_RESID_npn))
cov_edges <- reshape2::melt(huge_cov) %>%
{set_names(.$value, paste0(.$Var1, "|", .$Var2))}


#### 2e. ASSEMBLE ####
all_graphs <- list(Discovery=SpiecEasi::adj2igraph(huge_select_tr$refit),
                   Replication= SpiecEasi::adj2igraph(huge_fit_te$path[[1]]),
                   Final=final_graph); V(all_graphs$Discovery)$name <-
  V(all_graphs$Replication)$name <- 
  V(all_graphs$Replication)$name <- colnames(net_dat_RESID_npn)

##################################
# 3. LOUVAIN COMMUNITY DETECTION #
##################################

set.seed(25)
louvain <- cluster_louvain(final_graph)
df_louvain <- data.frame(cluster_id=as.character(louvain$membership),
                         v=louvain$names) %>%
  mutate(v_type=ifelse(grepl("^[a-z]_", v), "OTU","Gene")) %>%
  group_by(cluster_id) %>%
  mutate(c_type=case_when(mean(v_type=="OTU")==1 ~ "B",
                          mean(v_type=="Gene")==1 ~ "C",
                          TRUE ~ "A"),
         c_size=n()) %>%
  ungroup %>%
  split(.$c_type) %>%
  map(arrange, desc(c_size)) %>%
  map(~mutate(., c_name=paste0(c_type, 
                               cumsum(!duplicated(cluster_id))))) %>%
  bind_rows()

cluster_membership <- set_names(df_louvain$c_name, df_louvain$v)

#save all network results
save(all_graphs,huge_cov,cov_edges, louvain, df_louvain, cluster_membership,
     pseq_network,file = "network.Rda")
