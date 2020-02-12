library(ggnetwork)
library(patchwork)
library(intergraph)
library(forcats)

load("network.Rda")

#############################
# 1. Assemble Data for Plot #
#############################

plot_graph <- all_graphs$Final

#attach shapes
V(plot_graph)$shape<- ifelse(grepl("^[a-z]__", V(plot_graph)$name), 
                             "OTU", "Gene")

#attach louvain clusters
V(plot_graph)$cluster <- cluster_membership[V(plot_graph)$name]

#attach covariances
E(plot_graph)$cov <- formatC(cov_edges[attr(E(plot_graph),"vnames")],digits=2,format="f")
E(plot_graph)$sign_cov <- ifelse(sign(as.numeric(E(plot_graph)$cov))==-1,"neg","pos")

#attach degree
V(plot_graph)$degree <- igraph::degree(plot_graph)

#create labels
df_labels <- data.frame(cl=V(plot_graph)$cluster, 
                        deg=V(plot_graph)$degree, 
                        type=V(plot_graph)$shape,
                        name=V(plot_graph)$name)%>%
  #create taxonomically meaninful labels for OTUs
  left_join(tax_table(pseq_network)@.Data %>% 
              gsub("[a-z]__", "", .) %>%
              as.data.frame %>% 
              tibble::rownames_to_column("name") %>%
              mutate(label = case_when(Species != "" ~ paste0(Genus, " ", Species),
                                       Genus != "" ~ paste0(Genus, " sp."),
                                       Family != "" ~ paste0("Family ", Family, " sp."),
                                       Order != "" ~ paste0("Order ", Order, " sp."))) %>%
              select(name, label)) %>%
  group_by(label) %>%
  #create roman numerals for duplicate named OTUs
  mutate(num=row_number(),
         rom=tolower(as.character(utils::as.roman(num))),
         n=n(), add_roman=n>1 & type=="OTU") %>%
  ungroup %>%
  mutate(label=ifelse(add_roman, paste0(label, " (", rom ,")"), label)) %>%
  #give OTUs numbers for graph labeling
  mutate(label_num=as.numeric(fct_reorder(label, deg, .desc=TRUE))) %>%
  #carry over gene name from 'name'
  mutate(label=ifelse(is.na(label), name, label)) %>%
  arrange(cl,desc(deg)) %>%
  select(-num:-add_roman) %>%
  set_rownames(.$name) 


#add labels to graph
V(plot_graph)$label = df_labels[V(plot_graph)$name, "label",drop=TRUE]
V(plot_graph)$label_num = df_labels[V(plot_graph)$name, "label_num",drop=TRUE]
V(plot_graph)$label_gene = ifelse(V(plot_graph)$shape=="Gene", 
                                  df_labels[V(plot_graph)$name, "label",drop=TRUE], NA)

#Keep a label if (a) it is a gene, (b) it shares an edge with a gene, (c) it shares 
#an edge with a node that shares an edge with a gene, or (d) it has top degree
subgraphs <- map(paste0("A",1:5) %>% set_names(.,.), 
                 ~subgraph(plot_graph, which(V(plot_graph)$cluster==.)))

#is connected to a gene by distance less than 3?
d <- map(subgraphs, distances, weights=NA)
d.conn_to_gene <- imap(d, ~.x[,V(subgraphs[[.y]])$shape=="Gene"]) %>% map(apply, 1, function(x) any(x<3))
d.direct_to_gene <- imap(d, ~.x[,V(subgraphs[[.y]])$shape=="Gene"]) %>% map(apply, 1, function(x) any(x<2))
subgraphs <- imap(subgraphs, ~set_vertex_attr(.x, "conn_to_gene", 1:gorder(.x), d.conn_to_gene[[.y]])) %>%
  imap(~set_vertex_attr(.x, "direct_to_gene", 1:gorder(.x), d.direct_to_gene[[.y]])) %>%
  map(~set_vertex_attr(., "label_gene", 1:gorder(.), ifelse(V(.)$shape=="Gene", V(.)$label,NA))) %>%
  map(~set_vertex_attr(., "label_otu", 1:gorder(.), ifelse(V(.)$shape=="OTU" &
                                                             (V(.)$conn_to_gene | 
                                                                V(.)$label %in%
                                                                na.omit(df_labels$labels_keep)),
                                                           V(.)$label, NA))) %>%
  map(~set_vertex_attr(., "label_num", 1:gorder(.), as.numeric(factor(V(.)$label_otu,
                                                                      levels=na.omit(V(.)$label_otu[order(-V(.)$direct_to_gene)]))))) %>%
  map(~set_vertex_attr(., "label", 1:gorder(.), ifelse(is.na(V(.)$label_gene), 
                                                       V(.)$label_num, V(.)$label_gene)))

##################
# 2. Create Plot #
##################

#create footnotes for plot abbreviations
plot_abbreviations <- map(subgraphs, ~data.frame(n=na.omit(V(.)$label_num), 
                                                 l=na.omit(V(.)$label_otu))) %>%
  map(arrange,n) %>%
  map(mutate, b=paste0(n, ": ", l)) %>%
  map(pull) %>%
  map(paste, collapse="; ")

#function to create cluster colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cluster_colors <- gg_color_hue(6)[1:5] %>% set_names(paste0("A",1:5))


#Figure 4
figure4 <- imap(subgraphs,                      
                     ~ggplot(.x,
                             aes(x = x, y = y, xend = xend, yend = yend)) + 
                       geom_edges(aes(color=sign_cov), alpha=0.5) +
                       geom_nodes(aes(shape=shape,size=degree), fill=cluster_colors[.y]) +
                       geom_nodetext(aes(label=label_num, alpha=direct_to_gene)) +
                       geom_nodelabel_repel(aes(label=label_gene),
                                            min.segment.length = 0,
                                            alpha=0.7) +
                       scale_size_continuous(range=c(2,10)) +
                       scale_shape_manual(values=c(22,21)) +
                       scale_color_manual(values=c("red","black")) +
                       scale_alpha_manual(values=c(0.4,1)) +
                       theme_blank() +
                       theme(legend.position = "none", plot.caption=element_text(hjust = 0)) +
                       labs(caption=stringr::str_wrap(plot_abbreviations[[.y]],
                                                      width = ifelse(.y %in% c("A1","A3"),80,50))))

svg(filename = "svg/net_subplots.svg", height=10, width=10); (figure4$A1 | figure4$A3) / (figure4$A2 | figure4$A4 |figure4$A5) + plot_annotation(tag_levels = "1", tag_prefix = "A"); dev.off()
