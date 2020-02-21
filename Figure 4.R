library(ggnetwork)
library(patchwork)
library(intergraph)
library(forcats)

load("network.Rda")

#############################
# 1. Assemble Data for Plot #
#############################

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

svg(filename = "svg/Figure4.svg", height=10, width=10); (figure4$A1 | figure4$A3) / (figure4$A2 | figure4$A4 |figure4$A5) + plot_annotation(tag_levels = "1", tag_prefix = "A"); dev.off()
