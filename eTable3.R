load("network.Rda")
library(tidyverse)

#supplemental table 3: all vertices
graph_for_eTable3 <- plot_graph %>%
  set_vertex_attr("name", value=V(.)$label)

#e <- as_data_frame(plot_graph)

eTable3 <- as_data_frame(tbl_graph, "vertices") %>%
  arrange(label_num) %>%
  mutate(list_edges = map(name, function(v)   as_data_frame(tbl_graph) %>%
                            filter(from==v | to==v) %>% 
                            select(from, to) %>% 
                            unlist %>% 
                            unique %>% 
                            `[`(-which(.==v))),
         edges = map_chr(list_edges, paste, collapse=", "),
         shape = recode(shape, OTU="circle",Gene="square")) %>%
  select(-list_edges,-label,-label_gene, -degree, Number_in_figure=label_num)

write.csv(eTable3, file = "csv/eTable3.csv", quote=FALSE, row.names = FALSE)
