library(tidyverse)
load("network.Rda")

#link to a function that allows multiple color scales in one ggplot
source("https://gist.githubusercontent.com/eliocamp/eabafab2825779b88905954d84c82b32/raw/ab5a293508e632b1b7deb656eb7a0711d26f7716/new_aes.R")

#plot
set.seed(1)
net_plot <- plot_graph %>%
      ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_edges(aes(color=sign_cov), alpha=0.2) +
      scale_color_manual("Covariance", values=c("red","black")) +
      new_scale_color() +
      geom_nodes(aes(shape=shape, fill=cluster), size=7) +
      geom_nodetext(aes(label=label_num)) +
      geom_nodelabel(aes(label=label_gene, color=cluster),
                     alpha=0.85,
                     hjust="outward", vjust="outward") +
      scale_shape_manual(values=c(22,21), guide="none") +
      theme_blank()

svg(filename = "svg/eFigure2.svg", height = 12, width=12)
net_plot
dev.off()
