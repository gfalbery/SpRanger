
# Deriving multiple network stats ####

AllNetworkStats <- function(graph){

  library(igraph); library(tidyverse); library(ggregplot)

  list(

    Degree = degree(graph),

    Eigen = eigen_centrality(graph),

    #Between = betweenness(graph),

    #Closeness = closenness(graph),

    Trans = transitivity(graph),

    Components = components(graph)

    #Prev = Prev(graph)

  ) %>% return
}
