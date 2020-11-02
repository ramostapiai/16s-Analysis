##### CO-OCURRENCE Networks ####
#http://www.castrolab.org/isme/microbial_networks/microbial_networks.html#perfil-taxonomico
#Thanks to Katterinne Mendez for her help and doing this workflow for this.

library(SpiecEasi)
library(phyloseq)
setwd("path/")

#leer objetos phyloseq
readRDS("path/ps4.nasal.healthy.RDS") -> ps4.nasal.healthy
readRDS("path/ps4.nasal.asthmatic.RDS") -> ps4.nasal.asthmatic
readRDS("path/ps4.oral.healthy.RDS") -> ps4.oral.healthy
readRDS("path/ps4.oral.asthmatic.RDS") -> ps4.oral.asthmatic

# Infer the network using SpiecEasi from each phyloseq
# SpiecEasi has two inference models, which are defined with the 'method' argument: neighborhood selection ('mb') and inverse covariance selection ('mb')
# Use the 'ncores' argument to indicate the number of processors available (if you have more than one)
ps4.nasal.healthy_mb <- spiec.easi(ps4.nasal.healthy, method='mb', lambda.min.ratio=1e-2, nlambda=20,
                    pulsar.params=list(rep.num=50, ncores=10))

ps4.nasal.asthmatic_mb <- spiec.easi(ps4.nasal.asthmatic, method='mb',lambda.min.ratio=1e-2, nlambda=20,
                    pulsar.params=list(rep.num=50, ncores=10))

ps4.oral.healthy_mb <- spiec.easi(ps4.oral.healthy, method='mb',lambda.min.ratio=1e-2, nlambda=20,
                    pulsar.params=list(rep.num=50, ncores=10))

ps4.oral.asthmatic_mb <- spiec.easi(ps4.oral.asthmatic, method='mb',lambda.min.ratio=1e-2, nlambda=20,
                    pulsar.params=list(rep.num=50, ncores=10))

saveRDS(ps4.nasal.healthy_mb, "ps4.nasal.healthy_mb.RDS")
saveRDS(ps4.nasal.asthmatic_mb, "ps4.nasal.asthmatic_mb.RDS")
saveRDS(ps4.oral.healthy_mb, "ps4.oral.healthy_mb.RDS")
saveRDS(ps4.oral.asthmatic_mb, "ps4.oral.asthmatic_mb.RDS")

#Analysis
library(microbiome)
library(genefilter)
library(SpiecEasi)
library(seqtime)
library(igraph)
library(qgraph)
library(ggnet)
library(RColorBrewer)
library(tidyverse)
library(grid)
library(gridExtra)
library(network)
library(sna)
library(ggplot2)

#Read networks
readRDS("path/ps4.nasal.healthy_mb.RDS") -> ps4.nasal.healthy_mb
readRDS("path/ps4.nasal.asthmatic_mb.RDS") -> ps4.nasal.asthmatic_mb
readRDS("path/ps4.oral.healthy_mb.RDS") -> ps4.oral.healthy_mb
readRDS("path/ps4.oral.asthmatic_mb.RDS") -> ps4.oral.asthmatic_mb
#Read phyloseq object
readRDS("path/ps4.nasal.healthy.RDS") -> ps4.nasal.healthy
readRDS("path/ps4.nasal.asthmatic.RDS") -> ps4.nasal.asthmatic
readRDS("path/ps4.oral.healthy.RDS") -> ps4.oral.healthy
readRDS("path/ps4.oral.asthmatic.RDS") -> ps4.oral.asthmatic

# We use the getRefit () and adj2igraph () functions from the SpiecEasi package, to extract the refit matrix from the se_mb object and build the microbial network from it, respectively.
# Build network from 'sparse adjacency matrix' or 'refit matrix'
# Add OTU names to rows and columns
# Create igraph objects
# We use the ID of the taxa to name the vertices or nodes of the network

nh_mb_net <- adj2igraph(getRefit(ps4.nasal.healthy_mb), 
                            rmEmptyNodes = TRUE, diag = FALSE, 
                            vertex.attr = list(name = taxa_names(ps4.nasal.healthy)))

na_mb_net <- adj2igraph(getRefit(ps4.nasal.asthmatic_mb),
                            rmEmptyNodes = TRUE, diag = FALSE, 
                            vertex.attr = list(name = taxa_names(ps4.nasal.asthmatic)))

oh_mb_net <- adj2igraph(getRefit(ps4.oral.healthy_mb), 
                            rmEmptyNodes = TRUE, diag = FALSE, 
                            vertex.attr = list(name = taxa_names(ps4.oral.healthy)))

oa_mb_net <- adj2igraph(getRefit(ps4.oral.asthmatic_mb), 
                            rmEmptyNodes = TRUE, diag = FALSE, 
                            vertex.attr = list(name = taxa_names(ps4.oral.asthmatic)))


# Graph the networks
plot_network(nh_mb_net, ps4.nasal.healthy, type = "taxa", color = "Phylum", shape = "Kingdom", label = NULL)
plot_network(na_mb_net, ps4.nasal.asthmatic, type = "taxa", color = "Phylum", shape = "Kingdom", label = NULL)
plot_network(oh_mb_net, ps4.oral.healthy, type = "taxa", color = "Phylum", shape = "Kingdom", label = NULL)
plot_network(oa_mb_net, ps4.oral.asthmatic, type = "taxa", color = "Phylum", shape = "Kingdom", label = NULL)

# Extract adjacency matrix
nh_mb_net_class <- as_adjacency_matrix(nh_mb_net, type = "both")
na_mb_net_class <- as_adjacency_matrix(na_mb_net, type = "both")
oh_mb_net_class <- as_adjacency_matrix(oh_mb_net, type = "both")
oa_mb_net_class <- as_adjacency_matrix(oa_mb_net, type = "both")

# Generate class object 'network'
#nh
nh_mb_net_class <- network(as.matrix(nh_mb_net_class), 
                     vertex.attrnames = taxa_names(ps4.nasal.healthy), 
                     matrix.type = "adjacency", directed = F)
ggnet2(nh_mb_net_class)

#na
na_mb_net_class <- network(as.matrix(na_mb_net_class), 
                         vertex.attrnames = taxa_names(ps4.nasal.asthmatic), 
                         matrix.type = "adjacency", directed = F)
ggnet2(na_mb_net_class)

#oh
oh_mb_net_class <- network(as.matrix(oh_mb_net_class), 
                         vertex.attrnames = taxa_names(ps4.oral.healthy), 
                         matrix.type = "adjacency", directed = F)
ggnet2(oh_mb_net_class)

#oh
oa_mb_net_class <- network(as.matrix(oa_mb_net_class), 
                         vertex.attrnames = taxa_names(ps4.oral.asthmatic), 
                         matrix.type = "adjacency", directed = F)
ggnet2(oa_mb_net_class)

# Extract regression coefficients for mb objects
nh_mb_net_betaMat <- as.matrix(symBeta(getOptBeta(ps4.nasal.healthy_mb)))
na_mb_net_betaMat <- as.matrix(symBeta(getOptBeta(ps4.nasal.asthmatic_mb)))
oh_mb_net_betaMat <- as.matrix(symBeta(getOptBeta(ps4.oral.healthy_mb)))
oa_mb_net_betaMat <- as.matrix(symBeta(getOptBeta(ps4.oral.asthmatic_mb)))

# Calculate the number of positive and negative edges in networks
#nh
nh_positive <- length(nh_mb_net_betaMat[nh_mb_net_betaMat>0])/2 
nh_negative <- length(nh_mb_net_betaMat[nh_mb_net_betaMat<0])/2 
nh_total <- length(nh_mb_net_betaMat[nh_mb_net_betaMat!=0])/2
#na
na_positive <- length(na_mb_net_betaMat[na_mb_net_betaMat>0])/2 
na_negative <- length(na_mb_net_betaMat[na_mb_net_betaMat<0])/2 
na_total <- length(na_mb_net_betaMat[na_mb_net_betaMat!=0])/2
#oh
oh_positive <- length(oh_mb_net_betaMat[oh_mb_net_betaMat>0])/2 
oh_negative <- length(oh_mb_net_betaMat[oh_mb_net_betaMat<0])/2 
oh_total <- length(oh_mb_net_betaMat[oh_mb_net_betaMat!=0])/2
#oa
oa_positive <- length(oa_mb_net_betaMat[oa_mb_net_betaMat>0])/2 
oa_negative <- length(oa_mb_net_betaMat[oa_mb_net_betaMat<0])/2 
oa_total <- length(oa_mb_net_betaMat[oa_mb_net_betaMat!=0])/2

# We divide by 2 since the edge is represented by two entries in the matrix
saveRDS(ps4.nasal.healthy, "ps4.nasal.healthy.RDS")
saveRDS(ps4.nasal.asthmatic, "ps4.nasal.asthmatic.RDS")
saveRDS(ps4.oral.healthy, "ps4.oral.healthy.RDS")
saveRDS(ps4.oral.asthmatic, "ps4.oral.asthmatic.RDS")

# The first step is to extract the signs of the regression coefficients from the matrix of regression coefficients
#Change variables for each network (tax_ids, edges, net and betamat)
#nh
tax_ids <- taxa_names(ps4.nasal.healthy)
edges <- E(nh_mb_net)# edges
net <- nh_mb_net
betaMat <- nh_mb_net_betaMat
#Funcion aisgnar relaciones
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"forestgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red") # negative
  }
}
E(nh_mb_net)$color <- edge_colors

#na
tax_ids <- taxa_names(ps4.nasal.asthmatic)
edges <- E(na_mb_net)# edges
net <- na_mb_net
betaMat <- na_mb_net_betaMat
#Funcion aisgnar relaciones
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"forestgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red") # negative
  }
}
E(na_mb_net)$color <- edge_colors

#oh
tax_ids <- taxa_names(ps4.oral.healthy)
edges <- E(oh_mb_net)# edges
net <- oh_mb_net
betaMat <- oh_mb_net_betaMat
#Funcion aisgnar relaciones
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"forestgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red") # negative
  }
}
E(oh_mb_net)$color <- edge_colors

#oh
tax_ids <- taxa_names(ps4.oral.asthmatic)
edges <- E(oa_mb_net)# edges
net <- oa_mb_net
betaMat <- oa_mb_net_betaMat
#Funcion aisgnar relaciones
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"forestgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red") # negative
  }
}
E(oa_mb_net)$color <- edge_colors

#Generate the network class object and graph using ggnet2.
# Extract adjacency matrix
nh_net_class <- as_adjacency_matrix(nh_mb_net, type = "both")
na_net_class <- as_adjacency_matrix(na_mb_net, type = "both")
oh_net_class <- as_adjacency_matrix(oh_mb_net, type = "both")
oa_net_class <- as_adjacency_matrix(oa_mb_net, type = "both")

# Generate class object 'network'_mb
nh_net_class <- network(as.matrix(nh_net_class), 
                     vertex.attrnames = taxa_names(ps4.nasal.asthmatic), 
                     matrix.type = "adjacency", directed = F)
na_net_class <- network(as.matrix(na_net_class), 
                     vertex.attrnames = taxa_names(ps4.nasal.asthmatic), 
                     matrix.type = "adjacency", directed = F)
oh_net_class <- network(as.matrix(oh_net_class), 
                     vertex.attrnames = taxa_names(ps4.oral.healthy), 
                     matrix.type = "adjacency", directed = F)
oa_net_class <- network(as.matrix(oa_net_class), 
                     vertex.attrnames = taxa_names(ps4.oral.asthmatic), 
                     matrix.type = "adjacency", directed = F)
# We use the edge.color argument of the ggnet2 () function to color the edges.
# Graph network
ggnet2(nh_net_class, edge.color = E(nh_mb_net)$color)
ggnet2(na_net_class, edge.color = E(na_mb_net)$color)
ggnet2(oh_net_class, edge.color = E(oh_mb_net)$color)
ggnet2(oa_net_class, edge.color = E(oa_mb_net)$color)

#Add taxanomy to each node of the network
nh_mb_net2 <- nh_mb_net
na_mb_net2 <- na_mb_net
oh_mb_net2 <- oh_mb_net
oa_mb_net2 <- oa_mb_net

# Extract taxonomy table
nh_tax_tbl <- as.data.frame(ps4.nasal.healthy@tax_table@.Data)
na_tax_tbl <- as.data.frame(ps4.nasal.asthmatic@tax_table@.Data)
oh_tax_tbl <- as.data.frame(ps4.oral.healthy@tax_table@.Data)
oa_tax_tbl <- as.data.frame(ps4.oral.asthmatic@tax_table@.Data)
# Replace the name of each node with its phylum
nh_nodenames <- as.character(getTaxonomy(V(nh_mb_net)$name, nh_tax_tbl, level = "phylum", useRownames = TRUE))
na_nodenames <- as.character(getTaxonomy(V(na_mb_net)$name, na_tax_tbl, level = "phylum", useRownames = TRUE))
oh_nodenames <- as.character(getTaxonomy(V(oh_mb_net)$name, oh_tax_tbl, level = "phylum", useRownames = TRUE))
oa_nodenames <- as.character(getTaxonomy(V(oa_mb_net)$name, oa_tax_tbl, level = "phylum", useRownames = TRUE))

#Generate the network class object and graph using ggnet2.
# Extract adjacency matrix
nh_net_class <- as_adjacency_matrix(nh_mb_net, type = "both")
na_net_class <- as_adjacency_matrix(na_mb_net, type = "both")
oh_net_class <- as_adjacency_matrix(oh_mb_net, type = "both")
oa_net_class <- as_adjacency_matrix(oa_mb_net, type = "both")

# Generate class object 'network'_mb
nh_net_class <- network(as.matrix(nh_net_class), 
                        vertex.attrnames = taxa_names(ps4.nasal.healthy), 
                        matrix.type = "adjacency", directed = F)
na_net_class <- network(as.matrix(na_net_class), 
                        vertex.attrnames = taxa_names(ps4.nasal.asthmatic), 
                        matrix.type = "adjacency", directed = F)
oh_net_class <- network(as.matrix(oh_net_class), 
                        vertex.attrnames = taxa_names(ps4.oral.healthy), 
                        matrix.type = "adjacency", directed = F)
oa_net_class <- network(as.matrix(oa_net_class), 
                        vertex.attrnames = taxa_names(ps4.oral.asthmatic), 
                        matrix.type = "adjacency", directed = F)
# We use the edge.color argument of the ggnet2 () function to color the edges.
# Graph network
ggnet2(nh_net_class, color = nh_nodenames, edge.color = E(nh_mb_net)$color)
ggnet2(na_net_class, color = na_nodenames, edge.color = E(na_mb_net)$color)
ggnet2(oh_net_class, color = oh_nodenames, edge.color = E(oh_mb_net)$color)
ggnet2(oa_net_class, color = oa_nodenames, edge.color = E(oa_mb_net)$color)

#Fix the layout of the network.
# The coordinates must be in the form of an n vs. matrix. m, where n is the name of the nodes and m is the x and y coordinates of each node
# Assign coordinates to the network layout
nh_mb_net$layout <- array(1:40, dim = c(40, 8))
na_mb_net$layout <- array(1:40, dim = c(40, 8))
oh_mb_net$layout <- array(1:40, dim = c(40, 8))
oa_mb_net$layout <- array(1:40, dim = c(40, 8))
# Assign the layout as a fixed attribute of the network
nh_mb_net$layout <- layout.fruchterman.reingold(nh_mb_net)
na_mb_net$layout <- layout.fruchterman.reingold(na_mb_net)
oh_mb_net$layout <- layout.fruchterman.reingold(oh_mb_net)
oa_mb_net$layout <- layout.fruchterman.reingold(oa_mb_net)
# Graph network with fixed coordinates to be able to compare networks
ggnet2(nh_net_class, mode = nh_mb_net$layout)
ggnet2(na_net_class, mode = na_mb_net$layout)
ggnet2(oh_net_class, mode = oh_mb_net$layout)
ggnet2(oa_net_class, mode = oa_mb_net$layout)

#saving progress
#networks
saveRDS(nh_net_class, "path/nh_net_class.RDS")
saveRDS(na_net_class, "path/na_net_class.RDS")
saveRDS(oh_net_class, "path/oh_net_class.RDS")
saveRDS(oa_net_class, "path/oa_net_class.RDS")
#igraph
saveRDS(nh_mb_net, "path/nh_mb_net.RDS")
saveRDS(na_mb_net, "path/na_mb_net.RDS")
saveRDS(oh_mb_net, "path/oh_mb_net.RDS")
saveRDS(oa_mb_net, "path/oa_mb_net.RDS")

#Degree
#The degree value of a node represents the number of edges connected to the node in question, that is, the more relationships a particular node has with other nodes in the network, the greater degree.
#The distribution of degree values of the nodes that make up a network gives us an idea of the connectivity of the nodes in the network.

# Calculate node degree
#The degree value of a node represents the number of edges connected to the node in question
nh_deg <- igraph::degree(nh_mb_net, mode = "all")
na_deg <- igraph::degree(na_mb_net, mode = "all")
oh_deg <- igraph::degree(oh_mb_net, mode = "all")
oa_deg <- igraph::degree(oa_mb_net, mode = "all")
# Calcular degree distribution
nh_deg.dist <- degree_distribution(nh_mb_net, mode = "all", cumulative = F)
na_deg.dist <- degree_distribution(na_mb_net, mode = "all", cumulative = F)
oh_deg.dist <- degree_distribution(oh_mb_net, mode = "all", cumulative = F)
oa_deg.dist <- degree_distribution(oa_mb_net, mode = "all", cumulative = F)

# Graph degree distribution
#nh
plot(nh_deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(nh_deg.dist) 
#na
plot(na_deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(na_deg.dist)
#oh
plot(oh_deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(oh_deg.dist)
#oa
plot(oa_deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(oa_deg.dist)

#closeness
#Closeness refers to how central each node is to the entire network.
#To measure closeness, the network is randomly traversed and the frequency with which each node is visited is quantified.
#Nodes that are visited more frequently have a higher value of closeness.
nh_clos <- igraph::closeness(nh_mb_net, mode = "all")
na_clos <- igraph::closeness(na_mb_net, mode = "all")
oh_clos <- igraph::closeness(oh_mb_net, mode = "all")
oa_clos <- igraph::closeness(oa_mb_net, mode = "all")

#Transitivity, also called the clustering coefficient.
#Measures the probability that nodes adjacent to the node in question are connected to each other.

# We use the type = "global" argument to calculate the total transitivity value for the entire network
nh_clustering_coeff_global <- transitivity(nh_mb_net, type = "global")
na_clustering_coeff_global <- transitivity(na_mb_net, type = "global")
oh_clustering_coeff_global <- transitivity(oh_mb_net, type = "global")
oa_clustering_coeff_global <- transitivity(oa_mb_net, type = "global")
# We use the type = "local" argument to calculate the transitivity value of each node in the network
nh_clustering_coeff_local <- transitivity(nh_mb_net, type = "local")
na_clustering_coeff_local <- transitivity(na_mb_net, type = "local")
oh_clustering_coeff_local <- transitivity(oh_mb_net, type = "local")
oa_clustering_coeff_local <- transitivity(oa_mb_net, type = "local")

# Module detection
# This function tries to detect densely connected subnets, using 'random walks'
# 'random walks' refers to "walking" the network randomly
# It is assumed that "short routes" tend to stay in the same sub-network or module
nh_wt <- walktrap.community(nh_mb_net)
na_wt <- walktrap.community(na_mb_net)
oh_wt <- walktrap.community(oh_mb_net)
oa_wt <- walktrap.community(oa_mb_net)

# Visualize the hierarchical structure of the microbial community in a dendrogram
igraph::plot_dendrogram(nh_wt)
igraph::plot_dendrogram(na_wt)
igraph::plot_dendrogram(oh_wt)
igraph::plot_dendrogram(oa_wt)

# Calculate modularity
#Modularity is a good measure of how strong a network divides into modules or clusters.
#High modularity indicates that the network presents dense connections within certain groups of nodes (modules), and in turn, dispersed connections between groups of different nodes.
modularity(nh_mb_net, membership(nh_wt))
modularity(na_mb_net, membership(na_wt))
modularity(oh_mb_net, membership(oh_wt))
modularity(oa_mb_net, membership(oa_wt))

# Plot
plot(nh_wt, nh_mb_net)
plot(na_wt, na_mb_net)
plot(oh_wt, oh_mb_net)
plot(oa_wt, oa_mb_net)

# Extract adjacency matrix
nh_mb_net_class <- as_adjacency_matrix(nh_mb_net, type = "both")
na_mb_net_class <- as_adjacency_matrix(na_mb_net, type = "both")
oh_mb_net_class <- as_adjacency_matrix(oh_mb_net, type = "both")
oa_mb_net_class <- as_adjacency_matrix(oa_mb_net, type = "both")
# Generate class object 'network'
nh_mb_net_class <- network(as.matrix(nh_mb_net_class), vertex.attrnames = taxa_names(ps4.nasal.healthy), matrix.type = "adjacency", directed = F)
na_mb_net_class <- network(as.matrix(na_mb_net_class), vertex.attrnames = taxa_names(ps4.nasal.asthmatic), matrix.type = "adjacency", directed = F)
oh_mb_net_class <- network(as.matrix(oh_mb_net_class), vertex.attrnames = taxa_names(ps4.oral.healthy), matrix.type = "adjacency", directed = F)
oa_mb_net_class <- network(as.matrix(oa_mb_net_class), vertex.attrnames = taxa_names(ps4.oral.asthmatic), matrix.type = "adjacency", directed = F)

#plot
ggnet2(nh_mb_net_class, mode = nh_mb_net$layout, color = nh_wt$membership)
ggnet2(na_mb_net_class, mode = na_mb_net$layout, color = na_wt$membership)
ggnet2(oh_mb_net_class, mode = oh_mb_net$layout, color = oh_wt$membership)
ggnet2(oa_mb_net_class, mode = oa_mb_net$layout, color = oa_wt$membership)

# Name of the nodes
nh_nodenames
na_nodenames
oh_nodenames
oa_nodenames

# Search for keystone species (Node Degree + Node Centrality)
# degree of networks
nh_deg
na_deg
oh_deg
oa_deg

# Sort nodes according to degree from highest to lowest
nh_deg_sort <- sort(nh_deg, decreasing = TRUE)
na_deg_sort <- sort(na_deg, decreasing = TRUE)
oh_deg_sort <- sort(oh_deg, decreasing = TRUE)
oa_deg_sort <- sort(oa_deg, decreasing = TRUE)
#observar las primeras lineas 
head(nh_deg_sort)
head(na_deg_sort)
head(oh_deg_sort)
head(oa_deg_sort)

# Calculate betweenness of networks
#Betweenness is also a measure of centrality of the nodes that make up the network.
#The betweenness value of a node is calculated as the total number of shortest paths from all nodes to all other nodes that pass through the node in question.
nh_bn <- igraph::betweenness(nh_mb_net)
na_bn <- igraph::betweenness(na_mb_net)
oh_bn <- igraph::betweenness(oh_mb_net)
oa_bn <- igraph::betweenness(oa_mb_net)
# Sort nodes according to betweenness from highest to lowest
nh_bn_sort <- sort(nh_bn, decreasing = TRUE)
na_bn_sort <- sort(na_bn, decreasing = TRUE)
oh_bn_sort <- sort(oh_bn, decreasing = TRUE)
oa_bn_sort <- sort(oa_bn, decreasing = TRUE)

head(nh_bn_sort)
head(na_bn_sort)
head(oh_bn_sort)
head(oa_bn_sort)

#Transform to data.frame to plot the statistics
#Degree
nh_deg_sort_df <- as.data.frame(nh_deg_sort)
na_deg_sort_df <- as.data.frame(na_deg_sort)
oh_deg_sort_df <- as.data.frame(oh_deg_sort)
oa_deg_sort_df <- as.data.frame(oa_deg_sort)
#Betweenness
nh_bn_sort_df <- as.data.frame(nh_bn_sort)
na_bn_sort_df <- as.data.frame(na_bn_sort)
oh_bn_sort_df <- as.data.frame(oh_bn_sort)
oa_bn_sort_df <- as.data.frame(oa_bn_sort)

# Keystone species criteria
#nh
# We use the 'grid.arrange' function from the 'gridExtra' package to display both histograms together
pdeg_nh <- {ggplot(nh_deg_sort_df, aes(x = nh_deg_sort)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(aes(xintercept=mean(nh_deg_sort)),
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Degree", y = "Node count", title = "Network Nodes Degree nh")}
pbn_nh <- {ggplot(nh_bn_sort_df, aes(x = nh_bn_sort)) + 
  geom_histogram(binwidth = 30) +  
  geom_vline(aes(xintercept=mean(nh_bn_sort)), 
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness nh")}
grid.arrange(pdeg_nh, pbn_nh, ncol = 2)
#Degree > 4
#Betweenness > 100

#na
# We use the 'grid.arrange' function from the 'gridExtra' package to display both histograms together
pdeg_na <- {ggplot(na_deg_sort_df, aes(x = na_deg_sort)) + 
    scale_x_continuous(breaks=c(0,5,10,20,30,40,50,60)) + 
    geom_histogram(binwidth = 1) + 
    geom_vline(aes(xintercept=mean(na_deg_sort)),
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Degree", y = "Node count", title = "Network Nodes Degree na")}
pbn_na <- {ggplot(na_bn_sort_df, aes(x = na_bn_sort)) + 
    geom_histogram(binwidth = 30) +  
    geom_vline(aes(xintercept=mean(na_bn_sort)), 
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness na")}
grid.arrange(pdeg_na, pbn_na, ncol = 2)
#Degree > 4
#Betweenness > 100

#oh
# We use the 'grid.arrange' function from the 'gridExtra' package to display both histograms together
pdeg_oh <- {ggplot(oh_deg_sort_df, aes(x = oh_deg_sort)) + 
    scale_x_continuous(breaks=c(0,5,10,20,30,40,50,60)) + 
    geom_histogram(binwidth = 1) + 
    geom_vline(aes(xintercept=mean(oh_deg_sort)),
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Degree", y = "Node count", title = "Network Nodes Degree oh")}
pbn_oh <- {ggplot(oh_bn_sort_df, aes(x = oh_bn_sort)) + 
    geom_histogram(binwidth = 30) +  
    geom_vline(aes(xintercept=mean(oh_bn_sort)), 
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness oh")}
grid.arrange(pdeg_oh, pbn_oh, ncol = 2)
#Degree > 5
#Betweenness > 200

#oa
# We use the 'grid.arrange' function from the 'gridExtra' package to display both histograms together
pdeg_oa <- {ggplot(oa_deg_sort_df, aes(x = oa_deg_sort)) + 
    scale_x_continuous(breaks=c(0,5,10,20,30,40,50,60)) + 
    geom_histogram(binwidth = 1) + 
    geom_vline(aes(xintercept=mean(oa_deg_sort)),
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Degree", y = "Node count", title = "Network Nodes Degree oa")}
pbn_oa <- {ggplot(oa_bn_sort_df, aes(x = oa_bn_sort)) + 
    scale_x_continuous(breaks=c(0,100,200,300,400,500,1000)) + 
    geom_histogram(binwidth = 30) +  
    geom_vline(aes(xintercept=mean(oa_bn_sort)), 
               color="black", linetype="dashed", size=1) + 
    theme_minimal() + 
    labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness oa")}
grid.arrange(pdeg_oa, pbn_oa, ncol = 2)
#Degree > 5
#Betweenness > 500

#Degree
nh_deg_sort_df #Degree > 4
na_deg_sort_df #Degree > 4
oh_deg_sort_df #Degree > 5
oa_deg_sort_df #Degree > 5
#Betweenness
nh_bn_sort_df #Betweenness > 100
na_bn_sort_df #Betweenness > 100
oh_bn_sort_df #Betweenness > 200
oa_bn_sort_df #Betweenness > 500

#Add column with ASV (TaxID) Degree
nh_deg_sort_df$TaxID <- row.names(nh_deg_sort_df)
na_deg_sort_df$TaxID <- row.names(na_deg_sort_df)
oh_deg_sort_df$TaxID <- row.names(oh_deg_sort_df)
oa_deg_sort_df$TaxID <- row.names(oa_deg_sort_df)
#Add column with ASV (TaxID) Betweenness
nh_bn_sort_df$TaxID <- row.names(nh_bn_sort_df)
na_bn_sort_df$TaxID <- row.names(na_bn_sort_df)
oh_bn_sort_df$TaxID <- row.names(oh_bn_sort_df)
oa_bn_sort_df$TaxID <- row.names(oa_bn_sort_df)

# We filter the taxa with degree> X and save them in a new data frame
#Degree
nh_deg_high_df <- dplyr::filter(nh_deg_sort_df, nh_deg_sort > 4)
na_deg_high_df <- dplyr::filter(na_deg_sort_df, na_deg_sort > 4)
oh_deg_high_df <- dplyr::filter(oh_deg_sort_df, oh_deg_sort > 5)
oa_deg_high_df <- dplyr::filter(oa_deg_sort_df, oa_deg_sort > 5)
#Betweenness
nh_bn_high_df <- dplyr::filter(nh_bn_sort_df, nh_bn_sort > 100)
na_bn_high_df <- dplyr::filter(na_bn_sort_df, na_bn_sort > 100)
oh_bn_high_df <- dplyr::filter(oh_bn_sort_df, oh_bn_sort > 100)
oa_bn_high_df <- dplyr::filter(oa_bn_sort_df, oa_bn_sort > 500)

# We join the filtered degree and betweenness tables
nh_keystone <- merge(nh_deg_high_df, nh_bn_high_df, all.x = FALSE)
na_keystone <- merge(na_deg_high_df, na_bn_high_df, all.x = FALSE)
oh_keystone <- merge(oh_deg_high_df, oh_bn_high_df, all.x = FALSE)
oa_keystone <- merge(oa_deg_high_df, oa_bn_high_df, all.x = FALSE)

# We sort 'keystone' according to degree from highest to lowest
nh_keystone <- nh_keystone[order(nh_keystone$nh_deg_sort, decreasing = TRUE),]
na_keystone <- na_keystone[order(na_keystone$na_deg_sort, decreasing = TRUE),]
oh_keystone <- oh_keystone[order(oh_keystone$oh_deg_sort, decreasing = TRUE),]
oa_keystone <- oa_keystone[order(oa_keystone$oa_deg_sort, decreasing = TRUE),]
# Define the tax IDs as row names of the 'keystone' data frame
row.names(nh_keystone) <- nh_keystone$TaxID
row.names(na_keystone) <- na_keystone$TaxID
row.names(oh_keystone) <- oh_keystone$TaxID
row.names(oa_keystone) <- oa_keystone$TaxID

#Graph the co-occurrence networks including "keystone species"
ggnet2(nh_mb_net_class, mode = nh_mb_net$layout, color = nh_wt$membership)
ggnet2(na_mb_net_class, mode = na_mb_net$layout, color = na_wt$membership)
ggnet2(oh_mb_net_class, mode = oh_mb_net$layout, color = oh_wt$membership)
ggnet2(oa_mb_net_class, mode = oa_mb_net$layout, color = oa_wt$membership)

# Graficar red
# Definir paleta de colores por Phylum
colors1 <- c("Actinobacteria" = "#263238", 
             "Bacteroidetes" = "#03A9F4", 
             "Epsilonbacteraeota" = "#E91E63", 
             "Firmicutes" = "#673AB7", 
             "Fusobacteria" = "#FFC107",
             "Proteobacteria" = "#4CAF50")

#nh
ggnet2(nh_mb_net_class, mode = nh_mb_net$layout,
       node.alpha = 0.9,
       palette = colors1,
       color = nh_nodenames,
       edge.color = E(nh_mb_net)$color, 
       label = nh_keystone$TaxID, label.size = 4) + ggtitle("Nasal Healthy") -> final_nh
#na
ggnet2(na_mb_net_class, mode = na_mb_net$layout,
       node.alpha = 0.9,
       palette = colors1,
       color = na_nodenames,
       edge.color = E(na_mb_net)$color, 
       label = na_keystone$TaxID, label.size = 4) + ggtitle("Nasal Asthmatic") -> final_na
#oh
ggnet2(oh_mb_net_class, mode = oh_mb_net$layout,
       node.alpha = 0.9,
       palette = colors1,
       color = oh_nodenames,
       edge.color = E(oh_mb_net)$color, 
       label = oh_keystone$TaxID, label.size = 4) + ggtitle("Oral Healthy") -> final_oh
#oa
ggnet2(oa_mb_net_class, mode = oa_mb_net$layout,
       node.alpha = 0.9,
       palette = colors1,
       color = oa_nodenames,
       edge.color = E(oa_mb_net)$color, 
       label = oa_keystone$TaxID, label.size = 4) + ggtitle("Oral Asthmatic") -> final_oa

#Resumen plots
final_nh
final_na
final_oh
final_oa

multiplot(final_nh, final_oh, final_na, final_oa, cols=2)

#Names Keystone species
#nh
#ASV2059 "Prevotellaceae"     "Prevotella_2"            "conceptionensis"
#ASV2128 "Prevotellaceae"     "Prevotella_7"            "melaninogenica"

#na
#ASV27   "Leptotrichiaceae"   "Leptotrichia"            NA  
#AS60    "Leptotrichiaceae"   "Leptotrichia"            NA       
#ASV105  "Leptotrichiaceae"   "Leptotrichia"            "buccalis"
#ASV2007 "Porphyromonadaceae" "Porphyromonas"           NA
#ASV2081 "Prevotellaceae"     "Prevotella_6"            "salivae"
#ASV3054 "Neisseriaceae"      "Kingella"                NA          

#oh
#ASV203  "Streptococcaceae"   "Streptococcus"           NA   

#oa
#ASV203  "Streptococcaceae"   "Streptococcus"           NA       
#ASV296  "Carnobacteriaceae"  "Granulicatella"          "elegans"   
#ASV441  "Veillonellaceae"    "Veillonella"             NA

#“Determine if any co-occurrence relationship is consistent across ecosystems 
#graph.intersection.by.name – igraph function
#intersection edges

o_inter_net <- graph.intersection(oh_mb_net, oa_mb_net)
o_inter_net
n_inter_net <- graph.intersection(nh_mb_net, na_mb_net)
n_inter_net
ggnet2(n_inter_net)

##### CORE MICROBIOME #####
library(microbiome)

# Calculate compositional version of the data
# (relative abundances)
#Subset_samples in diferentes phyloseq_object
subset_samples(ps4, host_tissue_sampled%in%c("oral mucosa"))-> ps4.oral
subset_samples(ps4.nasal, host_disease%in%c("Healthy"))-> ps4.nasal.healthy
subset_samples(ps4.nasal, host_disease%in%c("Asthma"))-> ps4.nasal.asthmatic
subset_samples(ps4.oral, host_disease%in%c("Healthy"))-> ps4.oral.healthy
subset_samples(ps4.oral, host_disease%in%c("Asthma"))-> ps4.oral.asthmatic
sample_data(ps4.nasal)$Host_Age <- as.numeric(sample_data(ps4.nasal)$Host_Age)
subset_samples(ps4.nasal, sample_data(ps4.nasal)[["Host_Age"]]>14) -> ps4.nasal.teenagers
subset_samples(ps4.nasal, sample_data(ps4.nasal)[["Host_Age"]]<=14) -> ps4.nasal.children

subset_samples(ps4, collection_date%in%c("2017","2018"))-> ps4.case.control
subset_samples(ps4.case.control, host_tissue_sampled%in%c("oral mucosa"))-> ps4.case.control.oral
subset_samples(ps4.case.control.oral, host_disease%in%c("Healthy"))-> ps4.case.control.oral.h
subset_samples(ps4.case.control.oral, host_disease%in%c("Asthma"))-> ps4.case.control.oral.a
subset_samples(ps4.case.control, host_tissue_sampled%in%c("nasal mucosa"))-> ps4.case.control.nasal
subset_samples(ps4.case.control.nasal, host_disease%in%c("Healthy"))-> ps4.case.control.nasal.h
subset_samples(ps4.case.control.nasal, host_disease%in%c("Asthma"))-> ps4.case.control.nasal.a

nh_rel <- microbiome::transform(ps4.case.control.nasal.h, "compositional")
na_rel <- microbiome::transform(ps4.case.control.nasal.a, "compositional")
oh_rel <- microbiome::transform(ps4.case.control.oral.h, "compositional")
oa_rel <- microbiome::transform(ps4.case.control.oral.a, "compositional")

# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

###Plot nh ###
library(RColorBrewer)
p_nh <- plot_core(nh_rel, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = F)
p_nh + ggtitle("Core microbiome nasal healthy") -> p_nh
print(p_nh)

{# get the data used for plotting 
df <- p_nh$data 

# get the list of OTUs
list <- df$Taxa 

# check the OTU ids
# print(list) 

# get the taxonomy data
tax <- tax_table(ps4.case.control.nasal.h)
tax <- as.data.frame(tax)

# add the OTus to last column
tax$OTU <- rownames(tax)

# select taxonomy of only 
# those OTUs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

# head(tax2)

# We will merege all the column into one except the Doamin as all is bacteria in this case
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Genus"), sep = "_;", remove = TRUE)

tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

# add this new information into the plot data df

df$Taxa <- tax.unit$Taxa_level

# you can see now we have the taxonomic information
knitr::kable(head(df))

# replace the data in the plot object
p_nh$data <- df}

plot(p_nh + theme(axis.text.y = element_text(face="italic")))

###Plot na ###
library(RColorBrewer)
p_na <- plot_core(na_rel, plot.type = "heatmap", 
                  prevalences = prevalences,
                  detections = detections,
                  colours = rev(brewer.pal(5, "Spectral")),
                  min.prevalence = .2, horizontal = F)
p_na + ggtitle("Core microbiome nasal asthmatic") -> p_na
print(p_na)

# get the data used for plotting 
{df <- p_na$data 

# get the list of OTUs
list <- df$Taxa 

# check the OTU ids
# print(list) 

# get the taxonomy data
tax <- tax_table(ps4.case.control.nasal.a)
tax <- as.data.frame(tax)

# add the OTus to last column
tax$OTU <- rownames(tax)

# select taxonomy of only 
# those OTUs that are used in the plot
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

# head(tax2)

# We will merege all the column into one except the Doamin as all is bacteria in this case
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Genus"), sep = "_;", remove = TRUE)

tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)

# add this new information into the plot data df

df$Taxa <- tax.unit$Taxa_level

# you can see now we have the taxonomic information
knitr::kable(head(df))

# replace the data in the plot object
p_na$data <- df}

plot(p_na + theme(axis.text.y = element_text(face="italic")))

###Plot oh ###
library(RColorBrewer)
p_oh <- plot_core(oh_rel, plot.type = "heatmap", 
                  prevalences = prevalences,
                  detections = detections,
                  colours = rev(brewer.pal(5, "Spectral")),
                  min.prevalence = .2, horizontal = F)
p_oh + ggtitle("Core microbiome oral healthy") -> p_oh
print(p_oh)

{df <- p_oh$data 
  
  # get the list of OTUs
  list <- df$Taxa 
  
  # check the OTU ids
  # print(list) 
  
  # get the taxonomy data
  tax <- tax_table(ps4.case.control.oral.h)
  tax <- as.data.frame(tax)
  
  # add the OTus to last column
  tax$OTU <- rownames(tax)
  
  # select taxonomy of only 
  # those OTUs that are used in the plot
  tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
  
  # head(tax2)
  
  # We will merege all the column into one except the Doamin as all is bacteria in this case
  tax.unit <- tidyr::unite(tax2, Taxa_level,c("Genus"), sep = "_;", remove = TRUE)
  
  tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
  
  # add this new information into the plot data df
  
  df$Taxa <- tax.unit$Taxa_level
  
  # you can see now we have the taxonomic information
  knitr::kable(head(df))
  
  # replace the data in the plot object
  p_oh$data <- df}

plot(p_oh + theme(axis.text.y = element_text(face="italic")))

###Plot oa ###
library(RColorBrewer)
p_oa <- plot_core(oa_rel, plot.type = "heatmap", 
                  prevalences = prevalences,
                  detections = detections,
                  colours = rev(brewer.pal(5, "Spectral")),
                  min.prevalence = .2, horizontal = F)
p_oa + ggtitle("Core microbiome oral asthmatic") -> p_oa
print(p_oa)

{df <- p_oa$data 
  
  # get the list of OTUs
  list <- df$Taxa 
  
  # check the OTU ids
  # print(list) 
  
  # get the taxonomy data
  tax <- tax_table(ps4.case.control.oral.a)
  tax <- as.data.frame(tax)
  
  # add the OTus to last column
  tax$OTU <- rownames(tax)
  
  # select taxonomy of only 
  # those OTUs that are used in the plot
  tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
  
  # head(tax2)
  
  # We will merege all the column into one except the Doamin as all is bacteria in this case
  tax.unit <- tidyr::unite(tax2, Taxa_level,c("Genus"), sep = "_;", remove = TRUE)
  
  tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
  
  # add this new information into the plot data df
  
  df$Taxa <- tax.unit$Taxa_level
  
  # you can see now we have the taxonomic information
  knitr::kable(head(df))
  
  # replace the data in the plot object
  p_oa$data <- df}

plot(p_oa + theme(axis.text.y = element_text(face="italic")))

multiplot(p_nh, p_na, p_oh, p_oa,  cols = 2)

save.image(file = "path/co-ocurrence_analysis.RData")