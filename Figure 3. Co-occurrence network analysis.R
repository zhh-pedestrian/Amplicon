# Set the encoding format to UTF-8.

### 1 Correlation Analysis-----------------------
otu <- read.csv(file = "./0 otu.csv", row.names = 1)
columns_with_R0 <- grep("R0", names(otu), value = TRUE)
otu_rank0 <- (otu[, columns_with_R0])
otu_rank0 <- t(as.matrix(otu_rank0))

occor1 <- corAndPvalue(otu_rank0, y = NULL)
occor2 <- rcorr(otu_rank0)
occor_p <- occor1$p
occor_r <- occor2$r
occor_r[occor_p > 0.01 | abs(occor_r) < 0.7] = 0
occor_r[is.na(occor_r)] <- 0
occor_r[is.nan(occor_r)] <- 0
occor_r[is.infinite(occor_r)] <- 0
write.csv(occor_r,file="./Rank0.csv", fileEncoding = "UTF-8")

igraph <- graph_from_adjacency_matrix(occor_r, mode="undirected",weighted=TRUE,diag=FALSE)
bad.vs = V(igraph)[degree(igraph) == 0]
igraph = delete.vertices(igraph, bad.vs)
igraph.weight = E(igraph)$weight

## 2 Calculate network_properties-------------------
positive_correlations = sum(E(igraph)$weight > 0)
negative_correlations = sum(E(igraph)$weight < 0)
num_edges = length(E(igraph))
num_vertices = length(V(igraph))
connectance = edge_density(igraph, loops = FALSE)
average_degree = mean(degree(igraph))
if (any(E(igraph)$weight < 0)) {
  E(igraph)$weight <- abs(E(igraph)$weight)
}
average_path_length = mean_distance(igraph, weights = E(igraph)$weight)
diameter = diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
edge_connectivity = edge_connectivity(igraph)
clustering_coefficient = transitivity(igraph)
no_clusters = no.clusters(igraph)
centralization_betweenness = centralization.betweenness(igraph)$centralization
centralization_degree = centralization.degree(igraph)$centralization

# 创建数据框
network_properties <- data.frame(
  English = c("Positive correlations", "Negative correlations", "Number of edges", "Number of vertices", 
              "Connectance", "Average degree", "Average path length", "Diameter", 
              "Edge connectivity", "Clustering coefficient", "Number of clusters", 
  Value = c(positive_correlations, negative_correlations, num_edges, num_vertices, 
            connectance, average_degree, average_path_length, diameter, 
            edge_connectivity, clustering_coefficient, no_clusters, 
            centralization_betweenness, centralization_degree)
)

library(openxlsx)
write.xlsx(network_properties, file = "./rank0_network_properties.xlsx")


### 3 Calculate nodes and edges----------------------------------------
options(max.print = 3000000)
igraph = graph_from_adjacency_matrix(occor_r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph

set.seed(123)
k <- 3
cluster <- as.data.frame(kmeans(t(otu_rank0), centers = k)$cluster)
colnames(cluster)[1] <- "clu"
cluster$clu <- paste0("c_", cluster$clu)
nodes_df <- as.data.frame(V(igraph))
node <- merge(nodes_df, cluster, by = "row.names")
node$x <- NULL
colnames(node)[1] <- "ID"
node$Label <- node$ID
node <- node[c("ID", "Label", "clu")]
edge <- as_edgelist(igraph, names=T)
colnames(edge) <- c("Source", "Target")

write.csv(edge, file = "./rank0_igraph_edge.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(node, file = "./rank0_igraph_node.csv", row.names = FALSE, fileEncoding = "UTF-8")

### 4 Subsequently, import igraph for further beautification.---------------

