# Load necessary libraries
library(STRINGdb)
library(igraph)

# Initialize STRINGdb with Gallus gallus species ID
string_db <- STRINGdb$new(species=9031, version="11.5")

# Load DEGs data
deg_data <- read.csv("pooled_DEGs_names.csv")

# Iterate through each tissue column to analyze networks
for (tissue in colnames(deg_data)) {
  tissue_genes <- na.omit(deg_data[[tissue]])  # Get the genes for the tissue, removing NA values
  
  # Map genes to STRING IDs
  mapped_genes <- string_db$map(data.frame(Gene=tissue_genes), "Gene", removeUnmappedRows=TRUE)
  
  # Get interactions for mapped genes
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Create an igraph object from interactions
  g <- graph_from_data_frame(interactions[, c("from", "to")], directed=FALSE)
  
  # Plot the network
  plot(g, main = paste("Network for", tissue), vertex.label.cex = 0.7, vertex.size = 5, edge.width = 0.5)
}
