library(data.table)
library(dplyr)
library(ggplot2)
library(UpSetR)
library(ComplexUpset)
library(limma)
# Set the path to the edgeR results directory
# set working directory to current script's location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

edgeR_featureCounts_dir <- getwd()

# Function to list and read DEG files for specific comparisons
gather_degs <- function(dir_path, tissues, timepoints, comparisons) {
  degs_list <- list()
  
  for (tissue in tissues) {
    tissue_dir <- paste0(dir_path, "/Results_", tissue)
    
    for (timepoint in timepoints) {
      for (comparison in comparisons) {
        file_path <- paste0(tissue_dir, "/", timepoint, "/DEGs_", tissue, "_", timepoint, "_", comparison, ".csv")
        if (file.exists(file_path)) {
          deg_data <- fread(file_path)
          # cat("Successfully read file:", file_path, "\n")
          # cat("Number of DEGs in", tissue, timepoint, comparison, ":", nrow(deg_data), "\n")
          degs_list[[paste(tissue, timepoint, comparison, sep = "_")]] <- deg_data[[1]]  # Assuming first column has gene names
        } else {
          cat("File not found:", file_path, "\n")
          degs_list[[paste(tissue, timepoint, comparison, sep = "_")]] <- NULL
        }
      }
    }
  }
  
  return(degs_list)
}

# Define parameters
tissues_of_interest <- c("Bursa", "Intestine", "Kidney", "Liver")  # Replace with the tissues you are interested in
timepoints_to_pool <- c("24hpi", "48hpi")
comparisons_to_consider <- c("Low_vs_NC", "High_vs_NC")

# Gather DEGs for the specified tissues, timepoints, and comparisons
degs_list <- gather_degs(edgeR_featureCounts_dir, tissues_of_interest, timepoints_to_pool, comparisons_to_consider)

# Pool DEGs for 24hpi and 48hpi for each tissue before finding the union across comparisons
pooled_degs_list <- list()
for (tissue in tissues_of_interest) {
  tissue_degs <- list()
  
  for (timepoint in timepoints_to_pool) {
    for (comparison in comparisons_to_consider) {
      key <- paste(tissue, timepoint, comparison, sep = "_")
      if (!is.null(degs_list[[key]])) {
        tissue_degs <- union(tissue_degs, degs_list[[key]])
        cat("Number of DEGs in", key, length(tissue_degs), "\n")
      }
    }
  }
  pooled_degs_list[[tissue]] <- unique(tissue_degs)
  cat("\n --- Number of DEGs in", tissue, "after pooling timepoints:", length(pooled_degs_list[[tissue]]), "\n\n")
}

# Combine DEGs from all tissues (union across tissues)
all_degs <- unique(unlist(pooled_degs_list))
test_common_degs <- Reduce(intersect, pooled_degs_list)
cat("\nTotal unique DEGs after pooling across tissues:", length(all_degs), "\n\n")
cat("Total common DEGs across all tissues:", length(test_common_degs), "\n\n",
    "Common DEGs:", paste(test_common_degs, collapse = ", "), "\n\n")

# Create a list for UpSet plot, indicating presence/absence in each tissue
degs_upset_list <- list()
for (tissue in names(pooled_degs_list)) {
  degs_upset_list[[tissue]] <- pooled_degs_list[[tissue]]
}



upset_data <- fromList(degs_upset_list)

detach("package:UpSetR", unload = TRUE)

upset(
  upset_data,
  name = "",
  intersect = names(degs_upset_list),
  sort_intersections = "descending",
  base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE, fill = "purple", text = list(size = 6, family = "Times New Roman")) +
      theme(text = element_text(family = "Times New Roman", size = 24))
  ),
  queries = list(
    upset_query(
      intersect = c('Intestine', 'Liver', 'Kidney', 'Bursa'),
      color = 'red',
      fill = 'red',
      only_components = c('intersections_matrix', 'Intersection Size'),
      min_size = 1,
      stripes = "lightgray"
    )
  ),
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(text=element_text(size=24, family="Times New Roman")),
      'overall_sizes'=theme(text=element_text(family="Times New Roman", size=24))
    )
  )
) 


# save svg
ggsave("figures/upset_plot.svg", width = 12, height = 9, units = "in", dpi = 300)

# Function to write pooled DEGs to a wide format CSV with DEG names under each tissue column
write_pooled_degs_names_csv <- function(pooled_degs, output_path) {
  # Find the maximum number of DEGs for any tissue to set the number of rows
  max_degs <- max(sapply(pooled_degs, length))
  
  # Create a data frame with each tissue as a column and NA for empty cells
  degs_df <- data.table(matrix(NA, nrow = max_degs, ncol = length(pooled_degs)))
  setnames(degs_df, names(pooled_degs))
  
  # Fill each column with the respective DEGs, leaving excess rows as NA
  for (tissue in names(pooled_degs)) {
    degs_df[[tissue]][1:length(pooled_degs[[tissue]])] <- pooled_degs[[tissue]]
    cat("\nThe number of pooled DEGs for", tissue, length(pooled_degs[[tissue]]), "\n\n")
  }
  
  # Write the data frame to CSV
  fwrite(degs_df, file = output_path, row.names = FALSE)
  cat("Pooled DEGs written to CSV with DEG names for each tissue at:", output_path, "\n")
}

# Define the output CSV file path
output_csv_path <- file.path(edgeR_featureCounts_dir, "pooled_DEGs_names.csv")

# Write the pooled DEGs with names in each tissue column
write_pooled_degs_names_csv(pooled_degs_list, output_csv_path)






