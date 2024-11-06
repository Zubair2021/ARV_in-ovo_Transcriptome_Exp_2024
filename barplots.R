library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)

# set working directory to current script's location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the data
counts <- read.csv("final_counts_matrix.csv", row.names = 1)
metadata <- read.csv("Final_Metadata.csv")

# Align the sample names between the counts matrix and the metadata
colnames(counts) <- sub("X", "", colnames(counts))

# confirm if the sample names are the same in the counts matrix and metadata and give both if or if else statements
if (all(colnames(counts) %in% metadata$sample_name)) {
  print("Sample names in counts matrix match metadata.")
} else {
  stop("Please check the sample names in the counts matrix and metadata.")
}

original_wd <- getwd()

for (tissue in unique(metadata$Tissue)) {
  # Create a directory for each tissue
  tissue_dir <- paste0(original_wd, "/Results_", tissue)
  if (!dir.exists(tissue_dir)) {
    dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  for (time_hpi in unique(metadata$hpi)) {
    cat("Analyzing:", tissue, time_hpi, "\n")
    
    # Filter samples based on current tissue and time_hpi
    current_metadata <- metadata[metadata$Tissue == tissue & metadata$hpi == time_hpi,]
    if (nrow(current_metadata) == 0) {
      cat("No samples for this group. Skipping...\n")
      next
    }
    current_samples <- current_metadata$sample_name
    current_counts <- counts[, current_samples, drop=FALSE]
    
    if (length(unique(current_metadata$Treatment)) < 2) {
      cat("Not enough treatment groups present. Skipping...\n")
      next
    }
    
    # Prepare a DGEList object for edgeR analysis
    y <- DGEList(counts = current_counts)
    y <- calcNormFactors(y)
    
    # Set NC as the reference level for the treatment variable
    current_metadata$Treatment <- relevel(factor(current_metadata$Treatment), ref = "NC")
    
    # Create the design matrix including the treatment condition
    design <- model.matrix(~ Treatment, data=current_metadata)
    
    # Estimate dispersion
    y <- estimateDisp(y, design)
    
    # Fit the model
    fit <- glmQLFit(y, design)
    
    # Define contrasts
    contrast_matrix <- makeContrasts(
      Low_vs_NC = TreatmentLow_S1133,
      High_vs_NC = TreatmentHigh_S1133,
      High_vs_Low = TreatmentHigh_S1133 - TreatmentLow_S1133,
      levels = design
    )
    
    # Conduct the tests for differential expression for each contrast
    contrasts <- colnames(contrast_matrix)
    for (contrast in contrasts) {
      cat("Testing contrast:", contrast, "\n")
      
      qlf <- glmQLFTest(fit, contrast=contrast_matrix[, contrast])
      
      # Get differentially expressed genes
      de_genes <- topTags(qlf, n=Inf)$table
      
      # Filter for significant DEGs
      sig_genes <- de_genes[abs(de_genes$logFC) > 1 & de_genes$FDR < 0.05, ]
      
      # Construct the output directory path for the timepoint
      timepoint_dir <- paste0(tissue_dir, "/", time_hpi)
      if (!dir.exists(timepoint_dir)) {
        dir.create(timepoint_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Now, construct the full path for the file including the output directory
      output_filename <- paste0(timepoint_dir, "/DEGs_", tissue, "_", time_hpi, "_", contrast, ".csv")
      
      # Save the DEGs to the file in the specified path
      if (nrow(sig_genes) > 0) {
        write.csv(sig_genes, file = output_filename, row.names = TRUE)
      }
    }
  }
}

# Read in all DEG files and create stacked barplots
results <- list()
for (tissue in unique(metadata$Tissue)) {
  tissue_dir <- paste0(original_wd, "/Results_", tissue)
  for (time_hpi in unique(metadata$hpi)) {
    for (contrast in c("Low_vs_NC", "High_vs_NC", "High_vs_Low")) {
      file_path <- paste0(tissue_dir, "/", time_hpi, "/DEGs_", tissue, "_", time_hpi, "_", contrast, ".csv")
      if (file.exists(file_path)) {
        deg_data <- read.csv(file_path)
        if (nrow(deg_data) > 0) {
          deg_data$tissue <- tissue
          deg_data$time_hpi <- time_hpi
          deg_data$contrast <- contrast
          results[[length(results) + 1]] <- deg_data
        }
      }
    }
  }
}

# Combine all results into one data frame
if (length(results) > 0) {
  all_results <- bind_rows(results)
  # all_results <- all_results %>% mutate(Direction = ifelse(logFC > 0, "Upregulated", "Downregulated"))
  
  plot <- ggplot(all_results, aes(x = time_hpi, fill = contrast)) +
    geom_bar(position = "stack") +
    facet_wrap(~ tissue, scales = "free") +
    labs(x = "Timepoint (hpi)", y = "Number of DEGs", fill = "Contrast") +
    theme_linedraw() + 
    # panel facet fill
    theme(strip.background = element_rect(color = "black", fill = "white"),
          strip.text = element_text(size = 20, family = "Times New Roman", color = "black", face = "bold"),
          strip.text.x = element_text(size = 20, family = "Times New Roman"),
          strip.text.y = element_text(size = 20, family = "Times New Roman"),
          strip.placement = "outside") +
    
    theme(
      text = element_text(size = 20, family = "Times New Roman"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.spacing = unit(2, "lines"),
    )
  print(plot)
  ggsave("figure/unique_DEGs_tissues.svg", plot, width = 12, height = 10, dpi = 300)
} else {
  cat("No significant DEGs found for any tissue or timepoint.\n")
}

# Save the plot to a file
dev.off()

