# Load necessary libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(DESeq2)

# set working directory to current script's location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Load data
counts_matrix <- read.csv("final_counts_matrix.csv", row.names = 1)
head(counts_matrix)
metadata <- read.csv("Final_metadata.csv")
# rename levels in the metadata
metadata$Treatment <- factor(metadata$Treatment, levels = c("NC", "Low_S1133", "High_S1133"), labels = c("NC", "Low S1133", "High S1133"))
head(metadata)


dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = metadata, design = ~ 1)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
counts_matrix_t <- assay(vsd)


# Transpose counts matrix for PCA (samples as rows, genes as columns)
counts_matrix_t <- t(counts_matrix_t)
head(counts_matrix_t)
# Remove columns with zero or near-zero variance
counts_matrix_t <- counts_matrix_t[, apply(counts_matrix_t, 2, var) > 0]

# Run PCA
pca_result <- prcomp(counts_matrix_t, scale. = TRUE)

# Create a data frame with PCA results and metadata
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample_ID <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by.x = "Sample_ID", by.y = "sample_name")


ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = Treatment)) +
  geom_point(aes(size = factor(hpi))) +  # Set size inside geom_point to avoid interference with stat_ellipse
  theme_linedraw() +
  stat_ellipse(aes(group = Tissue), linetype = "dashed") +  # Set group specifically for ellipses
  theme(strip.text = element_text(size = 28, family = "Times New Roman")) +
  theme(panel.spacing = unit(2, "lines")) +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  theme(text = element_text(size = 20, family = "Times New Roman")) +
  theme(legend.key.size = unit(1.5, "cm")) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),
    shape = guide_legend(override.aes = list(size = 5))
  ) +
  scale_shape_discrete(name = "Treatment") +
  scale_size_manual(name = "Time \nPost-Inoculation", values = c("24hpi" = 2, "48hpi" = 3)) +
  scale_color_discrete(name = "Tissue") 


ggsave("PCA_plot_tissue_wise.svg", width = 12, height = 8)


# cell type faceted graphs
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Tissue)) +
  geom_point(size = 3) +
  theme_linedraw(base_family = "Times New Roman") +
  stat_ellipse(aes(group = Treatment), level = 0.90, linetype = "dashed") +
  facet_grid(hpi~Tissue, scales = "free") +
  # change font size of the grid text
  theme(strip.text = element_text(size = 24, family="Times New Roman")) +
  # increase panel spacing
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = "Principal Component 1", y = "Principal Component 2") +
  theme(text = element_text(size = 20, family="Times New Roman")) +
  # adjust legend key size
  theme(legend.key.size = unit(1, "cm")) +  # This controls spacing
  # control the size of points in the legend
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # Adjust color legend point size
    shape = guide_legend(override.aes = list(size = 5))  # Adjust shape legend point siz   # Adjust size legend point size
  ) +
  # color the treatments manually by specifying colors 
  scale_color_manual(name = "Treatment", values = c("NC" = "seagreen", "Low S1133" = "blue3", "High S1133" = "#800")) +
  # scale_color_discrete(name = "Treatment") +
  #specify the size of the dots in the legend
  # scale_size_discrete(name = "Time (hpi)") +
  scale_shape_discrete(name = "Tissue")

ggsave("PCA_plot_tissue_time_faceted.svg", width = 12, height = 8)

library(plotly)
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, size = 3, color = ~Tissue, symbol = ~Treatment, type = 'scatter3d', mode = 'markers')



# permanova
# Calculate Euclidean distance matrix
distance_matrix <- dist(counts_matrix_t, method = "euclidean")

sink("pca_permanova_results.txt")

# PERMANOVA for Tissue
permanova_Tissue <- adonis2(distance_matrix ~ Tissue, data = metadata, permutations = 999)
print(permanova_Tissue)

# PERMANOVA for hpi
permanova_hpi <- adonis2(distance_matrix ~ hpi, data = metadata, permutations = 999)
print(permanova_hpi)

# PERMANOVA for Treatment
permanova_treatment <- adonis2(distance_matrix ~ Treatment, data = metadata, permutations = 999)
print(permanova_treatment)


# PERMANOVA for Treatment
permanova_overall <- adonis2(distance_matrix ~ Treatment + hpi + Tissue, data = metadata, permutations = 999)
print(permanova_overall)

sink()
