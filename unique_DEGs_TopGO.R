# Load necessary libraries
library(readxl)
library(readr)
library(topGO)
library(org.Gg.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)

# Read the DEG data for each organ
deg_data <- read.csv("pooled_DEGs_names.csv")
tail(deg_data)

Bursa_genes <- deg_data$Bursa[!is.na(deg_data$Bursa)]
Intestine_genes <- deg_data$Intestine[!is.na(deg_data$Intestine)]
Kidney_genes <- deg_data$Kidney[!is.na(deg_data$Kidney)]
Liver_genes <- deg_data$Liver[!is.na(deg_data$Liver)]

# Read counts data and set up gene universe
counts_data <- read_csv("final_counts_matrix.csv")
# rename the first column as Gene IDs
colnames(counts_data)[1] <- "Gene_ID"
head(counts_data)
gene_universe <- counts_data$Gene_ID
length(gene_universe)

# Create named vectors for each organ
Bursa_list <- factor(as.integer(gene_universe %in% Bursa_genes), levels = c(0, 1), labels = c("Not Significant", "Significant"))
names(Bursa_list) <- gene_universe
Intestine_list <- factor(as.integer(gene_universe %in% Intestine_genes))
names(Intestine_list) <- gene_universe
Kidney_list <- factor(as.integer(gene_universe %in% Kidney_genes))
names(Kidney_list) <- gene_universe
Liver_list <- factor(as.integer(gene_universe %in% Liver_genes))
names(Liver_list) <- gene_universe

# Function to run topGO analysis
run_topGO <- function(gene_list, gene_selection_function) {
  GO_data <- new("topGOdata",
                 ontology = "BP",  # Change to "MF" or "CC" as needed
                 allGenes = gene_list,
                 geneSel = gene_selection_function,
                 annot = annFUN.org,
                 mapping = "org.Gg.eg.db",
                 ID = "symbol")
  resultFisher <- runTest(GO_data, algorithm = "classic", statistic = "fisher")
  GenTable(GO_data, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 10)
}

gene_selection <- function(x) x == "Significant"

# Run topGO for each organ
Bursa_GO_results <- run_topGO(Bursa_list, gene_selection)
Intestine_GO_results <- run_topGO(Intestine_list, gene_selection)
Kidney_GO_results <- run_topGO(Kidney_list, gene_selection)
Liver_GO_results <- run_topGO(Liver_list, gene_selection)

# Save results
write.csv(Bursa_GO_results, "Bursa_GO_results.csv")
write.csv(Intestine_GO_results, "Intestine_GO_results.csv")
write.csv(Kidney_GO_results, "Kidney_GO_results.csv")
write.csv(Liver_GO_results, "Liver_GO_results.csv")

# Plotting
# Load GO results
Bursa_GO_BP_results <- read.csv("Bursa_GO_results.csv")
Intestine_GO_BP_results <- read.csv("Intestine_GO_results.csv")
Kidney_GO_BP_results <- read.csv("Kidney_GO_results.csv")
Liver_GO_BP_results <- read.csv("Liver_GO_results.csv")

# Calculate ratios and bind data for plotting
calculate_observed_expected_ratio <- function(GO_results, gene_universe_size, significant_genes_count) {
  GO_results$Annotated <- as.numeric(GO_results$Annotated)
  GO_results$Significant <- as.numeric(GO_results$Significant)
  GO_results$GeneRatio <- GO_results$Significant / GO_results$Annotated
  return(GO_results)
}

gene_universe_size <- length(gene_universe)
combined_GO_results <- bind_rows(
  calculate_observed_expected_ratio(Bursa_GO_BP_results, gene_universe_size, length(Bursa_genes)) %>% mutate(organ = "Bursa"),
  calculate_observed_expected_ratio(Intestine_GO_BP_results, gene_universe_size, length(Intestine_genes)) %>% mutate(organ = "Intestine"),
  calculate_observed_expected_ratio(Kidney_GO_BP_results, gene_universe_size, length(Kidney_genes)) %>% mutate(organ = "Kidney"),
  calculate_observed_expected_ratio(Liver_GO_BP_results, gene_universe_size, length(Liver_genes)) %>% mutate(organ = "Liver")
)

combined_GO_results$Term <- str_to_sentence(combined_GO_results$Term)
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, " regulation ", " reg. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Regulation", "Reg.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Response", "Resp.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Defense", "Def.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Negative", "Neg.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "activit...", "act.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, " response ", " resp. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "dna", "DNA")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "plasma-membrane a...", "plasma-membrane")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, ", postsyna...", "")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Very long-chain fatty acid metabolic pro...", "Long-chain fatty acid metabolic process")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "via plasma memb...", "")

# write combined results 
write.csv(combined_GO_results, "combined_GO_results.csv")


# Create a new color column based on specific strings in the 'Term' column
combined_GO_results <- combined_GO_results %>%
  mutate(color_group = case_when(
    str_detect(Term, "Def.|immune|Immune") ~ "Immunity",
    str_detect(Term, "Coagulation| coagulation|Wound|wound| hemostasis|Hemostasis| blood") ~ "Coagulation/Blood",
    TRUE ~ "Misc."
  ))

# Assign specific colors for each group
colors <- c("Immunity" = "orange", "Coagulation/Blood" = "#c83718", "Misc." = "#10cabb")

# Plot with the new color mapping
ggplot(combined_GO_results, aes(x = str_wrap(Term, width = 50), y = GeneRatio, color = color_group, size = GeneRatio)) +
  geom_point() +
  facet_wrap(~ organ, scales = "free_y") +
  coord_flip() +
  scale_color_manual(values = colors) +
  labs(x = "Pathways Enriched (Gene Ontology Biological Processes)", y = "Gene Ratio (Significant / Annotated)", color = "Pathway Category", size = "Gene Ratio") +
  theme_linedraw() +
  guides(color = guide_legend(order = 1), size = guide_legend(order = 2)) +
  # rename legend text not title
  theme(text = element_text(size = 12, family = "Times New Roman"),
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 16, vjust = 0.3),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 22),
        legend.text = element_text(size = 18, family = "Times New Roman"),
        legend.title = element_text(size = 20),
        # space between lines in legend 
        legend.key.height = unit(2, "lines"),
  )



# Save the plot
ggsave("GO_BP_dotplot_colored.svg", width = 12, height = 8, units = "in", dpi = 300)
