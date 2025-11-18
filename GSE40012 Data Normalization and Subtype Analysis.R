# Step 1: Data Inspection and Preprocessing

# Load libraries
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# Load dataset (Simple Method)
gse <- getGEO("GSE40012", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])

# Check for missing values
cat("Missing values in expression data:", sum(is.na(expression_data)), "\n")

# Boxplot before normalization
png(file = file.path("GSE40012_boxplot_before_normalization.png"))
boxplot(expression_data, main = "GSE40012 - Before Normalization", las = 2, outline = FALSE)
dev.off()

# Apply Normalization
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

# Boxplot after normalization
png(file = file.path("GSE40012_boxplot_after_normalization.png"))
boxplot(expression_data, main = "GSE40012 - After Normalization", las = 2, outline = FALSE)
dev.off()

# Filter low-expressed genes
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("Dimensions after filtering:", dim(expression_data), "\n")

# Extract Subtypes (Influenza, Bacterial, SIRS, Healthy)
subtype <- dplyr::case_when(
  grepl("Influenza", metadata$title, ignore.case = TRUE) & !grepl("Mixed", metadata$title, ignore.case = TRUE) ~ "Influenza A",
  grepl("Bacterial", metadata$title, ignore.case = TRUE) & !grepl("Mixed", metadata$title, ignore.case = TRUE) ~ "Bacterial Pneumonia",
  grepl("SIRS", metadata$title, ignore.case = TRUE) ~ "SIRS (Non-Infectious)",
  grepl("Healthy", metadata$title, ignore.case = TRUE) ~ "Healthy Control",
  TRUE ~ "Other"
)
subtype <- as.factor(subtype)

cat("Subtype distribution:\n")
table(subtype)

# PCA plot
pca_result <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Subtype = subtype)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point() +
  ggtitle("GSE40012 PCA Plot: Influenza vs Bacterial vs SIRS") +
  theme_minimal()

ggsave(file.path("GSE40012_pca_plot.png"), plot = pca_plot)

# Save processed data
save(expression_data, metadata, subtype, file = file.path("GSE40012_processed_data_step1.RData"))
