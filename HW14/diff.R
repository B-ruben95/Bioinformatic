# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Load the counts matrix and design file
counts <- read.csv("res/counts-hisat.csv", row.names = 1)
design <- read.csv("design_diff.csv", stringsAsFactors = T)

# Ensure the design file is properly formatted
rownames(design) <- design$sample

design <- design[, -which(colnames(design) == "sample")] # Remove Sample column if present

# Check if the rownames of design match the colnames of the counts matrix
if (!all(colnames(counts) %in% rownames(design))) {
  cat("Sample names in the counts matrix:", colnames(counts), "\n")
  cat("Sample names in the design file:", rownames(design), "\n")
  
  # Find mismatched samples
  missing_in_design <- setdiff(colnames(counts), rownames(design))
  missing_in_counts <- setdiff(rownames(design), colnames(counts))
  
  if (length(missing_in_design) > 0) {
    cat("Samples in counts matrix but not in design file:", missing_in_design, "\n")
  }
  if (length(missing_in_counts) > 0) {
    cat("Samples in design file but not in counts matrix:", missing_in_counts, "\n")
  }
  
  # Stop execution if there are mismatches
  stop("Sample names in the counts matrix and design file do not match!")
}

# Reorder the design file to match the counts matrix
if (!all(colnames(counts) == rownames(design))) {
  design <- design[colnames(counts), ]
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = design,
                              design = ~ group) # Replace 'Condition' with the variable of interest

# Prefiltering to remove rows with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Sort results by adjusted p-value
res <- res[order(res$padj), ]

# Write results to a CSV file
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# Identify significant genes (adjusted p-value < 0.05)
significant_genes <- subset(res, padj < 0.05)

# Write significant genes to a CSV file
write.csv(as.data.frame(significant_genes), file = "significant_genes.csv")

# Generate MA plot
plotMA(dds, ylim = c(-2, 2))

# Generate a normalized counts matrix
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = "normalized_counts.csv")

# Merge results with normalized counts
res_with_counts <- as.data.frame(res)
res_with_counts$gene <- rownames(res_with_counts)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$gene <- rownames(normalized_counts_df)

# Merge by gene
merged_results <- merge(res_with_counts, normalized_counts_df, by = "gene")

# Write merged results to a CSV file
write.csv(merged_results, file = "merged_results_counts.csv")

# Perform PCA analysis
rld <- rlog(dds) # Regularized log transformation
pca_data <- plotPCA(rld, intgroup = "group", returnData = TRUE) # Replace "Condition" with your grouping variable
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot PCA
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA plot") +
  theme_minimal()

# Save PCA plot
pdf("pca_plot.pdf")
print(pca_plot)
dev.off()

# Generate heatmap for top significant genes
top_genes <- head(rownames(res[order(res$padj), ]), 50) # Select top 50 genes by adjusted p-value
heatmap_data <- assay(rld)[top_genes, ] # Extract rlog-normalized counts for top genes
pheatmap(heatmap_data, 
         annotation_col = design, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         color = colorRampPalette(c("green", "black", "red"))(50))

# Save heatmap
pdf("heatmap_top50.pdf")
pheatmap(heatmap_data, 
         annotation_col = design, 
         scale = "row", 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         color = colorRampPalette(c("green", "black", "red"))(50))
