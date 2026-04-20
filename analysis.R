# RNA-seq Differential Expression Analysis
# Author: Niranjan P B

library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load data
data <- read.table("counts.txt",
                   header = TRUE,
                   row.names = 1,
                   sep = "\t",
                   check.names = FALSE)

# Metadata
condition <- factor(c(rep("Control", 67), rep("Treated", 67)))
coldata <- data.frame(row.names = colnames(data), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = data,
                             colData = coldata,
                             design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds)

# Volcano prep
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Volcano plot
png("results/volcano_final.png", width = 1000, height = 800, res = 150)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.5) +
  theme_minimal()

dev.off()

# Heatmap
topgenes <- rownames(res)[1:30]
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[topgenes, ]

png("results/heatmap.png", width = 1000, height = 1000)
pheatmap(mat, scale = "row")
dev.off()