#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
count_file <- file.path(root, "work", "featurecounts", "feature.txt")
out_dir <- file.path(root, "data", "results")
fig_dir <- file.path(root, "figures", "differential-expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(count_file)) {
  stop("Missing featureCounts file: ", count_file,
       "\nGenerate it with featureCounts or place it at this path.")
}

fc <- read.delim(count_file, comment.char = "#", check.names = FALSE)
rownames(fc) <- fc$Geneid
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
counts <- as.matrix(fc[, setdiff(colnames(fc), annotation_columns), drop = FALSE])
colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))

samples <- c("C1", "C2", "C3", "T1", "T2", "T3")
missing_samples <- setdiff(samples, colnames(counts))
if (length(missing_samples)) stop("Missing samples: ", paste(missing_samples, collapse = ", "))
counts <- counts[, samples, drop = FALSE]
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

col_data <- data.frame(
  row.names = samples,
  condition = factor(c(rep("control", 3), rep("treatment", 3)),
                     levels = c("control", "treatment"))
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))
res <- res[order(res$padj), ]

result <- as.data.frame(res)
result$gene_id <- rownames(result)
norm <- as.data.frame(counts(dds, normalized = TRUE))
norm$gene_id <- rownames(norm)
result <- merge(result, norm, by = "gene_id", sort = FALSE)

write.table(result,
  file.path(out_dir, "deseq2_results_with_normalized_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

fdr <- subset(result, !is.na(padj) & padj < 0.05)
strict <- subset(result, !is.na(padj) & padj < 0.01 & abs(log2FoldChange) >= 1)
clean_id <- function(x) sub("^gene:", "", x)
writeLines(clean_id(subset(fdr, log2FoldChange > 0)$gene_id),
           file.path(out_dir, "de_genes_fdr05_up.txt"))
writeLines(clean_id(subset(fdr, log2FoldChange < 0)$gene_id),
           file.path(out_dir, "de_genes_fdr05_down.txt"))
writeLines(clean_id(subset(strict, log2FoldChange > 0)$gene_id),
           file.path(out_dir, "de_genes_strict_up.txt"))
writeLines(clean_id(subset(strict, log2FoldChange < 0)$gene_id),
           file.path(out_dir, "de_genes_strict_down.txt"))
writeLines(clean_id(result$gene_id[!is.na(result$padj)]),
           file.path(out_dir, "de_background_all_tested.txt"))

message(sprintf("FDR < 0.05: total=%d, up=%d, down=%d",
  nrow(fdr), sum(fdr$log2FoldChange > 0), sum(fdr$log2FoldChange < 0)))
message(sprintf("Strict GO set: total=%d, up=%d, down=%d",
  nrow(strict), sum(strict$log2FoldChange > 0), sum(strict$log2FoldChange < 0)))

png(file.path(fig_dir, "ma-plot.png"), width = 1100, height = 800, res = 150)
plotMA(res, main = "DESeq2 MA: treatment vs control", ylim = c(-3, 3))
abline(h = c(-1, 1), lty = 2, col = "grey70")
dev.off()

result$class <- ifelse(!is.na(result$padj) & result$padj < 0.01 &
  abs(result$log2FoldChange) >= 1,
  ifelse(result$log2FoldChange > 0, "Up", "Down"), "Not significant")
volcano <- ggplot(result, aes(log2FoldChange, -log10(padj), colour = class)) +
  geom_point(alpha = 0.7, size = 1, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey70") +
  geom_hline(yintercept = -log10(0.01), linetype = 2, colour = "grey70") +
  scale_colour_manual(values = c(Down = "#2C7BB6", `Not significant` = "grey45", Up = "#D7191C")) +
  labs(title = "DESeq2 volcano: treatment vs control", x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "volcano-plot.png"), volcano, width = 8, height = 6, dpi = 300)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
