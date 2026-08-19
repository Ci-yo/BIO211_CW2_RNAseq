#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
in_dir <- file.path(root, "data", "results")
fig_dir <- file.path(root, "figures", "differential-expression")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_ids <- function(name) scan(file.path(in_dir, name), what = "", quiet = TRUE)
up <- read_ids("de_genes_strict_up.txt")
down <- read_ids("de_genes_strict_down.txt")
background <- read_ids("de_background_all_tested.txt")

run_go <- function(genes, ontology) {
  enrichGO(
    gene = genes,
    universe = background,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = FALSE
  )
}

for (direction in c("up", "down")) {
  genes <- if (direction == "up") up else down
  for (ontology in c("BP", "MF", "CC")) {
    result <- run_go(genes, ontology)
    table <- as.data.frame(result)
    stem <- paste0("go_", direction, "_", tolower(ontology))
    write.table(table, file.path(in_dir, paste0(stem, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE)
    if (nrow(table)) {
      plot <- dotplot(result, showCategory = 10) +
        ggtitle(sprintf("GO %s (%s-regulated)", ontology, direction))
      ggsave(file.path(fig_dir, paste0(stem, ".png")), plot,
        width = 9, height = 7, dpi = 300)
    }
  }
}
