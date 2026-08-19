# Data notes

## Included

- `processed/alignment_rates.tsv`: STAR mapping summary for all six libraries.
- `processed/featurecounts_summary.tsv`: assignment categories reported by featureCounts.
- `results/deseq2_results_with_normalized_counts.tsv`: DESeq2 statistics plus normalized counts.
- `results/de_genes_fdr05_*.txt`: genes significant at `padj < 0.05`, split by fold-change direction.
- `results/de_genes_strict_*.txt`: GO input lists using `padj < 0.01` and `|log2FC| >= 1`.
- `results/de_background_all_tested.txt`: tested-gene background for enrichment.
- `results/sessionInfo.txt`: recorded R environment information from the coursework run.

## Not included

Raw FASTQ, BAM, TDF, reference FASTA/GTF, STAR indices, `.RData`, repeated FastQC archives and draft outputs are excluded. These files are large, reproducible or duplicated. The original local working directory was approximately 2.31 GB; the compressed raw-data bundle alone was approximately 977 MB.

## Sample design

| Sample | Condition | Role |
|---|---|---|
| C1 | control | biological replicate 1 |
| C2 | control | biological replicate 2 |
| C3 | control | biological replicate 3 |
| T1 | iron-deficiency treatment | biological replicate 1 |
| T2 | iron-deficiency treatment | biological replicate 2 |
| T3 | iron-deficiency treatment | biological replicate 3 |

The practical files are a teaching subset linked to GSE163190. Consult the source record for full experimental metadata before reusing the data beyond this coursework context.
