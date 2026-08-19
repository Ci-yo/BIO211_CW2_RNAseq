#!/usr/bin/env bash
set -euo pipefail

SAMPLE="${1:-C1}"
ROOT="${2:-$(pwd)}"
THREADS="${THREADS:-8}"
WORK="$ROOT/work"
REF="$ROOT/reference"
OUT="$WORK/assembly"

mkdir -p "$OUT" "$WORK/tdf"
test -s "$WORK/bam/${SAMPLE}.bam"
test -s "$REF/TAIR10.gtf"

stringtie "$WORK/bam/${SAMPLE}.bam" \
  -G "$REF/TAIR10.gtf" \
  -o "$OUT/${SAMPLE}.stringtie.gtf" \
  -A "$OUT/${SAMPLE}.gene_abundance.tsv" \
  -p "$THREADS"

# Compare reconstructed transcripts with the reference. Inspect class codes before
# calling a model novel; class code "u" is intergenic, while other codes can mark
# alternative or partially matching structures.
gffcompare -r "$REF/TAIR10.gtf" -o "$OUT/${SAMPLE}.gffcmp" \
  "$OUT/${SAMPLE}.stringtie.gtf"

# Optional IGV coverage tracks for all available BAM files.
for bam in "$WORK"/bam/*.bam; do
  sample="$(basename "$bam" .bam)"
  igvtools count -z 5 -w 10 -e 0 "$bam" "$WORK/tdf/${sample}.tdf" \
    "$REF/TAIR10.chrom.sizes"
done
