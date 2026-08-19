#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
THREADS="${THREADS:-8}"
RAW="$ROOT/data/raw"
REF="$ROOT/reference"
WORK="$ROOT/work"

mkdir -p "$WORK/fastqc_raw" "$WORK/trimmed" "$WORK/fastqc_trimmed" \
  "$WORK/star_index" "$WORK/bam" "$ROOT/data/processed"

samples=(C1 C2 C3 T1 T2 T3)
raw_files=()
for sample in "${samples[@]}"; do
  test -s "$RAW/${sample}.fastq" || {
    echo "Missing $RAW/${sample}.fastq" >&2
    exit 1
  }
  raw_files+=("$RAW/${sample}.fastq")
done
test -s "$REF/TAIR10.toplevel.fa"
test -s "$REF/TAIR10.gtf"

fastqc --threads "$THREADS" --outdir "$WORK/fastqc_raw" "${raw_files[@]}"

for sample in "${samples[@]}"; do
  trim_galore --quality 20 --length 20 --output_dir "$WORK/trimmed" \
    "$RAW/${sample}.fastq"
done

fastqc --threads "$THREADS" --outdir "$WORK/fastqc_trimmed" \
  "$WORK"/trimmed/*_trimmed.fq

STAR --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$WORK/star_index" \
  --genomeFastaFiles "$REF/TAIR10.toplevel.fa" \
  --sjdbGTFfile "$REF/TAIR10.gtf" \
  --genomeSAindexNbases 12

for sample in "${samples[@]}"; do
  STAR --runThreadN "$THREADS" \
    --genomeDir "$WORK/star_index" \
    --readFilesIn "$WORK/trimmed/${sample}_trimmed.fq" \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$WORK/bam/${sample}."
  mv "$WORK/bam/${sample}.Aligned.sortedByCoord.out.bam" "$WORK/bam/${sample}.bam"
  samtools index "$WORK/bam/${sample}.bam"
done

summary="$ROOT/data/processed/alignment_rates.tsv"
printf 'sample\tinput_reads\tuniq_mapped_%%\tmulti_mapped_%%\tmapped_total_%%\tunmapped_mismatch_%%\tunmapped_short_%%\tunmapped_other_%%\n' > "$summary"
for sample in "${samples[@]}"; do
  log="$WORK/bam/${sample}.Log.final.out"
  value() { awk -F '|' -v key="$1" '$1 ~ key {gsub(/[%[:space:]]/, "", $2); print $2}' "$log"; }
  input="$(value 'Number of input reads')"
  unique="$(value 'Uniquely mapped reads %')"
  multi="$(value '% of reads mapped to multiple loci')"
  mismatch="$(value '% of reads unmapped: too many mismatches')"
  short="$(value '% of reads unmapped: too short')"
  other="$(value '% of reads unmapped: other')"
  total="$(awk -v a="$unique" -v b="$multi" 'BEGIN {printf "%.2f", a+b}')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$input" "$unique" "$multi" "$total" "$mismatch" "$short" "$other" >> "$summary"
done

multiqc --outdir "$WORK/multiqc" "$WORK/fastqc_raw" "$WORK/fastqc_trimmed" "$WORK/bam"
