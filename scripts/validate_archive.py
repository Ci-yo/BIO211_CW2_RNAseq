#!/usr/bin/env python3
"""Validate the compact results shipped with this coursework archive."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def as_float(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


with (DATA / "processed" / "alignment_rates.tsv").open(
    encoding="utf-8", newline=""
) as handle:
    alignment = list(csv.DictReader(handle, delimiter="\t"))

assert [row["sample"] for row in alignment] == ["C1", "C2", "C3", "T1", "T2", "T3"]
mapped = [float(row["mapped_total_%"]) for row in alignment]
assert min(mapped) == 98.65 and max(mapped) == 99.00

up_fdr = down_fdr = up_strict = down_strict = 0
results_path = DATA / "results" / "deseq2_results_with_normalized_counts.tsv"
with results_path.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        padj = as_float(row["padj"])
        fold_change = as_float(row["log2FoldChange"])
        if math.isnan(padj) or math.isnan(fold_change):
            continue
        if padj < 0.05:
            up_fdr += fold_change > 0
            down_fdr += fold_change < 0
        if padj < 0.01 and abs(fold_change) >= 1:
            up_strict += fold_change > 0
            down_strict += fold_change < 0

assert (up_fdr, down_fdr) == (227, 226)
assert (up_strict, down_strict) == (137, 52)

expected_lists = {
    "de_genes_fdr05_up.txt": 227,
    "de_genes_fdr05_down.txt": 226,
    "de_genes_strict_up.txt": 137,
    "de_genes_strict_down.txt": 52,
}
for filename, expected_count in expected_lists.items():
    lines = (DATA / "results" / filename).read_text(encoding="utf-8").splitlines()
    assert len(lines) == expected_count, (filename, len(lines), expected_count)

print("Archive checks passed")
print(f"  samples: {len(alignment)}")
print(f"  total mapping range: {min(mapped):.2f}-{max(mapped):.2f}%")
print(f"  FDR < 0.05: {up_fdr + down_fdr} ({up_fdr} up, {down_fdr} down)")
print(f"  strict GO set: {up_strict + down_strict} ({up_strict} up, {down_strict} down)")
