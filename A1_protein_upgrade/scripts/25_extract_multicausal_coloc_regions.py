#!/usr/bin/env python3
"""Stream regional pQTL/GWAS records needed for multi-causal coloc."""

from __future__ import annotations

import gzip
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ROOT.parent
RESOURCE = Path(os.environ.get("A1_RESOURCE_ROOT", PROJECT_ROOT / "data" / "external"))
OUT = ROOT / "data_processed" / "multicausal_coloc_regions"

GENES = {
    "TREM2": (6, 40126786, 42126786),
    "IL6R": (1, 153377669, 155377669),
    "CFB": (6, 30913724, 32913724),
    "CFH": (1, 195621008, 197621008),
    "CFP": (23, 46413002, 48413002),
    "PCSK9": (1, 54505221, 56505221),
}

OUTCOMES = {
    "AD": "AD_Wightman_cleaned_hg19.tsv.gz",
    "Dry": "AMD_Dry_R12_cleaned_hg19.tsv.gz",
    "Wet": "AMD_Wet_R12_cleaned_hg19.tsv.gz",
    "Any": "AMD_H7_R12_cleaned_hg19.tsv.gz",
}


def extract_pqtl(gene: str, region: tuple[int, int, int]) -> int:
    chrom, start, end = region
    source = RESOURCE / "ukbppp_proteins" / "ukbppp_merged" / f"{gene}_merged.txt.gz"
    target = OUT / f"{gene}_pqtl_region.tsv.gz"
    count = 0
    with gzip.open(source, "rt", encoding="utf-8", errors="replace") as src, gzip.open(
        target, "wt", encoding="utf-8", newline=""
    ) as dst:
        header = src.readline()
        dst.write(header)
        for line in src:
            fields = line.split(" ", 4)
            if len(fields) < 3:
                continue
            try:
                variant_parts = fields[2].split(":")
                chromosome_text = variant_parts[0].upper().removeprefix("CHR")
                row_chr = 23 if chromosome_text == "X" else int(chromosome_text)
                row_pos = int(variant_parts[1])
            except ValueError:
                continue
            if row_chr == chrom and start <= row_pos <= end:
                dst.write(line)
                count += 1
    return count


def extract_gwas(outcome: str, filename: str) -> dict[str, int]:
    source = RESOURCE / "GWAS" / filename
    handles: dict[str, object] = {}
    counts = {gene: 0 for gene in GENES}
    try:
        with gzip.open(source, "rt", encoding="utf-8", errors="replace") as src:
            header = src.readline()
            for gene in GENES:
                handle = gzip.open(
                    OUT / f"{gene}_{outcome}_gwas_region.tsv.gz",
                    "wt",
                    encoding="utf-8",
                    newline="",
                )
                handle.write(header)
                handles[gene] = handle
            for line in src:
                fields = line.split("\t", 4)
                if len(fields) < 3:
                    continue
                try:
                    chromosome_text = fields[1].upper().removeprefix("CHR")
                    row_chr = 23 if chromosome_text == "X" else int(chromosome_text)
                    row_pos = int(fields[2])
                except ValueError:
                    continue
                for gene, (chrom, start, end) in GENES.items():
                    normalized_chr = 23 if row_chr == 23 else row_chr
                    if normalized_chr == chrom and start <= row_pos <= end:
                        handles[gene].write(line)
                        counts[gene] += 1
    finally:
        for handle in handles.values():
            handle.close()
    return counts


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for gene, region in GENES.items():
        rows.append(("pQTL", gene, "protein", extract_pqtl(gene, region)))
    for outcome, filename in OUTCOMES.items():
        counts = extract_gwas(outcome, filename)
        rows.extend(("GWAS", gene, outcome, count) for gene, count in counts.items())
    manifest = OUT / "regional_extraction_counts.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        handle.write("source_type\tgene\ttrait\tn_records\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    if any(row[3] < 50 and row[1] != "CFP" for row in rows):
        low = [row for row in rows if row[3] < 50 and row[1] != "CFP"]
        raise RuntimeError(f"Insufficient regional records: {low}")
    print(f"Wrote {len(rows)} regional extracts to {OUT}")


if __name__ == "__main__":
    main()
