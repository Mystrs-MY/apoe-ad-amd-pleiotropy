#!/usr/bin/env python3
"""Select unique gene-name-matched UKB-PPP assays for literature targets."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALIAS_FILE = ROOT / "config" / "literature_target_aliases.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--local-dir", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    provenance = read_tsv(args.provenance)
    inventory = read_tsv(args.inventory)
    primary_rows = [row for row in provenance if row.get("evidence_tier") == "Tier1"]
    primary_genes = sorted({
        row["gene_symbol"].strip().upper()
        for row in primary_rows if row.get("gene_symbol", "").strip()
    })
    local_genes = {
        path.stem.split("_")[0].upper()
        for directory in args.local_dir
        for path in directory.glob("*.tar")
    }
    files = [row for row in inventory if row.get("entity_type", "").endswith("FileEntity")]
    aliases: dict[str, set[str]] = {}
    if ALIAS_FILE.exists():
        for row in read_tsv(ALIAS_FILE):
            gene = row["standardized_gene_symbol"].strip().upper()
            aliases[gene] = {
                value.strip().upper() for value in row.get("alias_symbols", "").split(";")
                if value.strip().upper() not in {"", "NA"}
            }
    output: list[dict[str, str]] = []
    for gene in primary_genes:
        candidate_symbols = {gene} | aliases.get(gene, set())
        patterns = [re.compile(rf"^{re.escape(symbol)}_", re.IGNORECASE) for symbol in candidate_symbols]
        matches = [row for row in files if any(pattern.match(row.get("name", "")) for pattern in patterns)]
        if gene in local_genes and len(matches) == 1:
            status = "verified_provider_authorized_existing"
            reason = "A unique exact gene-prefix Olink assay is already available in the provider-authorized project resources."
        elif len(matches) == 1:
            status = "selected_for_download"
            reason = "A unique exact gene-prefix Olink assay is available; literature-to-assay confidence is evaluated downstream."
        elif len(matches) == 0:
            status = "assay_unavailable"
            reason = "No approved-symbol or official-alias gene-prefix Olink assay was found; missingness is not treated as zero."
        else:
            status = "mapping_unresolved"
            reason = "Multiple exact gene-prefix assays require assay-level review; none is selected automatically."
        match = matches[0] if len(matches) == 1 else {}
        output.append({
            "gene_symbol": gene,
            "synapse_id": match.get("synapse_id", ""),
            "synapse_file_name": match.get("name", ""),
            "selection_status": status,
            "selection_reason": reason,
        })

    fields = ["gene_symbol", "synapse_id", "synapse_file_name", "selection_status", "selection_reason"]
    write_tsv(args.output, output, fields)
    selected = sum(row["selection_status"] == "selected_for_download" for row in output)
    unavailable = sum(row["selection_status"] == "assay_unavailable" for row in output)
    unresolved = sum(row["selection_status"] == "mapping_unresolved" for row in output)
    existing = sum(row["selection_status"] == "verified_provider_authorized_existing" for row in output)
    print(f"Primary genes={len(output)}; selected={selected}; existing={existing}; unavailable={unavailable}; unresolved={unresolved}")


if __name__ == "__main__":
    main()
