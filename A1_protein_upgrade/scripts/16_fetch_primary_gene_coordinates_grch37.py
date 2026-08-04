#!/usr/bin/env python3
"""Fetch auditable GRCh37 coordinates for primary-panel genes from Ensembl REST."""

from __future__ import annotations

import csv
import time
from datetime import date
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
OUTPUT = ROOT / "config" / "protein_gene_coordinates_grch37.tsv"
BASE = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens"


def main() -> None:
    with PROVENANCE.open(encoding="utf-8", newline="") as handle:
        evidence = list(csv.DictReader(handle, delimiter="\t"))
    genes = sorted({
        row["gene_symbol"].strip().upper()
        for row in evidence
        if row.get("evidence_tier") == "Tier1" and row.get("gene_symbol", "").strip()
    })
    rows = []
    for gene in genes:
        url = f"{BASE}/{gene}?content-type=application/json"
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
        if response.status_code == 200:
            record = response.json()
            rows.append({
                "gene_symbol": gene,
                "Ensembl_gene_ID": record.get("id", "mapping_unresolved"),
                "chromosome": record.get("seq_region_name", "mapping_unresolved"),
                "gene_start_grch37": record.get("start", "mapping_unresolved"),
                "gene_end_grch37": record.get("end", "mapping_unresolved"),
                "strand": record.get("strand", "mapping_unresolved"),
                "cis_window_bp": 1_000_000,
                "source_url": url,
                "verification_date": date.today().isoformat(),
                "mapping_status": "verified_ensembl_grch37",
            })
        else:
            rows.append({
                "gene_symbol": gene,
                "Ensembl_gene_ID": "mapping_unresolved",
                "chromosome": "mapping_unresolved",
                "gene_start_grch37": "mapping_unresolved",
                "gene_end_grch37": "mapping_unresolved",
                "strand": "mapping_unresolved",
                "cis_window_bp": 1_000_000,
                "source_url": url,
                "verification_date": date.today().isoformat(),
                "mapping_status": f"http_{response.status_code}_requires_manual_verification",
            })
        time.sleep(0.1)

    fields = [
        "gene_symbol", "Ensembl_gene_ID", "chromosome", "gene_start_grch37",
        "gene_end_grch37", "strand", "cis_window_bp", "source_url",
        "verification_date", "mapping_status",
    ]
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    unresolved = sum(row["mapping_status"] != "verified_ensembl_grch37" for row in rows)
    print(f"Coordinates written: {len(rows)}; unresolved: {unresolved}; {OUTPUT}")


if __name__ == "__main__":
    main()
