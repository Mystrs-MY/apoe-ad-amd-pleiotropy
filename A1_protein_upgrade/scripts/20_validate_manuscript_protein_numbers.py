#!/usr/bin/env python3
"""Check stable protein-layer results against the Chinese canonical manuscript."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VERSION = "V1.0"
MANUSCRIPT = ROOT.parent / "manuscripts" / f"Article1_完整手稿_中文_整合初稿_{VERSION}.md"
MEDIATION = ROOT / "tables" / "APOE_linkable_two_step_mediation.tsv"
OUTPUT = ROOT / "logs" / "manuscript_protein_numeric_QA.tsv"


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    main_text = text.split("# Supplementary Information", 1)[0]
    mediation = pd.read_csv(MEDIATION, sep="\t")
    totals = mediation[mediation["row_type"] == "total"].copy()
    protein_rows = mediation[mediation["row_type"] == "protein"].copy()
    eligible_proteins = int(protein_rows["gene_symbol"].nunique())
    mediation_paths = int(len(protein_rows))
    checks: list[dict[str, str]] = []
    for _, row in totals.iterrows():
        token = f"{100 * row['mediated_proportion']:.2f}%"
        checks.append({
            "check": f"{row['variant']}_{row['outcome']}_mediated_percent",
            "expected": token,
            "observed": token if token in main_text else "not_found",
            "pass": str(token in main_text).lower(),
        })
    for label, token in [
        ("primary_target_entries", "41 个"),
        ("eligible_assay_proteins", f"{eligible_proteins} 个"),
        ("mediation_paths", f"{mediation_paths} 条"),
        ("main_beta_signal", "CSF2→湿性 AMD"),
    ]:
        checks.append({
            "check": label, "expected": token,
            "observed": token if token in main_text else "not_found",
            "pass": str(token in main_text).lower(),
        })
    forbidden = {
        "old_fixed_29_threshold": r"P\s*<\s*0\.05/29",
        "unsupported_triangular_pool": r"β_pooled",
        "unsupported_cmap_result": r"CMap 分析提示",
        "old_29_protein_primary_claim": r"29 种非 APOE 蛋白合计仅解释",
    }
    for label, pattern in forbidden.items():
        found = bool(re.search(pattern, main_text, flags=re.IGNORECASE))
        checks.append({
            "check": label, "expected": "absent",
            "observed": "found" if found else "absent",
            "pass": str(not found).lower(),
        })
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "expected", "observed", "pass"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)
    if not all(row["pass"] == "true" for row in checks):
        raise SystemExit(f"Manuscript numeric QA failed; inspect {OUTPUT}")
    print(f"Manuscript protein numeric QA passed: {OUTPUT}")


if __name__ == "__main__":
    main()
