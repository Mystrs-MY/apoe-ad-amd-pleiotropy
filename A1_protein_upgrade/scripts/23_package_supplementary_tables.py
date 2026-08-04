#!/usr/bin/env python3
"""Materialize submission-facing Supplementary Tables S18-S29 with provenance hashes."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tables" / "supplementary_name_match_revision"

COPIES = {
    "TableS18_Literature_Study_Screening.tsv": ROOT / "literature" / "study_screening.tsv",
    "TableS19_Protein_Provenance_Master.tsv": ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv",
    "TableS20_Cross_Platform_Mapping.tsv": ROOT / "tables" / "protein_cross_platform_mapping.tsv",
    "TableS21_Protein_to_UKB_PPP_Assay_Mapping.tsv": ROOT / "tables" / "APOE_linkable_target_assay_mapping.tsv",
    "TableS22_APOE_Variant_to_Protein_Alpha.tsv": ROOT / "tables" / "APOE_variant_to_literature_proteins_alpha.tsv",
    "TableS23_Same_Assay_Protein_Beta.tsv": ROOT / "tables" / "literature_panel_beta_results.tsv",
    "TableS24_APOE_Linkable_Eligibility_Flow.tsv": ROOT / "tables" / "APOE_linkable_subset_flow.tsv",
    "TableS25_Expanded_Primary_Two_Step_Mediation.tsv": ROOT / "tables" / "APOE_linkable_two_step_mediation.tsv",
    "TableS26_Cis_Only_Two_Step_Mediation.tsv": ROOT / "tables" / "APOE_linkable_two_step_mediation_cis_sensitivity.tsv",
    "TableS27_Literature_vs_Biology_Guided_Comparison.tsv": (
        ROOT / "tables" / "literature_vs_biology_guided_panel_comparison.tsv"
    ),
    "TableS28_Excluded_Name_Matched_Panel_Genes.tsv": ROOT / "tables" / "name_matched_panel_excluded_genes.tsv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for output_name, source in COPIES.items():
        if not source.exists():
            raise FileNotFoundError(source)
        destination = OUT / output_name
        shutil.copy2(source, destination)
        manifest_rows.append({
            "supplementary_table": output_name.split("_")[0],
            "output_file": output_name,
            "source_file": str(source.relative_to(ROOT)),
            "row_count": row_count(destination),
            "bytes": destination.stat().st_size,
            "sha256": sha256(destination),
            "status": "materialized_from_verified_machine_readable_source",
        })

    figure_source = ROOT / "figures" / "source_data"
    for source in sorted(figure_source.glob("Figure_5*.csv")):
        manifest_rows.append({
            "supplementary_table": "TableS29",
            "output_file": source.name,
            "source_file": str(source.relative_to(ROOT)),
            "row_count": row_count(source),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
            "status": "Figure_5_source_data",
        })
    for source in [ROOT / "logs" / "name_match_revision_QA.tsv"]:
        manifest_rows.append({
            "supplementary_table": "TableS29",
            "output_file": source.name,
            "source_file": str(source.relative_to(ROOT)),
            "row_count": row_count(source),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
            "status": "Figure_or_analysis_QA",
        })
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "TableS29_Figure5_Source_Data_and_QA_Manifest.tsv", sep="\t", index=False)
    assert len(manifest[manifest["supplementary_table"].eq("TableS29")]) >= 5
    print(f"Supplementary Tables S18-S29 materialized in: {OUT}")


if __name__ == "__main__":
    main()
