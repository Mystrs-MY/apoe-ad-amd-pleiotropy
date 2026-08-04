#!/usr/bin/env python3
"""Build the assay-level deCODE feasibility gate from the streamed supplement audit."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit" / "decode_2023_supplement_stream_audit.json"
SOURCE = ROOT / "data_raw" / "decode_public" / "41586_2023_6563_MOESM3_ESM.xlsx"
OUTPUT = ROOT / "tables" / "decode_same_platform_feasibility_gate.tsv"
FULL_ALPHA = ROOT / "tables" / "APOE_variant_to_decode_somascan_alpha.tsv"
FULL_BETA = ROOT / "tables" / "decode_smp_beta_results.tsv"
MEDIATION_SUMMARY = ROOT / "tables" / "decode_smp_mediation_summary.tsv"
TARGETS = ("CLN5", "COL10A1", "PLOD2", "SDF2", "TMEM106B", "VTN", "BRD2", "IL20RB")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


data = json.loads(AUDIT.read_text(encoding="utf-8"))
full_alpha: dict[tuple[str, str], dict[str, str]] = {}
if FULL_ALPHA.exists():
    with FULL_ALPHA.open("r", encoding="utf-8-sig", newline="") as handle:
        for alpha in csv.DictReader(handle, delimiter="\t"):
            full_alpha[(alpha["assay_target_ID"], alpha["variant"])] = alpha
beta_counts: dict[str, int] = {}
if FULL_BETA.exists():
    with FULL_BETA.open("r", encoding="utf-8-sig", newline="") as handle:
        for beta in csv.DictReader(handle, delimiter="\t"):
            if beta["method_role"] == "primary" and beta["beta_status"] == "reestimated":
                beta_counts[beta["assay_target_ID"]] = beta_counts.get(beta["assay_target_ID"], 0) + 1
assays: list[dict[str, str]] = []
for match in data["matches"]:
    if match["sheet_name"] != "ST2_somascan_assays":
        continue
    cells = match["cells"]
    gene = cells.get(f"E{match['row_number']}", "")
    if gene not in TARGETS:
        continue
    seqid_display = cells.get(f"A{match['row_number']}", "")
    seqid = seqid_display.removeprefix("SeqId.").replace("-", "_")
    assays.append(
        {
            "gene_symbol": gene,
            "SomaScan_SeqId": seqid,
            "SomaScan_display_ID": seqid_display,
            "UniProt_ID": cells.get(f"B{match['row_number']}", "NA"),
            "assay_target_name": cells.get(f"C{match['row_number']}", "NA"),
            "assay_target_full_name": cells.get(f"D{match['row_number']}", "NA"),
            "mapping_source_sheet": match["sheet_name"],
            "mapping_source_row": str(match["row_number"]),
            "mapping_confidence": "exact_aptamer_UniProt_gene",
        }
    )

alpha_rows: dict[tuple[str, str, str], dict[str, str]] = {}
for match in data["matches"]:
    if match["sheet_name"] != "ST19_somascan_pQTLs_smpnorm":
        continue
    row_number = match["row_number"]
    cells = match["cells"]
    gene = cells.get(f"C{row_number}", "")
    variant = cells.get(f"J{row_number}", "")
    seqid = cells.get(f"AT{row_number}", "")
    if gene in TARGETS and variant in {"rs429358", "rs7412"}:
        alpha_rows[(gene, seqid, variant)] = {
            "effect_allele": cells.get(f"Q{row_number}", "NA"),
            "other_allele": cells.get(f"R{row_number}", "NA"),
            "beta_adj": cells.get(f"Y{row_number}", "NA"),
            "minus_log10_p_adj": cells.get(f"AA{row_number}", "NA"),
            "joint_beta": cells.get(f"AV{row_number}", "NA"),
            "joint_p": cells.get(f"AW{row_number}", "NA"),
            "source_row": str(row_number),
        }

rows: list[dict[str, str]] = []
for assay in sorted(assays, key=lambda row: (row["gene_symbol"], row["SomaScan_SeqId"])):
    row = dict(assay)
    for variant in ("rs429358", "rs7412"):
        alpha = alpha_rows.get((row["gene_symbol"], row["SomaScan_SeqId"], variant))
        prefix = variant
        if alpha:
            row[f"{prefix}_alpha_in_significant_pQTL_supplement"] = "true"
            row[f"{prefix}_effect_allele"] = alpha["effect_allele"]
            row[f"{prefix}_other_allele"] = alpha["other_allele"]
            row[f"{prefix}_beta_adj"] = alpha["beta_adj"]
            row[f"{prefix}_minus_log10_p_adj"] = alpha["minus_log10_p_adj"]
            row[f"{prefix}_joint_beta"] = alpha["joint_beta"]
            row[f"{prefix}_joint_p"] = alpha["joint_p"]
            row[f"{prefix}_source_row"] = alpha["source_row"]
        else:
            row[f"{prefix}_alpha_in_significant_pQTL_supplement"] = "false"
            row[f"{prefix}_effect_allele"] = "NA"
            row[f"{prefix}_other_allele"] = "NA"
            row[f"{prefix}_beta_adj"] = "NA"
            row[f"{prefix}_minus_log10_p_adj"] = "NA"
            row[f"{prefix}_joint_beta"] = "NA"
            row[f"{prefix}_joint_p"] = "NA"
            row[f"{prefix}_source_row"] = "NA"
        direct = full_alpha.get((row["SomaScan_SeqId"], variant))
        row[f"{prefix}_direct_alpha_in_full_GWAS"] = "true" if direct else "false"
        row[f"{prefix}_full_GWAS_effect_allele"] = direct["harmonized_effect_allele"] if direct else "NA"
        row[f"{prefix}_full_GWAS_other_allele"] = direct["harmonized_other_allele"] if direct else "NA"
        row[f"{prefix}_full_GWAS_alpha"] = direct["alpha"] if direct else "NA"
        row[f"{prefix}_full_GWAS_SE"] = direct["alpha_SE"] if direct else "NA"
        row[f"{prefix}_full_GWAS_P"] = direct["alpha_P"] if direct else "NA"
        row[f"{prefix}_full_GWAS_N"] = direct["alpha_N"] if direct else "NA"
    archive_complete = all((row["SomaScan_SeqId"], variant) in full_alpha for variant in ("rs429358", "rs7412"))
    beta_complete = beta_counts.get(row["SomaScan_SeqId"], 0) == 4
    row.update(
        {
            "complete_aptamer_GWAS_available_in_authorized_archive": str(archive_complete).lower(),
            "same_platform_beta_reestimation_ready": str(beta_complete).lower(),
            "full_summary_data_requirement": "satisfied_provider_authorized_SMP_full_GWAS" if archive_complete else "pending",
            "gate_status": "same_platform_alpha_beta_mediation_complete" if archive_complete and beta_complete else "incomplete",
            "reestimated_outcome_beta_count": str(beta_counts.get(row["SomaScan_SeqId"], 0)),
            "absence_interpretation": "full_GWAS_direct_alpha_extracted_without_significance_selection",
            "source_file": SOURCE.name,
            "source_file_sha256": sha256(SOURCE),
            "full_alpha_table": str(FULL_ALPHA),
            "full_beta_table": str(FULL_BETA),
            "mediation_summary_table": str(MEDIATION_SUMMARY),
        }
    )
    rows.append(row)

if len(rows) != 9:
    raise AssertionError(f"Expected 9 assay rows for 8 genes (two BRD2 aptamers), found {len(rows)}")
if sum(row["rs7412_alpha_in_significant_pQTL_supplement"] == "true" for row in rows) != 1:
    raise AssertionError("Expected one directly reported rs7412 target-assay alpha")
if any(row["rs429358_alpha_in_significant_pQTL_supplement"] == "true" for row in rows):
    raise AssertionError("Unexpected rs429358 target-assay alpha in significant-pQTL supplement")

fields = list(rows[0])
with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"deCODE gate rows: {len(rows)}")
print(f"Output: {OUTPUT}")
