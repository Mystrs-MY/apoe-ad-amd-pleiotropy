#!/usr/bin/env python3
"""Apply study decisions and build conservative protein mapping/exclusion tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "literature_source_registry.tsv"
SCREENING = ROOT / "literature" / "study_screening.tsv"
AVAILABILITY = ROOT / "literature" / "full_text_availability.tsv"
PROVENANCE = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
MAPPING = ROOT / "tables" / "protein_cross_platform_mapping.tsv"
EXCLUSIONS = ROOT / "tables" / "excluded_proteins_with_reasons.tsv"
STUDY_STATUS = ROOT / "literature" / "included_study_evidence_registry.tsv"


def clean_items(values: pd.Series) -> list[str]:
    excluded = {"", "NA", "not_reported", "mapping_unresolved", "requires_source_mapping",
                "mapping_requires_source_pQTL_table", "gene-level label only"}
    return sorted({str(value).strip() for value in values if str(value).strip() not in excluded})


def join_items(values: pd.Series) -> str:
    items = clean_items(values)
    return ";".join(items) if items else "mapping_unresolved"


def apply_screening_decisions() -> None:
    registry = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("NA")
    screening = pd.read_csv(SCREENING, sep="\t", dtype=str).fillna("")
    availability = pd.read_csv(AVAILABILITY, sep="\t", dtype=str).fillna("")

    registry_ids = set(registry["record_id"])
    missing = registry_ids - set(screening["record_id"])
    if missing:
        raise ValueError(f"Registry records absent from study_screening.tsv: {sorted(missing)}")

    availability_map = availability.drop_duplicates("record_id").set_index("record_id").to_dict("index")
    registry_map = registry.set_index("record_id").to_dict("index")
    for index, row in screening.iterrows():
        record_id = row["record_id"]
        if record_id not in registry_map:
            continue
        decision = registry_map[record_id]
        available = availability_map.get(record_id, {})
        screening.at[index, "full_text_status"] = available.get("full_text_status", "not_assessed")
        screening.at[index, "full_text_decision"] = decision["full_text_decision"]
        screening.at[index, "full_text_exclusion_reason"] = (
            decision["decision_reason"] if decision["full_text_decision"] in {"exclude_primary", "pending_manual_download"}
            else ""
        )
        screening.at[index, "supplement_status"] = available.get("supplement_status", "not_assessed")
        screening.at[index, "duplicate_group"] = decision["duplicate_group"]
        screening.at[index, "evidence_independence"] = decision["evidence_independence"]
        screening.at[index, "manual_verification_status"] = decision["manual_verification_status"]
        screening.at[index, "notes"] = (
            f"{decision['evidence_role']}; {decision['extraction_status']}; {decision['decision_reason']}"
        )
    screening.to_csv(SCREENING, sep="\t", index=False)

    merged = registry.merge(
        availability[["record_id", "PMID", "DOI", "PMCID", "full_text_status", "supplement_status", "supplement_files", "errors"]],
        on="record_id", how="left",
    )
    merged.to_csv(STUDY_STATUS, sep="\t", index=False)


def build_mapping() -> None:
    evidence = pd.read_csv(PROVENANCE, sep="\t", dtype=str).fillna("NA")
    groups = []
    group_keys = ["gene_symbol", "protein_name", "protein_form_or_isoform"]
    for key, group in evidence.groupby(group_keys, dropna=False):
        gene, protein, form = key
        olink = group[group["proteomic_platform"].str.contains("Olink", case=False, na=False)]
        soma = group[group["proteomic_platform"].str.contains("Soma", case=False, na=False)]
        decode = group[group["proteomic_platform"].str.contains("deCODE", case=False, na=False)]
        all_targets = clean_items(group["assay_target_ID"])
        isoform_specific = "isoform" in str(form).lower() or "isoform" in str(protein).lower()
        one_many = len(all_targets) > 1
        complex_locus = (group["high_complexity_locus"] == "true").any()
        cross_reactivity = join_items(group["assay_cross_reactivity_risk"])
        confidence_values = set(group["mapping_confidence"])
        confidence = "high" if confidence_values == {"high"} else "moderate" if "high" in confidence_values or "moderate" in confidence_values else "low"
        if isoform_specific or one_many or complex_locus:
            merge = "false"
        else:
            merge = "conditional_only_after_assay_unit_and_direction_review"
        groups.append({
            "gene_symbol": gene,
            "protein_name": protein,
            "UniProt_ID": join_items(group["UniProt_ID"]),
            "Olink_target": join_items(olink["assay_target_ID"]) if len(olink) else "mapping_unresolved",
            "SomaScan_aptamer": join_items(soma["assay_target_ID"]) if len(soma) else "mapping_unresolved",
            "deCODE_protein_label": join_items(decode["assay_target_name"]) if len(decode) else "mapping_unresolved",
            "protein_form_or_isoform": form,
            "one_to_one_or_one_to_many_mapping": "one_to_many_in_extracted_evidence" if one_many else "one_to_one_within_extracted_evidence_only",
            "assay_cross_reactivity": cross_reactivity,
            "mapping_confidence": confidence,
            "allow_cross_platform_merge": merge,
            "high_complexity_locus": "true" if complex_locus else "false",
            "source_PMIDs": join_items(group["PMID"]),
            "source_assay_target_IDs": ";".join(all_targets) if all_targets else "mapping_unresolved",
            "manual_verification_status": "requires_manual_review_before_cross_platform_pooling",
            "notes": "Gene-symbol agreement alone does not authorize assay pooling.",
        })
    mapping = pd.DataFrame(groups).sort_values(["gene_symbol", "protein_name", "protein_form_or_isoform"])
    mapping.to_csv(MAPPING, sep="\t", index=False)


def build_exclusions() -> None:
    evidence = pd.read_csv(PROVENANCE, sep="\t", dtype=str).fillna("NA")
    excluded = evidence[evidence["inclusion_status"] != "primary_high_confidence_panel"].copy()
    columns = [
        "record_id", "PMID", "protein_name", "gene_symbol", "assay_target_ID", "proteomic_platform",
        "disease", "disease_subtype", "evidence_tier", "inclusion_status", "exclusion_reason",
        "high_complexity_locus", "mapping_confidence", "PP_H4", "replication_status",
        "source_table_file", "source_sheet", "source_row", "manual_verification_status", "notes",
    ]
    excluded[columns].drop_duplicates().sort_values(
        ["disease", "disease_subtype", "gene_symbol", "PMID", "assay_target_ID"]
    ).to_csv(EXCLUSIONS, sep="\t", index=False)


def main() -> None:
    apply_screening_decisions()
    build_mapping()
    build_exclusions()
    print(f"Updated: {SCREENING}")
    print(f"Created: {STUDY_STATUS}")
    print(f"Created: {MAPPING}")
    print(f"Created: {EXCLUSIONS}")


if __name__ == "__main__":
    main()
