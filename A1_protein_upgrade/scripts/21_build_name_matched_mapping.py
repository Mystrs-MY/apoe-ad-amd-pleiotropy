#!/usr/bin/env python3
"""Build auditable literature-target to same-assay UKB-PPP mappings."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "config" / "resources.yml").read_text(encoding="utf-8"))
PROVENANCE = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
ALIASES = ROOT / "config" / "literature_target_aliases.tsv"
TARGET_OUTPUT = ROOT / "tables" / "APOE_linkable_target_assay_mapping.tsv"
CROSS_PLATFORM_OUTPUT = ROOT / "tables" / "protein_cross_platform_mapping.tsv"
MANIFEST_OUTPUT = ROOT / "data_processed" / "ukbppp_local_assay_manifest.tsv"

COMPLEX_GENES = {"APOE", "TOMM40", "APOC1", "CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5"}
FORM_TERMS = re.compile(r"isoform|soluble|fragment|complex|cleaved|precursor|mature", re.IGNORECASE)
TAR_PATTERN = re.compile(
    r"^(?P<gene>[^_]+)_(?P<uniprot>[^_]+)_(?P<oid>OID\d+)_v(?P<version>\d+)_(?P<panel>.+)\.tar$",
    re.IGNORECASE,
)


def clean_set(series: pd.Series) -> set[str]:
    invalid = {"", "NA", "NOT_REPORTED", "MAPPING_UNRESOLVED", "NAN"}
    return {str(value).strip().upper() for value in series if str(value).strip().upper() not in invalid}


def join_values(series: pd.Series) -> str:
    values = sorted(clean_set(series))
    return ";".join(values) if values else "not_reported"


def build_manifest() -> pd.DataFrame:
    resource_dirs = [Path(path) for path in CONFIG["paths"]["ukb_ppp_tar_dirs"]]
    rows: list[dict[str, object]] = []
    seen_paths: set[str] = set()
    for directory in resource_dirs:
        for path in sorted(directory.glob("*.tar")):
            resolved = str(path.resolve()).lower()
            if resolved in seen_paths:
                continue
            seen_paths.add(resolved)
            match = TAR_PATTERN.match(path.name)
            if not match:
                rows.append({
                    "gene_symbol": "mapping_unresolved", "UniProt_ID": "mapping_unresolved",
                    "Olink_target_ID": "mapping_unresolved", "assay_version": "mapping_unresolved",
                    "Olink_panel": "mapping_unresolved", "tar_path": str(path),
                    "mapping_status": "filename_parse_failed", "resource_scope": "unresolved",
                })
                continue
            values = match.groupdict()
            rows.append({
                "gene_symbol": values["gene"].upper(), "UniProt_ID": values["uniprot"].upper(),
                "Olink_target_ID": values["oid"].upper(), "assay_version": f"v{values['version']}",
                "Olink_panel": values["panel"], "tar_path": str(path),
                "mapping_status": "parsed_from_source_filename",
                "resource_scope": "synapse_name_matched_panel" if "syn51365303" in resolved else "legacy_local_resource",
            })
    manifest = pd.DataFrame(rows).sort_values(["gene_symbol", "Olink_target_ID", "tar_path"])
    manifest.to_csv(MANIFEST_OUTPUT, sep="\t", index=False)
    return manifest


def main() -> None:
    provenance = pd.read_csv(PROVENANCE, sep="\t", dtype=str).fillna("NA")
    primary = provenance[provenance["inclusion_status"].eq("primary_high_confidence_panel")].copy()
    key = ["gene_symbol", "protein_name", "protein_form_or_isoform"]
    targets = primary[key].drop_duplicates().sort_values(key).reset_index(drop=True)
    manifest = build_manifest()
    assays = manifest[manifest["mapping_status"].eq("parsed_from_source_filename")].copy()
    alias_table = pd.read_csv(ALIASES, sep="\t", dtype=str).fillna("NA") if ALIASES.exists() else pd.DataFrame()
    alias_by_gene = alias_table.set_index("standardized_gene_symbol").to_dict("index") if len(alias_table) else {}

    rows: list[dict[str, object]] = []
    for index, target in targets.iterrows():
        mask = pd.Series(True, index=primary.index)
        for column in key:
            mask &= primary[column].eq(target[column])
        source = primary.loc[mask]
        gene = str(target["gene_symbol"]).upper()
        aliases = alias_by_gene.get(gene, {})
        alias_symbols = {
            item.strip().upper() for item in str(aliases.get("alias_symbols", "NA")).split(";")
            if item.strip().upper() not in {"", "NA"}
        }
        candidates = assays[assays["gene_symbol"].eq(gene)].copy()
        synonym_candidates = assays[assays["gene_symbol"].isin(alias_symbols)].copy()
        candidates = pd.concat([candidates, synonym_candidates], ignore_index=True).drop_duplicates("Olink_target_ID")
        candidate_count = len(candidates)
        selected = candidates.iloc[0] if candidate_count == 1 else None

        source_oids = clean_set(source["assay_target_ID"])
        source_uniprots = clean_set(source["UniProt_ID"])
        source_platforms = join_values(source["proteomic_platform"])
        target_name = str(target["protein_name"]).strip()
        form = str(target["protein_form_or_isoform"]).strip()
        canonical_name = target_name.upper() == gene
        form_ambiguous = bool(FORM_TERMS.search(f"{target_name} {form}"))
        exact_oid = bool(selected is not None and str(selected["Olink_target_ID"]).upper() in source_oids)
        exact_uniprot = bool(selected is not None and str(selected["UniProt_ID"]).upper() in source_uniprots)
        gene_match = bool(selected is not None and str(selected["gene_symbol"]).upper() == gene)
        synonym_match = bool(selected is not None and str(selected["gene_symbol"]).upper() in alias_symbols)

        if candidate_count == 0:
            confidence = "not_mapped"
            ambiguity = "assay_unavailable_after_standard_symbol_and_alias_audit"
            eligible_primary = eligible_strict = False
            reason = "No approved-symbol, official-alias, or UniProt match was found in the verified UKB-PPP inventory."
            exclusion = "UKB-PPP_assay_unavailable_after_alias_audit"
        elif candidate_count > 1:
            confidence = "low_ambiguous"
            ambiguity = "multiple_candidate_UKB_assays"
            eligible_primary = eligible_strict = False
            reason = "Multiple plausible UKB-PPP assays were found; no assay was selected automatically."
            exclusion = "one_to_many_assay_mapping_unresolved"
        elif form_ambiguous and not (exact_oid or exact_uniprot):
            confidence = "low_ambiguous"
            ambiguity = "literature_protein_form_not_resolved_by_canonical_Olink_assay"
            eligible_primary = eligible_strict = False
            reason = "The literature entry is form- or isoform-specific, whereas the unique Olink assay is canonical and not form-specific."
            exclusion = "protein_form_mapping_ambiguous"
        elif exact_oid or (exact_uniprot and canonical_name):
            confidence = "high"
            ambiguity = "none"
            eligible_primary = eligible_strict = True
            reason = "Unique gene-name match with exact literature OID or exact UniProt and an unambiguous target name."
            exclusion = "NA"
        elif gene_match and canonical_name:
            confidence = "moderate"
            ambiguity = "originating_literature_assay_annotation_incomplete_or_cross_platform"
            eligible_primary = True
            eligible_strict = False
            reason = "Unique standardized gene/protein-name match; alpha and beta will both be re-estimated from this same UKB-PPP assay."
            exclusion = "NA"
        elif synonym_match and not form_ambiguous:
            confidence = "moderate"
            ambiguity = "official_synonym_match"
            eligible_primary = True
            eligible_strict = False
            reason = "A unique official HGNC synonym matched the UKB-PPP assay; same-assay alpha and beta are required."
            exclusion = "NA"
        else:
            confidence = "low_ambiguous"
            ambiguity = "target_name_or_gene_symbol_inconsistency"
            eligible_primary = eligible_strict = False
            reason = "A unique assay exists, but the literature target name cannot be mapped without material ambiguity."
            exclusion = "target_name_mapping_ambiguous"

        rows.append({
            "literature_target_entry": f"L2T{index + 1:03d}",
            "literature_protein_name": target_name,
            "standardized_gene_symbol": gene,
            "literature_platform": source_platforms,
            "literature_assay_ID": join_values(source["assay_target_ID"]),
            "literature_UniProt": join_values(source["UniProt_ID"]),
            "literature_protein_form": form,
            "UKB_PPP_target_name": selected["gene_symbol"] if selected is not None else "mapping_unresolved",
            "UKB_PPP_gene_symbol": gene if selected is not None else "mapping_unresolved",
            "UKB_PPP_OID": selected["Olink_target_ID"] if selected is not None else "mapping_unresolved",
            "UKB_PPP_UniProt": selected["UniProt_ID"] if selected is not None else "mapping_unresolved",
            "number_of_candidate_UKB_assays": candidate_count,
            "selected_assay": selected["Olink_target_ID"] if selected is not None else "none",
            "selection_reason": reason,
            "name_match": str(canonical_name and selected is not None).lower(),
            "gene_symbol_match": str(gene_match).lower(),
            "synonym_match": str(synonym_match).lower(),
            "exact_assay_match": str(exact_oid).lower(),
            "exact_UniProt_match": str(exact_uniprot).lower(),
            "protein_form_match": "unresolved" if form_ambiguous and not (exact_oid or exact_uniprot) else "compatible",
            "mapping_confidence": confidence,
            "ambiguity_type": ambiguity,
            "high_complexity_locus": str(gene in COMPLEX_GENES or gene.startswith("HLA") or gene.startswith("IGH")).lower(),
            "eligible_for_primary": str(eligible_primary).lower(),
            "eligible_for_strict_sensitivity": str(eligible_strict).lower(),
            "exclusion_reason": exclusion,
            "manual_verification_status": "rule_based_verified_requires_manuscript_level_review" if eligible_primary else "reviewed_not_eligible",
            "notes": "Literature prioritizes the target; it does not supply beta for mediation. Same-assay UKB-PPP alpha and beta are required.",
            # Backward-compatible names consumed by existing downstream scripts.
            "target_entry_id": f"L2T{index + 1:03d}",
            "gene_symbol": gene,
            "protein_name": target_name,
            "protein_form_or_isoform": form,
            "local_Olink_target_ID": selected["Olink_target_ID"] if selected is not None else "mapping_unresolved",
            "local_Olink_UniProt_ID": selected["UniProt_ID"] if selected is not None else "mapping_unresolved",
            "mapping_status": f"name_matched_{confidence}",
            "mapping_basis": reason,
            "eligible_for_alpha_beta_triangulation": str(eligible_primary).lower(),
        })

    mapping = pd.DataFrame(rows)
    mapping.to_csv(TARGET_OUTPUT, sep="\t", index=False)
    cross = mapping[[
        "literature_target_entry", "standardized_gene_symbol", "literature_protein_name",
        "literature_platform", "literature_assay_ID", "literature_UniProt", "literature_protein_form",
        "UKB_PPP_target_name", "UKB_PPP_gene_symbol", "UKB_PPP_OID", "UKB_PPP_UniProt",
        "number_of_candidate_UKB_assays", "mapping_confidence", "ambiguity_type",
        "eligible_for_primary", "eligible_for_strict_sensitivity", "selection_reason",
        "manual_verification_status", "notes",
    ]].copy()
    cross["merge_allowed"] = cross["eligible_for_primary"].map({"true": "same_assay_reestimation_only", "false": "false"})
    cross.to_csv(CROSS_PLATFORM_OUTPUT, sep="\t", index=False)

    gene_summary = mapping.groupby("standardized_gene_symbol", as_index=False).agg(
        candidate_assays=("number_of_candidate_UKB_assays", "max"),
        primary_eligible=("eligible_for_primary", lambda x: any(value == "true" for value in x)),
        strict_eligible=("eligible_for_strict_sensitivity", lambda x: any(value == "true" for value in x)),
    )
    assert len(mapping) == 41, f"Expected 41 target entries, found {len(mapping)}"
    assert len(gene_summary) == 33, f"Expected 33 genes, found {len(gene_summary)}"
    print(f"Target entries={len(mapping)}; genes={len(gene_summary)}")
    print(f"Primary name-matched genes={gene_summary.primary_eligible.sum()}; strict genes={gene_summary.strict_eligible.sum()}")
    print(f"Output: {TARGET_OUTPUT}")


if __name__ == "__main__":
    main()
