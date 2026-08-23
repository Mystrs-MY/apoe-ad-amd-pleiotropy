#!/usr/bin/env python3
"""Integrate newly retrieved priority studies without promoting incompatible evidence."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "literature" / "fulltext_supplements"
TABLE_DIR = ROOT / "tables" / "supplementary_name_match_revision"
SCIENCE = ROOT / "data_processed" / "priority_study_supplements" / "PMID42384774" / "adx4852_tables_s1_to_s30.xlsx"
PRIORITY_PMIDS = {"34381170", "40397384", "40452368", "42384774"}
EXPECTED_BASE_ROWS = 284
EXPECTED_PRIORITY_ROWS = 61
EXPECTED_TOTAL_ROWS = 345


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_screening() -> None:
    path = TABLE_DIR / "TableS18_Literature_Study_Screening.tsv"
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    updates = {
        "42384774": ("publisher_full_text_and_supplement_verified", "include_layer1_tier2",
                     "publisher_supplement_downloaded", "full_text_and_source_tables_verified",
                     "Bonferroni-significant cis-PWAS candidates extracted; effect is a PWAS Z score, not a standardized MR beta."),
        "40452368": ("publisher_full_text_and_supplement_verified", "include_layer1_tier2",
                     "publisher_supplement_downloaded", "full_text_and_source_tables_verified",
                     "deCODE and UKB-PPP protein-outcome tables and colocalization tables extracted; GRN/CR1 overlap prespecified panel targets."),
        "40397384": ("publisher_full_text_and_supplement_verified", "include_layer1_tier2",
                     "publisher_supplement_downloaded", "full_text_and_source_tables_verified",
                     "Bidirectional UKB/deCODE results extracted; no same-assay APOE-variant alpha supplied."),
        "34381170": ("publisher_full_text_and_supplement_verified", "include_layer1_tier2",
                     "publisher_supplement_downloaded", "full_text_and_source_tables_verified",
                     "Brain and blood evidence separated; only circulating-protein rows retained in the study source universe."),
    }
    for pmid, values in updates.items():
        mask = frame["PMID"].eq(pmid)
        if not mask.any():
            raise RuntimeError(f"Priority PMID missing from screening table: {pmid}")
        frame.loc[mask, "full_text_status"] = values[0]
        frame.loc[mask, "full_text_decision"] = values[1]
        frame.loc[mask, "full_text_exclusion_reason"] = ""
        frame.loc[mask, "supplement_status"] = values[2]
        frame.loc[mask, "manual_verification_status"] = values[3]
        frame.loc[mask, "notes"] = values[4]
        frame.loc[mask, "proteome_scale"] = "proteome_wide_or_large_scale_verified"
        frame.loc[mask, "genetic_causal_method"] = "verified_genetically_predicted_protein_analysis"
    frame.to_csv(path, sep="\t", index=False)
    frame.to_csv(ROOT / "literature" / "study_screening.tsv", sep="\t", index=False)

    registry_path = ROOT / "literature" / "included_study_evidence_registry.tsv"
    registry = pd.read_csv(registry_path, sep="\t", dtype=str).fillna("")
    supplement_files = {
        "42384774": "PMID42384774_sm.pdf;PMID42384774_mdar.pdf;PMID42384774_tables_s1_to_s30.zip",
        "40452368": "PMID40452368_supp_methods.docx;PMID40452368_supp3.xlsx;PMID40452368_supp4.xlsx;PMID40452368_supp5.xlsx;PMID40452368_supp6.xlsx;PMID40452368_supp7.xlsx;PMID40452368_supp8.xlsx;PMID40452368_supp9.xlsx;PMID40452368_supp10.xlsx",
        "40397384": "PMID40397384_supp_methods.docx;PMID40397384_supp2.xlsx",
        "34381170": "PMID34381170_supplement.xlsx",
    }
    for pmid, values in updates.items():
        mask = registry["PMID"].eq(pmid)
        if not mask.any():
            raise RuntimeError(f"Priority PMID missing from evidence registry: {pmid}")
        registry.loc[mask, "full_text_decision"] = "include_external_tier2"
        registry.loc[mask, "evidence_role"] = "literature_source_universe_tier2"
        registry.loc[mask, "manual_verification_status"] = values[3]
        registry.loc[mask, "extraction_status"] = "extracted_to_priority_study_target_hits"
        registry.loc[mask, "decision_reason"] = values[4]
        registry.loc[mask, "full_text_status"] = "publisher_material_cached"
        registry.loc[mask, "supplement_status"] = "publisher_supplements_cached"
        registry.loc[mask, "supplement_files"] = supplement_files[pmid]
        registry.loc[mask, "errors"] = ""
    registry.to_csv(registry_path, sep="\t", index=False)


def empty_provenance_row(columns: list[str]) -> dict[str, object]:
    return {column: "not_reported" for column in columns}


def append_provenance() -> None:
    path = TABLE_DIR / "TableS19_Protein_Provenance_Master.tsv"
    master = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    columns = master.columns.tolist()
    if "record_id" not in columns or "PMID" not in columns or "evidence_tier" not in columns:
        raise RuntimeError("Table S19 is missing required provenance fields")
    master = master[~master["PMID"].isin(PRIORITY_PMIDS)].copy()
    if len(master) != EXPECTED_BASE_ROWS:
        raise RuntimeError(f"Expected {EXPECTED_BASE_ROWS} base evidence rows, found {len(master)}")
    if not master["record_id"].is_unique:
        duplicated = master.loc[master["record_id"].duplicated(keep=False), "record_id"].unique()[:10]
        raise RuntimeError(f"Base provenance record_id values are not unique: {duplicated.tolist()}")
    new_rows: list[dict[str, object]] = []

    pwas = pd.read_excel(SCIENCE, sheet_name="ST 1. AD PWAS Results", skiprows=6)
    pwas = pwas[pwas["bonsig"].eq(True)].copy()
    crosswalk = pd.read_excel(SCIENCE, sheet_name="ST 30. UKB SomaScan  v. Olink", skiprows=7)
    crosswalk = crosswalk.drop_duplicates("seqid")
    crosswalk["aptamer_norm"] = crosswalk["seqid"].astype(str).str.replace(".", "_", regex=False)
    for _, row in pwas.iterrows():
        rec = empty_provenance_row(columns)
        gene = str(row["Gene symbol"])
        aptamer = str(row["Aptamer ID"])
        match = crosswalk[crosswalk["aptamer_norm"].eq(aptamer)]
        rec.update({
            "record_id": f"PMID42384774_{aptamer}", "protein_name": gene, "gene_symbol": gene,
            "UniProt_ID": match["uniprot"].iloc[0] if len(match) and pd.notna(match["uniprot"].iloc[0]) else "mapping_unresolved",
            "assay_target_name": match["target"].iloc[0] if len(match) and pd.notna(match["target"].iloc[0]) else gene,
            "assay_target_ID": aptamer, "proteomic_platform": "SomaScan cis-PWAS",
            "protein_form_or_isoform": "aptamer_specific_target", "disease": "AD", "disease_subtype": "AD",
            "outcome_definition": "AD dementia GWAS", "source_study_first_author": "Walker KA",
            "publication_year": "2026", "journal": "Science Translational Medicine", "PMID": "42384774",
            "DOI": "10.1126/scitranslmed.adx4852", "publication_status": "peer_reviewed",
            "pQTL_source": "published plasma cis-pQTL prediction models", "ancestry": "European American",
            "cis_trans_status": "cis_prediction_model", "MR_method": "proteome-wide association study",
            "P_value": row["PWAS.P"], "multiple_testing_method": "study Bonferroni threshold",
            "corrected_significance": "true", "effect_direction": "positive_PWAS_Z" if row["PWAS.Z"] > 0 else "negative_PWAS_Z",
            "outcome_GWAS": "published AD GWAS", "overlap_with_primary_outcome_GWAS": "likely_overlap",
            "same_outcome_reanalysis": "false", "evidence_tier": "Tier2",
            "high_complexity_locus": "true" if gene in {"APOE", "PVRL2", "APOC1", "IGHG1|IGHG2|IGHG3|IGHG4|IGK@|IGL@"} else "false",
            "assay_cross_reactivity_risk": "multi_gene_or_complex_mapping" if "|" in gene else "SomaScan_cross_platform_review_required",
            "inclusion_status": "literature_source_universe_tier2",
            "exclusion_reason": "PWAS_Z_not_standardized_beta;same_assay_alpha_beta_not_available",
            "APOE_rs429358_alpha_available": "false", "APOE_rs7412_alpha_available": "false",
            "alpha_source": "not_available_for_same_SomaScan_assay", "eligible_for_two_step_MR": "false",
            "manual_verification_status": "full_text_and_source_table_verified",
            "notes": "Bonferroni-significant PWAS candidate; retained as external evidence and not interpreted as a de novo discovery in the present analysis.",
            "source_table_file": str(SCIENCE), "source_sheet": "ST 1. AD PWAS Results",
            "source_row": int(row.name) + 8, "beta_source_type": "published_PWAS_Z_not_beta",
            "mapping_confidence": "high" if len(match) and pd.notna(match["oid"].iloc[0]) else "moderate",
            "evidence_independence": "partially_overlapping", "duplicate_group": "AD_plasma_cis_PWAS_2026",
            "allele_harmonization_status": "not_applicable_to_PWAS_Z",
            "derived_field_notes": f"PWAS.Z={row['PWAS.Z']}; no beta/SE conversion performed.",
        })
        new_rows.append(rec)

    bidirectional = SUPP / "PMID40397384_supp2.xlsx"
    for sheet, pqtl_source, outcome_name in [("SuppTable4", "UKB-PPP", "AD-IGAP"), ("SuppTable6", "deCODE SomaScan", "AD-EADB")]:
        data = pd.read_excel(bidirectional, sheet_name=sheet, skiprows=1)
        data = data[(data["Method"] == "Inverse variance weighted") & (pd.to_numeric(data["Padj-BH"], errors="coerce") < 0.05)]
        for _, row in data.iterrows():
            gene = str(row["Exposure"])
            rec = empty_provenance_row(columns)
            rec.update({
                "record_id": f"PMID40397384_{sheet}_{gene}", "protein_name": gene, "gene_symbol": gene,
                "assay_target_name": gene, "proteomic_platform": pqtl_source, "protein_form_or_isoform": "assay_specific_target",
                "disease": "AD", "disease_subtype": "AD", "outcome_definition": outcome_name,
                "source_study_first_author": "Li Y", "publication_year": "2025", "journal": "Journal of Alzheimer's Disease",
                "PMID": "40397384", "DOI": "10.1177/13872877251345151", "publication_status": "peer_reviewed",
                "pQTL_source": pqtl_source, "ancestry": "European", "cis_trans_status": "mixed_or_not_explicitly_cis_only",
                "instrument_definition": "study-selected genome-wide instruments", "MR_method": "IVW",
                "beta_or_logOR": np.log(float(row["OR"])), "SE": row["SE"], "OR": row["OR"],
                "CI_lower": row["Lower95"], "CI_upper": row["Upper95"], "P_value": row["Pval"],
                "multiple_testing_method": "Benjamini-Hochberg FDR", "corrected_significance": "true",
                "effect_direction": "higher_protein_lower_risk" if float(row["OR"]) < 1 else "higher_protein_higher_risk",
                "outcome_GWAS": outcome_name, "overlap_with_primary_outcome_GWAS": "partially_overlapping_or_related",
                "same_outcome_reanalysis": "false", "evidence_tier": "Tier2", "high_complexity_locus": "false",
                "assay_cross_reactivity_risk": "platform_mapping_required", "inclusion_status": "literature_source_universe_tier2",
                "exclusion_reason": "no_colocalization_or_independent_replication;published_beta_not_same_assay_primary_analysis_beta",
                "APOE_rs429358_alpha_available": "false", "APOE_rs7412_alpha_available": "false",
                "eligible_for_two_step_MR": "false", "manual_verification_status": "full_text_and_source_table_verified",
                "notes": str(row.get("Note", "")), "source_table_file": str(bidirectional), "source_sheet": sheet,
                "source_row": int(row.name) + 3, "beta_source_type": "published_beta",
                "mapping_confidence": "moderate", "evidence_independence": "partially_overlapping",
                "duplicate_group": f"AD_{pqtl_source}_2025", "allele_harmonization_status": "not_applicable_to_published_beta",
            })
            new_rows.append(rec)

    replicated = SUPP / "PMID40452368_supp5.xlsx"
    data = pd.read_excel(replicated, skiprows=2)
    for _, row in data.iterrows():
        gene = str(row["Protein"])
        source = str(row["Source"])
        lower, upper = str(row["95% CI"]).split("-", 1)
        rec = empty_provenance_row(columns)
        rec.update({
            "record_id": f"PMID40452368_{source}_{gene}", "protein_name": gene, "gene_symbol": gene,
            "assay_target_name": gene, "proteomic_platform": "deCODE SomaScan" if source == "deCODE" else "UKB-PPP Olink",
            "protein_form_or_isoform": "assay_specific_target", "disease": "AD", "disease_subtype": "AD",
            "outcome_definition": "AD GWAS", "source_study_first_author": "Yu K", "publication_year": "2025",
            "journal": "Journal of Alzheimer's Disease", "PMID": "40452368", "DOI": "10.1177/13872877251344572",
            "publication_status": "peer_reviewed", "pQTL_source": source, "ancestry": "European",
            "cis_trans_status": "mixed_or_not_explicitly_cis_only", "MR_method": "study-reported MR",
            "beta_or_logOR": np.log(float(row["OR"])), "OR": row["OR"], "CI_lower": lower, "CI_upper": upper,
            "P_value": row["p"], "multiple_testing_method": "study Bonferroni threshold", "corrected_significance": "true",
            "effect_direction": "higher_protein_lower_risk" if float(row["OR"]) < 1 else "higher_protein_higher_risk",
            "outcome_GWAS": "published AD GWAS", "overlap_with_primary_outcome_GWAS": "partially_overlapping_or_related",
            "same_outcome_reanalysis": "false", "evidence_tier": "Tier2", "high_complexity_locus": "true" if gene in {"BIN1", "CR1"} else "false",
            "assay_cross_reactivity_risk": "platform_mapping_required", "inclusion_status": "literature_source_universe_tier2",
            "exclusion_reason": "published_cross_platform_beta_not_same_assay_primary_analysis_beta;highest_tier_requires_colocalization_or_replication_review",
            "APOE_rs429358_alpha_available": "false", "APOE_rs7412_alpha_available": "false", "eligible_for_two_step_MR": "false",
            "manual_verification_status": "full_text_and_source_table_verified",
            "notes": "Study-Bonferroni-significant protein-outcome association; same dataset reports GRN/CR1 colocalization support.",
            "source_table_file": str(replicated), "source_sheet": "Supplemental Table 3", "source_row": int(row.name) + 4,
            "beta_source_type": "published_beta", "mapping_confidence": "moderate", "evidence_independence": "partially_overlapping",
            "duplicate_group": f"AD_{source}_2025", "allele_harmonization_status": "not_applicable_to_published_beta",
        })
        new_rows.append(rec)

    ou = SUPP / "PMID34381170_supplement.xlsx"
    for sheet, outcome_label in [("Table S14", "AD_GWAS_1"), ("Table S15", "AD_GWAS_2")]:
        raw = pd.read_excel(ou, sheet_name=sheet, header=None)
        ace = raw[raw.apply(lambda row: any(str(value).strip().upper() == "ACE" for value in row), axis=1)]
        for row_index, row in ace.iterrows():
            values = row.dropna().tolist()
            rec = empty_provenance_row(columns)
            rec.update({
                "record_id": f"PMID34381170_{sheet.replace(' ', '_')}_ACE", "protein_name": "ACE", "gene_symbol": "ACE",
                "assay_target_name": "ACE", "proteomic_platform": "AGES blood proteome",
                "protein_form_or_isoform": "assay_specific_target", "disease": "AD", "disease_subtype": "AD",
                "outcome_definition": outcome_label, "source_study_first_author": "Ou YN", "publication_year": "2021",
                "journal": "Molecular Psychiatry", "PMID": "34381170", "DOI": "10.1038/s41380-021-01251-6",
                "publication_status": "peer_reviewed", "pQTL_source": "AGES blood proteome", "ancestry": "European",
                "cis_trans_status": "single_instrument_as_reported", "instrument_definition": "reported Wald-ratio instrument",
                "MR_method": str(values[2]), "beta_or_logOR": values[4], "SE": values[5], "P_value": values[6],
                "multiple_testing_method": "study-wide threshold", "corrected_significance": "true",
                "effect_direction": "higher_protein_higher_risk" if float(values[4]) > 0 else "higher_protein_lower_risk",
                "outcome_GWAS": outcome_label, "overlap_with_primary_outcome_GWAS": "possible", "same_outcome_reanalysis": "false",
                "evidence_tier": "Tier2", "high_complexity_locus": "false", "assay_cross_reactivity_risk": "cross_platform_mapping_required",
                "inclusion_status": "literature_source_universe_tier2",
                "exclusion_reason": "single_instrument_published_beta;same_assay_APOE_alpha_unavailable",
                "APOE_rs429358_alpha_available": "false", "APOE_rs7412_alpha_available": "false", "eligible_for_two_step_MR": "false",
                "manual_verification_status": "full_text_and_source_table_verified",
                "notes": "Circulating blood-protein row retained separately from brain-protein evidence.",
                "source_table_file": str(ou), "source_sheet": sheet, "source_row": int(row_index) + 1,
                "beta_source_type": "published_beta", "mapping_confidence": "moderate",
                "evidence_independence": "partially_overlapping", "duplicate_group": "AD_AGES_blood_2021",
                "allele_harmonization_status": "not_applicable_to_published_beta",
            })
            new_rows.append(rec)

    priority = pd.DataFrame(new_rows, columns=columns)
    if len(priority) != EXPECTED_PRIORITY_ROWS:
        raise RuntimeError(f"Expected {EXPECTED_PRIORITY_ROWS} priority-study rows, found {len(priority)}")
    if set(priority["PMID"]) != PRIORITY_PMIDS:
        raise RuntimeError(f"Unexpected priority-study PMIDs: {sorted(set(priority['PMID']))}")
    if not priority["record_id"].is_unique:
        duplicated = priority.loc[priority["record_id"].duplicated(keep=False), "record_id"].unique()[:10]
        raise RuntimeError(f"Priority-study record_id values are not unique: {duplicated.tolist()}")

    combined = pd.concat([master, priority], ignore_index=True)
    if len(combined) != EXPECTED_TOTAL_ROWS or not combined["record_id"].is_unique:
        raise RuntimeError("Integrated provenance failed the 345-row unique-record assertion")
    tier_counts = combined["evidence_tier"].value_counts().to_dict()
    if tier_counts.get("Tier1", 0) != 52 or tier_counts.get("Tier2", 0) != 293:
        raise RuntimeError(f"Unexpected evidence-tier counts: {tier_counts}")
    combined = combined.sort_values(
        ["disease", "disease_subtype", "gene_symbol", "PMID", "assay_target_ID", "record_id"]
    ).reset_index(drop=True)
    combined.to_csv(path, sep="\t", index=False)
    combined.to_csv(ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv", sep="\t", index=False)


def write_crosswalk_and_boundaries() -> None:
    crosswalk = pd.read_excel(SCIENCE, sheet_name="ST 30. UKB SomaScan  v. Olink", skiprows=7)
    genes = ["TREM2", "ACE", "BCAM", "CD55", "LILRB1", "LILRB5", "SCARA5"]
    out = crosswalk[crosswalk["gene_name"].isin(genes)].copy()
    out.insert(0, "source_PMID", "42384774")
    out["epitope_status"] = "not_reported_by_Olink_or_source_supplement"
    out["protein_form_interpretation"] = "assay target annotation only; soluble/cleaved form not assumed"
    out["analysis_role"] = np.where(
        out["gene_name"].isin(["TREM2", "ACE"]),
        "supports_cross_platform_mapping_of_existing_assay",
        "external_candidate_pending_available_summary_statistics",
    )
    out.to_csv(ROOT / "tables" / "moderate_target_assay_crosswalk.tsv", sep="\t", index=False)

    alt = pd.DataFrame([
        ["CLN5", "deCODE SomaScan; AGES blood proteome", "PMID40452368; PMID34381170", "yes", "no", "Published beta only; no rs429358/rs7412 alpha for the same assay."],
        ["COL10A1", "deCODE SomaScan", "PMID40452368; PMID40397384", "yes", "no", "No corrected-significant AD beta in retrieved tables; alpha unavailable."],
        ["PLOD2", "deCODE SomaScan", "PMID40452368; PMID40397384", "yes", "no", "No corrected-significant AD beta in retrieved tables; alpha unavailable."],
        ["SDF2", "deCODE SomaScan; brain pQTL", "PMID40452368; PMID40397384; PMID34381170", "yes", "no", "Circulating and brain evidence kept separate; alpha unavailable."],
        ["TMEM106B", "deCODE SomaScan; AGES blood proteome; brain pQTL", "PMID40452368; PMID40397384; PMID34381170", "yes", "no", "Corrected-significant published AD beta exists on deCODE, but same-assay alpha is unavailable."],
        ["VTN", "deCODE SomaScan; brain pQTL", "PMID40452368; PMID40397384; PMID34381170", "yes", "no", "Circulating published beta is not corrected-significant; alpha unavailable."],
    ], columns=["gene_symbol", "alternative_pQTL_platform", "verified_sources", "platform_evidence_available", "eligible_for_APOE_mediation", "boundary_reason"])
    alt.to_csv(ROOT / "tables" / "six_assay_unavailable_genes_alternative_pQTL_status.tsv", sep="\t", index=False)

def update_flow_and_manifest() -> None:
    flow_path = TABLE_DIR / "TableS24_APOE_Linkable_Eligibility_Flow.tsv"
    flow = pd.read_csv(flow_path, sep="\t", dtype=str).fillna("")
    extra = pd.DataFrame([
        ["priority_studies_newly_fulltext_verified", "4", "studies", "All four previously pending priority studies now have verified publisher supplements."],
        ["new_Bonferroni_PWAS_candidates_external_only", "19", "assay-level PWAS rows", "PWAS Z statistics are not standardized beta estimates and were not inserted into pooled mediation."],
        ["prespecified_five_protein_crosswalk_extension_completed", "5", "assays", "BCAM, CD55, LILRB1, LILRB5 and SCARA5 have complete readable provider-authorized assay archives and completed strict-QC alpha/beta/cis/mediation re-estimation as a sensitivity-only extension; they do not backfill the prespecified 25-protein primary panel."],
    ], columns=flow.columns)
    flow = pd.concat([flow[~flow["stage"].isin(extra["stage"])], extra], ignore_index=True)
    flow.to_csv(flow_path, sep="\t", index=False)

    manifest_path = TABLE_DIR / "TableS29_Figure5_Source_Data_and_QA_Manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    additions = []
    for label, path in {
        "PriorityStudyDownloadManifest": ROOT / "logs" / "priority_study_download_manifest.tsv",
        "PriorityStudyTargetHits": ROOT / "data_processed" / "priority_study_supplements" / "priority_study_target_hits.tsv",
        "ModerateTargetAssayCrosswalk": ROOT / "tables" / "moderate_target_assay_crosswalk.tsv",
        "AlternativePQTLStatus": ROOT / "tables" / "six_assay_unavailable_genes_alternative_pQTL_status.tsv",
        "MultiCausalColocStatus": ROOT / "tables" / "multicausal_coloc_status.tsv",
        "MultiCausalColocSignals": ROOT / "tables" / "multicausal_coloc_signal_pairs.tsv",
        "CovarianceMappingBootstrap": ROOT / "tables" / "covariance_mapping_bootstrap_sensitivity.tsv",
    }.items():
        if not path.exists():
            raise RuntimeError(f"Manifest source missing: {path}")
        if path.suffix.lower() in {".csv", ".tsv", ".txt", ".md", ".r", ".py"}:
            rows = max(0, sum(1 for _ in path.open(encoding="utf-8", errors="replace")) - 1)
        else:
            rows = "NA"
        try:
            source_file = str(path.relative_to(ROOT))
        except ValueError:
            source_file = str(path.relative_to(ROOT.parent))
        additions.append({
            "supplementary_table": label, "output_file": path.name, "source_file": source_file,
            "row_count": rows, "bytes": path.stat().st_size, "sha256": sha256(path), "status": "verified_machine_readable_source",
        })
    manifest = manifest[~manifest["supplementary_table"].isin([row["supplementary_table"] for row in additions])]
    manifest = pd.concat([manifest, pd.DataFrame(additions)], ignore_index=True)
    manifest.to_csv(manifest_path, sep="\t", index=False)


def main() -> None:
    update_screening()
    append_provenance()
    write_crosswalk_and_boundaries()
    update_flow_and_manifest()
    print("Integrated four priority studies and updated S18/S19/S24/S29 plus boundary tables")


if __name__ == "__main__":
    main()
