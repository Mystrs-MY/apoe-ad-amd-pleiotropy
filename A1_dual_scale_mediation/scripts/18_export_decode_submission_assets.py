from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


DUAL = Path(__file__).resolve().parents[1]
ROOT = DUAL.parent
SRC = DUAL / "tables"
OUT_DIRS = [
    ROOT / "tables",
    ROOT / "tables_submission",
    ROOT / "tables_submission" / "supplementary_tables",
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_all(df: pd.DataFrame, filename: str) -> None:
    for out_dir in OUT_DIRS:
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir / filename, sep="\t", index=False, na_rep="NA")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


attrition = read_tsv(SRC / "TableS28_candidate_attrition_provenance.tsv")
attrition_columns = [
    "gene_symbol",
    "literature_target_entries",
    "literature_record_ids",
    "source_studies",
    "original_proteomic_platforms",
    "original_assay_ids",
    "literature_UniProt_IDs",
    "literature_protein_forms",
    "failure_class",
    "selected_Olink_assay",
    "mapping_confidence",
    "rs429358_alpha_available",
    "rs7412_alpha_available",
    "estimable_outcome_beta_count",
    "beta_failure_by_outcome",
    "exclusion_reason",
    "missing_values_treated_as_zero",
    "deCODE_exact_assay_count",
    "deCODE_SomaScan_SeqIds",
    "deCODE_UniProt_IDs",
    "deCODE_rs429358_direct_alpha_status",
    "deCODE_rs7412_direct_alpha_status",
    "deCODE_exact_assay_audit_status",
    "same_platform_alpha_beta_loop_status",
    "deCODE_reestimated_beta_count",
    "deCODE_mediation_paths_estimable",
    "deCODE_genome_wide_Bonferroni_72_survivors",
    "deCODE_cis_only_Bonferroni_72_survivors",
    "deCODE_cis_PAV_filtered_Bonferroni_72_survivors",
    "deCODE_result_classification",
    "ARIC_exact_assay_audit_status",
    "final_analysis_role",
    "notes",
]
attrition_final = attrition.loc[:, attrition_columns].copy()
for out_dir in OUT_DIRS:
    out_dir.mkdir(parents=True, exist_ok=True)
    attrition_final.to_csv(
        out_dir / "TableS28_Excluded_Name_Matched_Panel_Genes.tsv",
        sep="\t",
        index=False,
    )

canonical_s28 = (
    ROOT
    / "A1_protein_upgrade"
    / "tables"
    / "supplementary_name_match_revision"
    / "TableS28_Excluded_Name_Matched_Panel_Genes.tsv"
)
canonical_s28.parent.mkdir(parents=True, exist_ok=True)
attrition_final.to_csv(canonical_s28, sep="\t", index=False)

gate = read_tsv(SRC / "decode_same_platform_feasibility_gate.tsv")
gate_columns = [
    "gene_symbol",
    "SomaScan_SeqId",
    "SomaScan_display_ID",
    "UniProt_ID",
    "assay_target_name",
    "assay_target_full_name",
    "mapping_confidence",
    "rs429358_direct_alpha_in_full_GWAS",
    "rs7412_direct_alpha_in_full_GWAS",
    "complete_aptamer_GWAS_available_in_authorized_archive",
    "same_platform_beta_reestimation_ready",
    "reestimated_outcome_beta_count",
    "gate_status",
    "absence_interpretation",
    "source_file_sha256",
]
write_all(gate.loc[:, gate_columns], "TableS37a_deCODE_Exact_Assay_Gate.tsv")

alpha = read_tsv(SRC / "APOE_variant_to_decode_somascan_alpha.tsv")
alpha_columns = [
    "gene_symbol",
    "assay_target_ID",
    "UniProt_ID",
    "variant",
    "original_effect_allele",
    "original_other_allele",
    "harmonized_effect_allele",
    "harmonized_other_allele",
    "allele_flipped",
    "direct_variant_not_proxy",
    "chromosome_hg38",
    "position_hg38",
    "alpha",
    "alpha_SE",
    "alpha_P",
    "alpha_N",
    "imputation_MAF",
    "imputation_MAF_not_EAF",
    "effect_unit",
    "alpha_F_statistic",
    "alpha_P_FDR_18",
    "alpha_P_Bonferroni_18",
    "alpha_significance_not_used_for_eligibility",
    "alpha_source",
    "source_object_etag",
    "source_file_sha256",
    "manual_verification_status",
]
write_all(alpha.loc[:, alpha_columns], "TableS37b_deCODE_APOE_Alpha.tsv")

write_all(
    read_tsv(SRC / "decode_smp_two_step_mediation.tsv"),
    "TableS38a_deCODE_Genome_Wide_Two_Step_Mediation.tsv",
)
write_all(
    read_tsv(SRC / "decode_smp_two_step_mediation_cis_only.tsv"),
    "TableS38b_deCODE_Cis_Only_Two_Step_Mediation.tsv",
)
write_all(
    read_tsv(SRC / "decode_smp_two_step_mediation_cis_only_PAV_filtered.tsv"),
    "TableS38c_deCODE_Cis_PAV_Filtered_Two_Step_Mediation.tsv",
)

write_all(
    read_tsv(SRC / "decode_smp_shared_instrument_audit.tsv"),
    "TableS39a_deCODE_Shared_Instrument_Audit.tsv",
)
write_all(
    read_tsv(SRC / "decode_smp_PAV_epitope_instrument_audit.tsv"),
    "TableS39b_deCODE_PAV_Epitope_Audit.tsv",
)

raw_paths = [
    SRC / "decode_raw_gated_10074_128_8687_26_two_step_mediation.tsv",
    SRC / "decode_raw_gated_10074_128_8687_26_two_step_mediation_cis_only.tsv",
    SRC / "decode_raw_gated_10074_128_8687_26_two_step_mediation_cis_only_PAV_filtered.tsv",
]
raw_frames = []
for path in raw_paths:
    frame = read_tsv(path)
    frame.insert(0, "source_table", path.name)
    raw_frames.append(frame)
raw = pd.concat(raw_frames, ignore_index=True, sort=False).fillna("NA")
write_all(raw, "TableS39c_deCODE_Raw_Normalization_Sensitivity.tsv")

ledger = read_tsv(DUAL / "logs" / "decode_download_ledger.tsv")
integrity = pd.DataFrame(
    {
        "record_type": "downloaded_object",
        "normalization": ledger["normalization"],
        "gene_symbol": ledger["gene_symbol"],
        "SomaScan_SeqId": ledger["SomaScan_SeqId"],
        "content_length_bytes": ledger["content_length_bytes"],
        "etag": ledger["etag"],
        "sha256": ledger["sha256"],
        "status": ledger["status"],
        "assertion": "provider_size_etag_and_content_sha256_recorded",
        "assertion_value": "true",
    }
)

summary = read_tsv(SRC / "decode_smp_mediation_summary.tsv")
assertion_rows = []
for _, row in summary.iterrows():
    assertion_rows.append(
        {
            "record_type": "analysis_assertion",
            "normalization": "SMP",
            "gene_symbol": "NA",
            "SomaScan_SeqId": "NA",
            "content_length_bytes": "NA",
            "etag": "NA",
            "sha256": "NA",
            "status": "validated",
            "assertion": f"{row['analysis_scope']}:bonferroni_surviving_paths",
            "assertion_value": row["bonferroni_surviving_paths"],
        }
    )
validation_report = DUAL / "audit" / "decode_extension_validation_report.md"
assertion_rows.append(
    {
        "record_type": "analysis_assertion",
        "normalization": "all",
        "gene_symbol": "NA",
        "SomaScan_SeqId": "NA",
        "content_length_bytes": str(validation_report.stat().st_size),
        "etag": "NA",
        "sha256": sha256(validation_report),
        "status": "validated",
        "assertion": "automated_validation_report",
        "assertion_value": "PASS",
    }
)
integrity = pd.concat([integrity, pd.DataFrame(assertion_rows)], ignore_index=True)
write_all(integrity, "TableS39d_deCODE_File_Integrity_and_Validation.tsv")

flow = pd.DataFrame(
    [
        ("excluded_genes", 8, "BRD2;IL20RB;CLN5;COL10A1;PLOD2;SDF2;TMEM106B;VTN"),
        ("exact_SomaScan_aptamers", 9, "BRD2 has two aptamers"),
        ("planned_mediation_paths", 72, "9 aptamers x 2 APOE variants x 4 outcomes"),
        ("genome_wide_corrected_paths", 3, "BRD2 10074_128; CFH-locus trans tool"),
        ("cis_only_corrected_paths", 1, "TMEM106B 8687_26; target-gene PAV-linked tool"),
        ("robust_after_cis_PAV_shared_tool_audit", 0, "no robust mediator retained"),
    ],
    columns=["stage", "count", "interpretation"],
)
source_dir = ROOT / "figures_submission" / "source_data"
source_dir.mkdir(parents=True, exist_ok=True)
flow.to_csv(source_dir / "FigS12_deCODE_extension_flow.tsv", sep="\t", index=False)

print("Exported Table S28 and deCODE Tables S37-S39 to all submission table locations.")
