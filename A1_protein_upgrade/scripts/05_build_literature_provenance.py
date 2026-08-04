#!/usr/bin/env python3
"""Build a traceable literature-protein evidence table from verified supplements."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FULLTEXT = ROOT / "data_raw" / "literature_fulltext"
PROCESSED = ROOT / "data_processed" / "literature_supplements_extracted"
OUTPUT = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
SUMMARY_OUTPUT = ROOT / "tables" / "literature_protein_evidence_summary.tsv"

REQUIRED_COLUMNS = [
    "record_id", "protein_name", "gene_symbol", "UniProt_ID", "assay_target_name",
    "assay_target_ID", "proteomic_platform", "protein_form_or_isoform", "disease",
    "disease_subtype", "outcome_definition", "source_study_first_author",
    "publication_year", "journal", "PMID", "DOI", "publication_status", "pQTL_source",
    "pQTL_sample_size", "ancestry", "cis_trans_status", "instrument_definition",
    "instrument_SNPs", "MR_method", "beta_or_logOR", "SE", "OR", "CI_lower",
    "CI_upper", "P_value", "multiple_testing_method", "corrected_significance",
    "colocalization_method", "PP_H4", "SMR_P", "HEIDI_P", "replication_dataset",
    "replication_status", "effect_direction", "outcome_GWAS", "overlap_with_primary_outcome_GWAS",
    "same_outcome_reanalysis", "evidence_tier", "high_complexity_locus",
    "assay_cross_reactivity_risk", "inclusion_status", "exclusion_reason",
    "APOE_rs429358_alpha_available", "APOE_rs7412_alpha_available", "alpha_source",
    "eligible_for_two_step_MR", "manual_verification_status", "notes",
]

EXTRA_COLUMNS = [
    "source_table_file", "source_sheet", "source_row", "beta_source_type",
    "mapping_confidence", "evidence_independence", "duplicate_group",
    "allele_harmonization_status", "derived_field_notes",
]

SOURCE_STUDY_IDS = {
    "PMID38903171", "PMID41064534", "PMID40327005", "PMID39408566",
    "PMID40037332", "PMID40050982", "PMID41350740", "PMID41630303", "PMID36510323",
}


def load_study_metadata() -> dict[str, dict[str, object]]:
    search_path = ROOT / "literature" / "search_results_deduplicated.tsv"
    search = pd.read_csv(search_path, sep="\t", dtype=str).fillna("NA")
    search = search[search["record_id"].isin(SOURCE_STUDY_IDS)].drop_duplicates("record_id")
    metadata: dict[str, dict[str, object]] = {}
    for _, value in search.iterrows():
        metadata[value["record_id"]] = {
            "first_author": value["first_author"],
            "year": value["publication_year"],
            "journal": value["journal"],
            "pmid": value["PMID"],
            "doi": value["DOI"],
        }
    missing_ids = SOURCE_STUDY_IDS - set(metadata)
    if missing_ids:
        raise ValueError(f"Study metadata missing from deduplicated search records: {sorted(missing_ids)}")
    return metadata


STUDIES = load_study_metadata()

GENE_ALIASES = {
    "TIMD3": "HAVCR2", "ACADV": "ACADVL", "Apo E": "APOE", "Apo E2": "APOE",
    "Siglec-3": "CD33", "a2-Antiplasmin": "SERPINF2",
}


def missing(value: object) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or str(value).strip() == ""


def text_value(value: object, default: str = "NA") -> str:
    return default if missing(value) else str(value).strip()


def num(value: object) -> float:
    return float(pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0])


def direction(beta: object) -> str:
    value = num(beta)
    if math.isnan(value):
        return "not_reported"
    return "higher_protein_higher_risk" if value > 0 else "higher_protein_lower_risk" if value < 0 else "null"


def is_complex(gene: str) -> bool:
    gene = gene.upper().strip()
    return bool(
        gene in {"APOE", "TOMM40", "APOC1", "NECTIN2"}
        or gene.startswith("CFH")
        or gene.startswith("HLA")
        or gene.startswith("IGH")
    )


def base_row(record_id: str, **updates: object) -> dict[str, object]:
    study = STUDIES[record_id]
    row: dict[str, object] = {column: "NA" for column in REQUIRED_COLUMNS + EXTRA_COLUMNS}
    row.update(
        record_id=record_id,
        source_study_first_author=study["first_author"],
        publication_year=study["year"],
        journal=study["journal"],
        PMID=study["pmid"],
        DOI=study["doi"],
        publication_status="peer_reviewed",
        APOE_rs429358_alpha_available="not_assessed",
        APOE_rs7412_alpha_available="not_assessed",
        alpha_source="not_assessed",
        eligible_for_two_step_MR="false_pending_alpha_and_reestimated_beta",
        manual_verification_status="full_text_and_source_table_verified",
        beta_source_type="published_beta",
        allele_harmonization_status="not_applicable_to_published_beta",
    )
    row.update(updates)
    gene = text_value(row.get("gene_symbol"), "mapping_unresolved")
    row["high_complexity_locus"] = "true" if gene != "mapping_unresolved" and is_complex(gene) else "false"
    row["effect_direction"] = direction(row.get("beta_or_logOR"))
    return row


def assign_tier(row: dict[str, object], coloc_threshold: float = 0.8) -> None:
    corrected = str(row["corrected_significance"]).lower() == "true"
    cis = str(row["cis_trans_status"]).startswith("cis") or "top_cis" in str(row["cis_trans_status"])
    mapping = row["mapping_confidence"] == "high"
    pp = num(row["PP_H4"])
    coloc = not math.isnan(pp) and pp >= coloc_threshold
    replication = str(row["replication_status"]).startswith("independent_")
    complex_locus = row["high_complexity_locus"] == "true"
    if corrected and cis and mapping and (coloc or replication) and not complex_locus:
        row["evidence_tier"] = "Tier1"
        row["inclusion_status"] = "primary_high_confidence_panel"
        row["exclusion_reason"] = "NA"
    else:
        reasons = []
        if not corrected:
            reasons.append("not_multiple_testing_corrected")
        if not cis:
            reasons.append("not_cis_only")
        if not mapping:
            reasons.append("mapping_not_high_confidence")
        if not (coloc or replication):
            reasons.append("no_strong_colocalization_or_independent_replication")
        if complex_locus:
            reasons.append("high_complexity_locus_requires_downgrade")
        row["evidence_tier"] = "Tier2"
        row["inclusion_status"] = "literature_source_universe_tier2"
        row["exclusion_reason"] = ";".join(reasons) or "tier2_by_prespecified_rule"
    if str(row["gene_symbol"]).upper() == "APOE":
        row["inclusion_status"] = "protein_outcome_evidence_excluded_from_APOE_mediation"
        row["exclusion_reason"] = "APOE_protein_default_exclusion_due_to_cis_circularity"


def parse_or_ci(value: object) -> tuple[float, float, float]:
    match = re.search(r"([0-9.eE+-]+)\s*\(([0-9.eE+-]+)\s*[-–,]\s*([0-9.eE+-]+)\)", text_value(value, ""))
    if not match:
        return math.nan, math.nan, math.nan
    return tuple(float(match.group(i)) for i in range(1, 4))


def extract_pu() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID38903171_PMC11187347" / "supplements" / "Table_1.XLSX"
    rows: list[dict[str, object]] = []
    specs = [
        ("deCODE", "Supplementary Table 2", "Supplementary Table 10", "any_AMD", "FinnGen_R10_H7_AMD"),
        ("deCODE", "Supplementary Table 3", "Supplementary Table 11", "dry_AMD", "FinnGen_R10_DRY_AMD"),
        ("deCODE", "Supplementary Table 4", "Supplementary Table 12", "wet_AMD", "FinnGen_R10_WET_AMD"),
        ("UKB-PPP_Olink", "Supplementary Table 6", "Supplementary Table 13", "any_AMD", "FinnGen_R10_H7_AMD"),
        ("UKB-PPP_Olink", "Supplementary Table 7", "Supplementary Table 14", "dry_AMD", "FinnGen_R10_DRY_AMD"),
        ("UKB-PPP_Olink", "Supplementary Table 8", "Supplementary Table 15", "wet_AMD", "FinnGen_R10_WET_AMD"),
    ]
    for platform, mr_sheet, coloc_sheet, subtype, outcome in specs:
        mr = pd.read_excel(path, sheet_name=mr_sheet, header=1)
        mr["p_SMR_fdr"] = pd.to_numeric(mr["p_SMR_fdr"], errors="coerce")
        mr["p_HEIDI"] = pd.to_numeric(mr["p_HEIDI"], errors="coerce")
        mr = mr[(mr["p_SMR_fdr"] < 0.05) & (mr["p_HEIDI"] > 0.01)].copy()
        coloc = pd.read_excel(path, sheet_name=coloc_sheet, header=1)
        coloc.columns = [str(c).strip() for c in coloc.columns]
        protein_col = "proteins" if "proteins" in coloc.columns else "protein"
        pp_map = dict(zip(coloc[protein_col].astype(str).str.upper(), pd.to_numeric(coloc["PP.H4.abf"], errors="coerce")))
        for index, value in mr.iterrows():
            gene = text_value(value.get("Gene"), "mapping_unresolved").upper()
            pp = pp_map.get(gene, math.nan)
            row = base_row(
                "PMID38903171", protein_name=gene, gene_symbol=gene, UniProt_ID="not_reported",
                assay_target_name=gene, assay_target_ID=text_value(value.get("probeID")),
                proteomic_platform=platform, protein_form_or_isoform="not_reported",
                disease="AMD", disease_subtype=subtype,
                outcome_definition=f"FinnGen R10 {subtype}", pQTL_source=platform,
                pQTL_sample_size="35559" if platform == "deCODE" else "54219",
                ancestry="European", cis_trans_status="top_cis_pQTL",
                instrument_definition="top cis-pQTL within +/-1 Mb; P<5e-8; SMR",
                instrument_SNPs=text_value(value.get("topSNP")), MR_method="SMR",
                beta_or_logOR=value.get("b_SMR"), SE=value.get("se_SMR"),
                OR=np.exp(num(value.get("b_SMR"))), CI_lower="not_reported", CI_upper="not_reported",
                P_value=value.get("p_SMR"), multiple_testing_method="Benjamini-Hochberg FDR",
                corrected_significance="true", colocalization_method="coloc.abf single-causal-variant",
                PP_H4=pp, SMR_P=value.get("p_SMR"), HEIDI_P=value.get("p_HEIDI"),
                replication_dataset="cross-platform and FinnGen replication tables reported",
                replication_status="partially_overlapping_cross_platform",
                outcome_GWAS=outcome, overlap_with_primary_outcome_GWAS="same_FinnGen_resource_different_release_R10_vs_R12",
                same_outcome_reanalysis="true_different_release", assay_cross_reactivity_risk="requires_platform_review",
                notes="Only FDR<0.05 and HEIDI>0.01 rows retained; colocalization joined by reported gene label.",
                source_table_file=str(path.relative_to(ROOT)), source_sheet=mr_sheet,
                source_row=int(index) + 3, mapping_confidence="high", evidence_independence="partially_overlapping",
                duplicate_group="AMD_FinnGen_UKB_deCODE", derived_field_notes="OR=exp(reported SMR beta)",
            )
            assign_tier(row, coloc_threshold=0.75)
            rows.append(row)
    return rows


def extract_belbasis() -> list[dict[str, object]]:
    path = PROCESSED / "PMID40037332" / "awaf018_supplementary_data" / "brain-2024-00441-File011.xlsx"
    mr = pd.read_excel(path, sheet_name="Table S2", header=2)
    mr["FDR"] = pd.to_numeric(mr["FDR"], errors="coerce")
    mr = mr[(mr["Outcome"] == "Alzheimer's disease") & (mr["FDR"] < 0.05)]
    coloc = pd.read_excel(path, sheet_name="Table S5", header=2)
    coloc = coloc[coloc["Outcome"] == "Alzheimer's disease"]
    pp_assay = dict(zip(coloc["Assay ID"].astype(str), pd.to_numeric(coloc["PP.H4.abf"], errors="coerce")))
    rows = []
    for index, value in mr.iterrows():
        gene = text_value(value["Protein"], "mapping_unresolved").upper()
        assay = text_value(value["Assay ID"])
        row = base_row(
            "PMID40037332", protein_name=gene, gene_symbol=gene, UniProt_ID=text_value(value["UniProt ID"]),
            assay_target_name=gene, assay_target_ID=assay, proteomic_platform=text_value(value["Platform"]),
            protein_form_or_isoform="assay_specific_target", disease="AD", disease_subtype="AD",
            outcome_definition="clinically diagnosed Alzheimer disease GWAS", pQTL_source=text_value(value["Platform"]),
            pQTL_sample_size="UKB-PPP 54219 or deCODE 35559 according to platform",
            ancestry="European", cis_trans_status="cis_lead_pQTL",
            instrument_definition="lead cis-pQTL per assay", instrument_SNPs=text_value(value["rsID"]),
            MR_method="Wald ratio using lead cis-pQTL", beta_or_logOR=value["Beta"], SE=value["SE"],
            OR=np.exp(num(value["Beta"])), CI_lower=np.exp(num(value["Beta"]) - 1.96 * num(value["SE"])),
            CI_upper=np.exp(num(value["Beta"]) + 1.96 * num(value["SE"])), P_value=value["P"],
            multiple_testing_method="Benjamini-Hochberg FDR", corrected_significance="true",
            colocalization_method="coloc.abf single-causal-variant", PP_H4=pp_assay.get(assay, math.nan),
            SMR_P="not_reported", HEIDI_P="not_reported", replication_dataset="cross-platform comparison reported",
            replication_status="partially_overlapping_cross_platform", outcome_GWAS="AD GWAS reported in Table S1",
            overlap_with_primary_outcome_GWAS="likely_overlap_with_Wightman_2021_components",
            same_outcome_reanalysis="possible", assay_cross_reactivity_risk="assay_specific_review_required",
            notes="FDR-significant AD row; colocalization joined by assay ID.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="Table S2", source_row=int(index) + 4,
            mapping_confidence="high", evidence_independence="partially_overlapping",
            duplicate_group="AD_UKB_deCODE_cross_platform", derived_field_notes="OR and CI derived from reported beta and SE",
        )
        assign_tier(row)
        rows.append(row)
    return rows


def extract_zhan() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID40050982_PMC11884171" / "supplements" / "12967_2025_6317_MOESM1_ESM.xlsx"
    discovery = pd.read_excel(path, sheet_name="Table S2", header=3)
    discovery = discovery.rename(columns={"Unnamed: 2": "Outcome", "Unnamed: 8": "FDR"})
    discovery["Protein_group"] = discovery["Protein"].ffill()
    discovery["UniProt_group"] = discovery["UniProt_ID"].ffill()
    discovery["Outcome_group"] = discovery["Outcome"].ffill()
    discovery["FDR"] = pd.to_numeric(discovery["FDR"], errors="coerce")
    primary_methods = {"Inverse variance weighted", "Wald ratio"}
    primary = discovery[
        discovery["Method"].isin(primary_methods)
        & (discovery["Outcome_group"] == "AD")
        & (discovery["FDR"] < 0.05)
    ].copy()

    replication = pd.read_excel(path, sheet_name="Table S4", header=3)
    replication = replication.rename(columns={"Unnamed: 0": "Protein", "Unnamed: 1": "Outcome"})
    replication["Protein_group"] = replication["Protein"].ffill()
    replication["Outcome_group"] = replication["Outcome"].ffill()
    replication_primary = replication[
        replication["Method"].isin(primary_methods) & (replication["Outcome_group"] == "AD")
    ].copy()
    replication_map = {str(value["Protein_group"]): value for _, value in replication_primary.iterrows()}

    coloc = pd.read_excel(path, sheet_name="Table S5", header=3)
    coloc = coloc.rename(columns={"Unnamed: 0": "Protein"})
    coloc_map = {str(value["Protein"]): value for _, value in coloc.iterrows() if not missing(value["Protein"])}

    rows = []
    for index, value in primary.iterrows():
        protein = text_value(value["Protein_group"])
        uniprot = text_value(value["UniProt_group"])
        or_value, ci_lower, ci_upper = parse_or_ci(value["OR(95%CI)"])
        beta = math.log(or_value) if not math.isnan(or_value) and or_value > 0 else math.nan
        replication_row = replication_map.get(protein)
        repl_p = num(replication_row["Pval"]) if replication_row is not None else math.nan
        repl_or, _, _ = parse_or_ci(replication_row["OR(95%CI)"]) if replication_row is not None else (math.nan, math.nan, math.nan)
        same_direction = not math.isnan(repl_or) and not math.isnan(or_value) and np.sign(math.log(repl_or)) == np.sign(beta)
        replicated = repl_p < 0.00208 and same_direction
        coloc_row = coloc_map.get(protein)
        pp_discovery = num(coloc_row["PP.H4"]) if coloc_row is not None else math.nan
        pp_replication = num(coloc_row["PP.H4.1"]) if coloc_row is not None else math.nan
        pp = min(pp_discovery, pp_replication) if not math.isnan(pp_discovery) and not math.isnan(pp_replication) else pp_discovery
        instrument_rows = discovery[discovery["Protein_group"] == protein]
        snps = [str(snp) for snp in instrument_rows["SNP"].dropna() if str(snp).startswith("rs")]
        row = base_row(
            "PMID40050982", protein_name=protein, gene_symbol=protein.upper(), UniProt_ID=uniprot,
            assay_target_name=protein, assay_target_ID=uniprot, proteomic_platform="UKB-PPP Olink",
            protein_form_or_isoform="UniProt-mapped protein", disease="AD", disease_subtype="AD",
            outcome_definition="Bellenguez et al. 2022 AD discovery; FinnGen R11 replication",
            pQTL_source="UKB-PPP cis-pQTL", pQTL_sample_size="34557 European participants",
            ancestry="European", cis_trans_status="cis_pQTL",
            instrument_definition="cis-pQTL; discovery MR with Bellenguez AD; FinnGen R11 outcome replication",
            instrument_SNPs=";".join(snps), MR_method=text_value(value["Method"]),
            beta_or_logOR=beta, SE="not_reported_directly", OR=or_value, CI_lower=ci_lower, CI_upper=ci_upper,
            P_value=value["Pval"], multiple_testing_method="Benjamini-Hochberg FDR in discovery; P<0.00208 replication threshold",
            corrected_significance="true", colocalization_method="coloc.abf in discovery and replication",
            PP_H4=pp, SMR_P="not_reported", HEIDI_P="not_reported", replication_dataset="FinnGen R11 AD",
            replication_status="independent_outcome_replication_corrected_same_direction" if replicated else "replication_not_corrected_or_direction_inconsistent",
            outcome_GWAS="Bellenguez 2022 discovery; FinnGen R11 replication",
            overlap_with_primary_outcome_GWAS="Bellenguez_partially_overlaps_Wightman_components; FinnGen_independent_of_primary_AD_outcome_GWAS",
            same_outcome_reanalysis="partially_overlapping_discovery_with_independent_outcome_replication",
            assay_cross_reactivity_risk="Olink assay review required",
            notes=f"Replication P={repl_p if not math.isnan(repl_p) else 'NA'}; replication PP.H4={pp_replication if not math.isnan(pp_replication) else 'NA'}.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="Table S2", source_row=int(index) + 5,
            mapping_confidence="high", evidence_independence="independent_outcome_replication" if replicated else "partially_overlapping",
            duplicate_group="AD_UKB_PPP_Bellenguez_FinnGen", derived_field_notes="logOR derived from rounded reported OR; not eligible for mediation beta",
        )
        assign_tier(row)
        rows.append(row)
    return rows


def extract_zhou() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID41350740_PMC12681122" / "supplements" / "13195_2025_1901_MOESM1_ESM.xlsx"
    mr = pd.read_excel(path, sheet_name="Table S6", header=1)
    mr["FDR"] = pd.to_numeric(mr["FDR"], errors="coerce")
    selected = {
        "PILRA isoform FDF03-M14": "PILRA", "PILRA isoform FDF03-deltaTM": "PILRA",
        "GRN": "GRN", "ACE": "ACE", "TIMD3": "HAVCR2", "TREM2": "TREM2",
        "ACADV": "ACADVL", "OMGP": "OMGP", "BIN1": "BIN1",
    }
    mr = mr[mr["Exposure"].isin(selected) & (mr["FDR"] < 0.05)]
    sens_raw = pd.read_excel(path, sheet_name="Table S24", header=None)
    sens = sens_raw.iloc[3:, :6].copy()
    sens.columns = ["Exposure", "Outcome", "PPH4", "HEIDI_P", "HEIDI_nsnp", "SMR_P"]
    sens["Exposure"] = sens["Exposure"].ffill()
    pp_map = dict(zip(sens["Exposure"].astype(str), pd.to_numeric(sens["PPH4"], errors="coerce")))
    heidi_map = dict(zip(sens["Exposure"].astype(str), pd.to_numeric(sens["HEIDI_P"], errors="coerce")))
    paper_tier1 = {"PILRA isoform FDF03-M14", "PILRA isoform FDF03-deltaTM", "GRN", "ACE", "TIMD3", "TREM2"}
    rows = []
    for index, value in mr.iterrows():
        protein = text_value(value["Exposure"])
        gene = selected[protein]
        pp = pp_map.get(protein, math.nan)
        replication = "independent_external_replication_reported" if protein in paper_tier1 else "single_external_replication_reported"
        row = base_row(
            "PMID41350740", protein_name=protein, gene_symbol=gene, UniProt_ID="not_reported",
            assay_target_name=protein, assay_target_ID="deCODE protein label requires cross-platform mapping",
            proteomic_platform="deCODE SomaScan", protein_form_or_isoform=protein if "isoform" in protein else "not_reported",
            disease="AD", disease_subtype="AD", outcome_definition="European AD meta-analysis; AD GWAS 1",
            pQTL_source="deCODE", pQTL_sample_size="35559", ancestry="European", cis_trans_status="cis_pQTL",
            instrument_definition="cis-pQTL within +/-1 Mb; P<5e-8; MHC excluded; r2<0.001",
            instrument_SNPs="reported in Table S3", MR_method=text_value(value["Method"]),
            beta_or_logOR=value["b"], SE=value["se"], OR=value["or"], CI_lower=value["or_lci95"],
            CI_upper=value["or_uci95"], P_value=value["pval"], multiple_testing_method="Benjamini-Hochberg FDR",
            corrected_significance="true", colocalization_method="Bayesian coloc within +/-1 Mb",
            PP_H4=pp, SMR_P="reported in Table S24", HEIDI_P=heidi_map.get(protein, math.nan),
            replication_dataset="IGAP; FinnGen; ARIC; INTERVAL as available",
            replication_status=replication, outcome_GWAS="AD GWAS 1 (398058 participants)",
            overlap_with_primary_outcome_GWAS="partially_overlaps_or_reuses_components_of_Wightman_2021",
            same_outcome_reanalysis="possible", assay_cross_reactivity_risk="SomaScan assay mapping review required",
            notes="Source-paper tier was not copied automatically; raw FDR, colocalization and replication fields drive study-specific evidence tiering.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="Table S6", source_row=int(index) + 3,
            mapping_confidence="high" if protein not in {"TIMD3", "ACADV"} else "moderate",
            evidence_independence="partially_overlapping", duplicate_group="AD_deCODE_multiomics_2025",
            derived_field_notes="none",
        )
        assign_tier(row, coloc_threshold=0.75)
        rows.append(row)
    return rows


def extract_hou() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID41064534_PMC12503973" / "supplements" / "3972293.f6.xlsx"
    data = pd.read_excel(path, sheet_name="Supplementary_table6", header=1)
    rows = []
    for index, value in data.iterrows():
        disc_fdr = num(value["FDR__discovery"])
        repl_fdr = num(value["FDR__replication"])
        disc_beta, repl_beta = num(value["BETA__discovery"]), num(value["BETA__replication"])
        corrected = disc_fdr < 0.05 and repl_fdr < 0.05 and np.sign(disc_beta) == np.sign(repl_beta)
        gene = text_value(value["Genename"]).upper()
        row = base_row(
            "PMID41064534", protein_name=gene, gene_symbol=gene, UniProt_ID="not_reported",
            assay_target_name=gene, assay_target_ID="mapping_requires_source_pQTL_table",
            proteomic_platform="UKB-PPP Olink discovery; deCODE where indicated",
            protein_form_or_isoform="not_reported", disease="AMD", disease_subtype="dry_AMD",
            outcome_definition="FinnGen R11 dry AMD discovery; MVP dry AMD replication",
            pQTL_source="UKB-PPP and/or deCODE", pQTL_sample_size="54219 and/or 35559",
            ancestry="European", cis_trans_status="cis_pQTL",
            instrument_definition="cis-pQTL within 1 Mb; P<5e-8; r2<0.1",
            instrument_SNPs="reported in instrument supplements", MR_method=text_value(value["Method_discovery"]),
            beta_or_logOR=disc_beta, SE=value["SE__discovery"], OR=np.exp(disc_beta),
            CI_lower=np.exp(disc_beta - 1.96 * num(value["SE__discovery"])),
            CI_upper=np.exp(disc_beta + 1.96 * num(value["SE__discovery"])), P_value=value["P__discovery"],
            multiple_testing_method="Benjamini-Hochberg FDR in discovery and replication",
            corrected_significance="true" if corrected else "false", colocalization_method="not_reported",
            PP_H4="not_reported", SMR_P="not_reported", HEIDI_P="not_reported",
            replication_dataset=text_value(value["Replication_dataset"]),
            replication_status="independent_outcome_replication_corrected_same_direction" if corrected else "replication_not_corrected_or_direction_inconsistent",
            outcome_GWAS="FinnGen R11 and MVP", overlap_with_primary_outcome_GWAS="same_FinnGen_resource_different_release_R11_vs_R12",
            same_outcome_reanalysis="true_discovery_different_release; independent_MVP_replication",
            assay_cross_reactivity_risk="requires_platform_review",
            notes=f"Replication beta={repl_beta}; replication FDR={repl_fdr}. Table title alone was not treated as proof of replication.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="Supplementary_table6", source_row=int(index) + 3,
            mapping_confidence="moderate", evidence_independence="independent_outcome_for_replication",
            duplicate_group="DryAMD_FinnGen_R11_MVP", derived_field_notes="OR and CI derived from discovery beta and SE",
        )
        assign_tier(row)
        rows.append(row)
    return rows


def extract_li_ocular() -> list[dict[str, object]]:
    path = PROCESSED / "PMID39408566" / "ijms-25-10236-s001" / "ijms-3215465-supplementary.xlsx"
    mr = pd.read_excel(path, sheet_name="S4_plasma_MR-IVW", header=1)
    mr["FDR.P n"] = pd.to_numeric(mr["FDR.P n"], errors="coerce")
    mr = mr[(mr["Ocular diseases"].astype(str).str.lower() == "age-related macular degeneration") & (mr["FDR.P n"] < 0.05)]
    coloc = pd.read_excel(path, sheet_name="S6_plasma_COLOC", header=1)
    coloc = coloc[coloc["Ocular diseases"].astype(str).str.lower() == "age-related macular degeneration"]
    pp_map = dict(zip(coloc["Protein"].astype(str), pd.to_numeric(coloc["PP.H4.abf e"], errors="coerce")))
    rows = []
    for index, value in mr.iterrows():
        gene = text_value(value["GENE b"]).replace("*", "").strip().upper()
        assay = text_value(value["Protein a"])
        row = base_row(
            "PMID39408566", protein_name=gene, gene_symbol=gene, UniProt_ID="not_reported",
            assay_target_name=gene, assay_target_ID=assay, proteomic_platform="deCODE SomaScan",
            protein_form_or_isoform="assay_specific_target", disease="AMD", disease_subtype="any_AMD",
            outcome_definition="age-related macular degeneration GWAS", pQTL_source="deCODE",
            pQTL_sample_size="35559", ancestry="European", cis_trans_status="requires_instrument_level_review",
            instrument_definition="genome-wide significant pQTLs; details in study methods",
            instrument_SNPs=f"{text_value(value['Nsnps k'])} SNPs; IDs in instrument table", MR_method="IVW",
            beta_or_logOR=value["BETA c"], SE=value["SE d"], OR=value["OR g"], CI_lower=value["ORLower h"],
            CI_upper=value["ORUpper i"], P_value=value["Pval j"], multiple_testing_method="Benjamini-Hochberg FDR",
            corrected_significance="true", colocalization_method="coloc.abf single-causal-variant",
            PP_H4=pp_map.get(assay, math.nan), SMR_P="not_reported", HEIDI_P="not_reported",
            replication_dataset="not_reported", replication_status="not_reported",
            outcome_GWAS="ocular disease GWAS listed in Table S1", overlap_with_primary_outcome_GWAS="requires_manual_verification",
            same_outcome_reanalysis="unclear", assay_cross_reactivity_risk="SomaScan assay review required",
            notes="AMD plasma MR row passing study FDR; colocalization joined by SeqId assay target.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="S4_plasma_MR-IVW", source_row=int(index) + 3,
            mapping_confidence="high", evidence_independence="partially_overlapping",
            duplicate_group="AMD_ocular_proteome_2024", derived_field_notes="none",
        )
        assign_tier(row)
        rows.append(row)
    return rows


def extract_chen() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID40327005_PMC12063708" / "supplements" / "tvst-14-5-8_s003.xlsx"
    specs = [
        ("Supplementary table 1", "deCODE", "any_AMD"), ("Supplementary table 2", "UKB-PPP Olink", "any_AMD"),
        ("Supplementary table 3", "deCODE", "early_AMD"), ("Supplementary table 4", "UKB-PPP Olink", "early_AMD"),
        ("Supplementary table 5", "deCODE", "dry_AMD"), ("Supplementary table 6", "UKB-PPP Olink", "dry_AMD"),
        ("Supplementary table 7", "deCODE", "wet_AMD"), ("Supplementary table 8", "UKB-PPP Olink", "wet_AMD"),
    ]
    rows = []
    for sheet, platform, subtype in specs:
        raw = pd.read_excel(path, sheet_name=sheet, header=None)
        headers = []
        grouped_headers = raw.iloc[1].ffill()
        for upper, lower in zip(grouped_headers, raw.iloc[2]):
            name = text_value(upper, "")
            sub = text_value(lower, "")
            headers.append(f"{name}_{sub}".strip("_") if sub else name)
        data = raw.iloc[3:].copy()
        data.columns = headers
        fdr_col = next(c for c in data.columns if c == "FDR" or c.endswith("_FDR"))
        heidi_col = next(c for c in data.columns if "p_HEIDI" in c)
        smr_p_col = next(c for c in data.columns if "SMR association_Pval" in c)
        smr_beta_col = next(c for c in data.columns if "SMR association_Beta" in c)
        smr_se_col = next(c for c in data.columns if "SMR association_Se" in c)
        data[fdr_col] = pd.to_numeric(data[fdr_col], errors="coerce")
        data[heidi_col] = pd.to_numeric(data[heidi_col], errors="coerce")
        data = data[(data[fdr_col] < 0.05) & (data[heidi_col] > 0.01)]
        for index, value in data.iterrows():
            gene = text_value(value["Gene"]).upper()
            beta = value[smr_beta_col]
            row = base_row(
                "PMID40327005", protein_name=gene, gene_symbol=gene, UniProt_ID="not_reported",
                assay_target_name=gene, assay_target_ID="gene-level label only", proteomic_platform=platform,
                protein_form_or_isoform="not_reported", disease="AMD", disease_subtype=subtype,
                outcome_definition=subtype, pQTL_source=platform,
                pQTL_sample_size="35559" if platform == "deCODE" else "54219", ancestry="European",
                cis_trans_status="top_cis_pQTL_SMR", instrument_definition="top cis-pQTL used in SMR",
                instrument_SNPs=text_value(value["topSNP"]), MR_method="SMR", beta_or_logOR=beta,
                SE=value[smr_se_col], OR=np.exp(num(beta)), CI_lower="not_reported", CI_upper="not_reported",
                P_value=value[smr_p_col], multiple_testing_method="Benjamini-Hochberg FDR",
                corrected_significance="true", colocalization_method="not_reported", PP_H4="not_reported",
                SMR_P=value[smr_p_col], HEIDI_P=value[heidi_col], replication_dataset="cross-platform comparison",
                replication_status="partially_overlapping_cross_platform", outcome_GWAS="study-specific AMD GWAS",
                overlap_with_primary_outcome_GWAS="requires_manual_verification", same_outcome_reanalysis="unclear",
                assay_cross_reactivity_risk="requires_platform_review",
                notes="Corrected SMR row; no colocalization table was reported.",
                source_table_file=str(path.relative_to(ROOT)), source_sheet=sheet, source_row=int(index) + 1,
                mapping_confidence="moderate", evidence_independence="partially_overlapping",
                duplicate_group="AMD_multi_platform_2025", derived_field_notes="OR=exp(reported SMR beta)",
            )
            assign_tier(row)
            rows.append(row)
    return rows


def extract_zhao() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID41630303_PMC12863894" / "supplements" / "medi-105-e45715-s001.xlsx"
    data = pd.read_excel(path, sheet_name="Table S3", header=1)
    rows = []
    for index, value in data.iterrows():
        protein = text_value(value["exposure"])
        gene = GENE_ALIASES.get(protein, protein).upper()
        row = base_row(
            "PMID41630303", protein_name=protein, gene_symbol=gene, UniProt_ID="not_reported",
            assay_target_name=protein, assay_target_ID=text_value(value["id.exposure"]),
            proteomic_platform="deCODE SomaScan", protein_form_or_isoform="assay_specific_target",
            disease="AD", disease_subtype="AD", outcome_definition="FinnGen AD",
            pQTL_source="deCODE", pQTL_sample_size="35559", ancestry="European",
            cis_trans_status="mixed_or_not_explicitly_cis_only", instrument_definition="P<5e-8; LD clumped instruments",
            instrument_SNPs=f"{text_value(value['nsnp'])} SNPs; IDs in Table S2", MR_method="IVW",
            beta_or_logOR=value["b"], SE=value["se"], OR=value["or"], CI_lower=value["or_lci95"],
            CI_upper=value["or_uci95"], P_value=value["pval"], multiple_testing_method="study-reported FDR<0.05",
            corrected_significance="true", colocalization_method="not_reported", PP_H4="not_reported",
            SMR_P="not_reported", HEIDI_P="not_reported", replication_dataset="not_reported",
            replication_status="not_reported", outcome_GWAS="FinnGen AD",
            overlap_with_primary_outcome_GWAS="different_outcome_GWAS_from_primary_analysis",
            same_outcome_reanalysis="false", assay_cross_reactivity_risk="SomaScan assay review required",
            notes="Paper states these 39 rows passed FDR; row-level adjusted P is absent, so retained in Tier2 only.",
            source_table_file=str(path.relative_to(ROOT)), source_sheet="Table S3", source_row=int(index) + 3,
            mapping_confidence="moderate", evidence_independence="partially_overlapping",
            duplicate_group="AD_deCODE_FinnGen_2026", derived_field_notes="none",
        )
        assign_tier(row)
        rows.append(row)
    return rows


def extract_yang() -> list[dict[str, object]]:
    path = FULLTEXT / "PMID36510323_PMC9746220" / "supplements" / "13073_2022_1140_MOESM2_ESM.xlsx"
    rows = []
    for sheet in ["S2", "S3", "S4", "S5", "S6", "S7"]:
        data = pd.read_excel(path, sheet_name=sheet)
        data.columns = [str(c).strip() for c in data.columns]
        outcome_col = next((c for c in data.columns if c.lower() in {"phenotype", "outcome", "trait"}), None)
        tissue_col = next((c for c in data.columns if c.lower() == "tissue"), None)
        if not outcome_col or not tissue_col:
            continue
        subset = data[
            data[outcome_col].astype(str).str.contains("Alzheimer", case=False, na=False)
            & data[tissue_col].astype(str).str.contains("plasma", case=False, na=False)
        ]
        for index, value in subset.iterrows():
            protein_col = next((c for c in data.columns if c.lower() in {"protein", "protein_name", "exposure"}), None)
            if not protein_col:
                continue
            protein = text_value(value[protein_col])
            gene = GENE_ALIASES.get(protein, protein).upper()
            beta_col = next((c for c in data.columns if c.lower() in {"beta", "b", "mr_beta"}), None)
            se_col = next((c for c in data.columns if c.lower() in {"se", "mr_se"}), None)
            p_col = next((c for c in data.columns if c.lower() in {"p", "pval", "p_value", "mr_p"}), None)
            pp_col = next((c for c in data.columns if "pp.h4" in c.lower()), None)
            row = base_row(
                "PMID36510323", protein_name=protein, gene_symbol=gene, UniProt_ID="not_reported",
                assay_target_name=protein, assay_target_ID="requires_source_mapping", proteomic_platform="SomaScan",
                protein_form_or_isoform="assay_specific_target", disease="AD", disease_subtype="AD",
                outcome_definition=text_value(value[outcome_col]), pQTL_source="WashU multi-tissue proteomics",
                pQTL_sample_size="529 plasma samples in discovery atlas", ancestry="European",
                cis_trans_status="cis_pQTL", instrument_definition="study-wide significant cis pQTL",
                instrument_SNPs="reported in supplement", MR_method="Wald ratio or IVW",
                beta_or_logOR=value[beta_col] if beta_col else "not_reported", SE=value[se_col] if se_col else "not_reported",
                OR="not_reported", CI_lower="not_reported", CI_upper="not_reported",
                P_value=value[p_col] if p_col else "not_reported", multiple_testing_method="FDR<0.05",
                corrected_significance="true", colocalization_method="coloc.abf and coloc.susie",
                PP_H4=value[pp_col] if pp_col else "not_reported", SMR_P="not_reported", HEIDI_P="not_reported",
                replication_dataset="previous plasma proteome-by-phenome study",
                replication_status="partially_overlapping_external_comparison", outcome_GWAS=text_value(value[outcome_col]),
                overlap_with_primary_outcome_GWAS="likely_overlap_with_primary_AD_GWAS_components",
                same_outcome_reanalysis="possible", assay_cross_reactivity_risk="SomaScan assay review required",
                notes="Extracted only from cis-pQTL workflow sheets; supplement schema requires final manual row audit.",
                source_table_file=str(path.relative_to(ROOT)), source_sheet=sheet, source_row=int(index) + 2,
                mapping_confidence="moderate", evidence_independence="partially_overlapping",
                duplicate_group="AD_multi_tissue_2022", derived_field_notes="none",
            )
            assign_tier(row)
            rows.append(row)
    return rows


def build_summary(frame: pd.DataFrame) -> pd.DataFrame:
    keys = ["gene_symbol", "protein_name", "disease", "disease_subtype"]
    summaries = []
    for key, group in frame.groupby(keys, dropna=False):
        tier1 = int((group["evidence_tier"] == "Tier1").sum())
        independent = group["evidence_independence"].astype(str).str.contains("independent").any()
        summaries.append({
            **dict(zip(keys, key)),
            "n_evidence_rows": len(group),
            "n_unique_studies": group["PMID"].nunique(),
            "n_tier1_rows": tier1,
            "best_evidence_tier": "Tier1" if tier1 else "Tier2",
            "evidence_independence_summary": "contains_independent_component" if independent else "overlap_or_independence_unclear",
            "source_PMIDs": ";".join(sorted(group["PMID"].astype(str).unique())),
            "platforms": ";".join(sorted(group["proteomic_platform"].astype(str).unique())),
            "notes": "Counts are not interpreted as independent replications without the independence flag.",
        })
    return pd.DataFrame(summaries)


def main() -> None:
    extractors = [extract_pu, extract_belbasis, extract_zhan, extract_zhou, extract_hou, extract_li_ocular, extract_chen, extract_zhao, extract_yang]
    rows: list[dict[str, object]] = []
    for extractor in extractors:
        extracted = extractor()
        print(f"{extractor.__name__}: {len(extracted)} rows")
        rows.extend(extracted)
    frame = pd.DataFrame(rows)
    for column in REQUIRED_COLUMNS + EXTRA_COLUMNS:
        if column not in frame:
            frame[column] = "NA"
    frame = frame[REQUIRED_COLUMNS + EXTRA_COLUMNS]
    frame = frame.replace({np.nan: "NA", "": "NA"})
    frame = frame.sort_values(["disease", "disease_subtype", "gene_symbol", "PMID", "assay_target_ID"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUTPUT, sep="\t", index=False)
    summary = build_summary(frame)
    summary.to_csv(SUMMARY_OUTPUT, sep="\t", index=False)
    print(f"Total evidence rows: {len(frame)}")
    print(f"Tier1 rows: {(frame['evidence_tier'] == 'Tier1').sum()}")
    print(f"Unique protein/gene labels: {frame['gene_symbol'].nunique()}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
