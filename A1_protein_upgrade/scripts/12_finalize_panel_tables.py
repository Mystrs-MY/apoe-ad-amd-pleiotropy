#!/usr/bin/env python3
"""Finalize alpha/beta flags and construct new-versus-legacy comparison tables."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
ALPHA = ROOT / "tables" / "APOE_variant_to_literature_proteins_alpha.tsv"
BETA = ROOT / "tables" / "literature_panel_beta_results.tsv"
MEDIATION = ROOT / "tables" / "APOE_linkable_two_step_mediation.tsv"
MEDIATION_CIS = ROOT / "tables" / "APOE_linkable_two_step_mediation_cis_sensitivity.tsv"
MEDIATION_STRICT = ROOT / "tables" / "APOE_linkable_two_step_mediation_strict_sensitivity.tsv"
TARGET_MAPPING = ROOT / "tables" / "APOE_linkable_target_assay_mapping.tsv"
TOTALS = ROOT / "tables" / "APOE_variant_total_effects_current_A1.tsv"
LEGACY_SELECTION = ROOT.parent / "tables" / "TableS4_Protein_Selection.csv"
LEGACY_ALPHA = ROOT.parent / "04_protein_mr" / "C6_rs429358_effects.csv"
LEGACY_BETA = ROOT.parent / "04_protein_mr" / "C6_ukbppp_mr_all.csv"
LEGACY_MEDIATION = ROOT.parent / "05_mediation" / "D_nvmr_mediation_all.csv"
COMPARISON = ROOT / "tables" / "new_vs_legacy_panel_comparison.tsv"
CURRENT_EXCLUSIONS = ROOT / "tables" / "name_matched_panel_excluded_genes.tsv"

OUTCOMES = ["AD", "any_AMD", "dry_AMD", "wet_AMD"]
LEGACY_OUTCOME_MAP = {"AD": "AD", "Any_AMD": "any_AMD", "Dry_AMD": "dry_AMD", "Wet_AMD": "wet_AMD"}


def semicolon(values: pd.Series) -> str:
    items = sorted({str(value) for value in values if str(value) not in {"", "NA", "nan"}})
    return ";".join(items) if items else "NA"


def update_provenance(
    provenance: pd.DataFrame,
    alpha: pd.DataFrame,
    eligible: set[str],
    target_mapping: pd.DataFrame,
) -> pd.DataFrame:
    key_columns = ["gene_symbol", "protein_name", "protein_form_or_isoform"]
    mapped_keys = {
        tuple(row[column] for column in key_columns)
        for _, row in target_mapping[
            target_mapping["eligible_for_alpha_beta_triangulation"].astype(str).str.lower().eq("true")
        ].iterrows()
    }
    provenance_keys = provenance[key_columns].apply(tuple, axis=1)
    scope_mask = provenance["inclusion_status"].isin([
        "primary_high_confidence_panel", "protein_outcome_evidence_excluded_from_APOE_mediation"
    ])
    provenance.loc[scope_mask, ["APOE_rs429358_alpha_available", "APOE_rs7412_alpha_available"]] = "false"
    provenance.loc[scope_mask, "alpha_source"] = "not_available_or_target_mapping_unresolved"
    provenance.loc[provenance["inclusion_status"] == "primary_high_confidence_panel", "eligible_for_two_step_MR"] = "false"
    for rsid, column in [("rs429358", "APOE_rs429358_alpha_available"), ("rs7412", "APOE_rs7412_alpha_available")]:
        available = set(alpha.loc[(alpha["variant"] == rsid) & (alpha["availability_status"] == "direct_variant_available"), "gene_symbol"])
        primary_mask = provenance["inclusion_status"].isin(["primary_high_confidence_panel", "protein_outcome_evidence_excluded_from_APOE_mediation"])
        row_available = provenance["gene_symbol"].isin(available) & provenance_keys.isin(mapped_keys)
        provenance.loc[primary_mask, column] = row_available.loc[primary_mask].map({True: "true", False: "false"})
    source_map = alpha.loc[alpha["availability_status"] == "direct_variant_available"].groupby("gene_symbol")["alpha_source"].agg(semicolon)
    mask = provenance["gene_symbol"].isin(source_map.index) & provenance_keys.isin(mapped_keys)
    provenance.loc[mask, "alpha_source"] = provenance.loc[mask, "gene_symbol"].map(source_map)
    primary_mask = provenance["inclusion_status"] == "primary_high_confidence_panel"
    row_eligible = provenance["gene_symbol"].isin(eligible) & provenance_keys.isin(mapped_keys)
    provenance.loc[primary_mask, "eligible_for_two_step_MR"] = row_eligible.loc[primary_mask].map({True: "true", False: "false"})
    provenance.to_csv(PROVENANCE, sep="\t", index=False)
    return provenance


def add_beta_placeholders(
    beta: pd.DataFrame, primary_targets: pd.DataFrame, target_mapping: pd.DataFrame
) -> pd.DataFrame:
    if "protein_name" not in beta:
        beta["protein_name"] = beta["gene_symbol"]
    else:
        beta.loc[beta["protein_name"].isin(["NA", "", "nan"]), "protein_name"] = beta["gene_symbol"]
    if "protein_form_or_isoform" not in beta:
        beta["protein_form_or_isoform"] = "not_reported"
    if "availability_scope" not in beta:
        beta["availability_scope"] = "reestimated_subset"
    existing = set(zip(beta["gene_symbol"], beta["outcome"]))
    primary_genes = primary_targets[["gene_symbol"]].drop_duplicates().sort_values("gene_symbol")
    eligible_mapping = target_mapping[
        target_mapping["eligible_for_primary"].astype(str).str.lower().eq("true")
    ]
    planned_family_size = eligible_mapping[["standardized_gene_symbol", "UKB_PPP_OID"]].drop_duplicates().shape[0]
    mapped_oids = eligible_mapping.groupby("standardized_gene_symbol")["UKB_PPP_OID"].first().to_dict()
    placeholders = []
    for _, target in primary_genes.iterrows():
        for outcome in OUTCOMES:
            if (target["gene_symbol"], outcome) in existing:
                continue
            mapped = target["gene_symbol"] in mapped_oids
            placeholders.append({
                "gene_symbol": target["gene_symbol"], "protein_name": target["gene_symbol"],
                "protein_form_or_isoform": "Olink assay-level target", "outcome": outcome,
                "assay_target_ID": mapped_oids.get(target["gene_symbol"], "mapping_unresolved"),
                "method": "not_estimable", "method_role": "primary", "nsnp": 0,
                "beta": "NA", "SE": "NA", "P_value": "NA", "beta_status": "not_reestimable",
                "exclusion_reason": (
                    "No harmonizable protein instruments were available for the current outcome GWAS."
                    if mapped else "No eligible UKB-PPP assay was found after approved-symbol and alias audit."
                ),
                "beta_source": "published_beta_retained_in_provenance_only",
                "planned_family_size": planned_family_size, "availability_scope": "planned_primary_gene",
            })
    combined = pd.concat([beta, pd.DataFrame(placeholders)], ignore_index=True, sort=False).fillna("NA")
    combined = combined.sort_values(["outcome", "gene_symbol", "protein_name", "method_role", "method"])
    combined.to_csv(BETA, sep="\t", index=False)
    return combined


def distribution_record(panel: str, variable: str, outcome: str, values: pd.Series, notes: str) -> dict[str, object]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "record_type": "distribution_summary", "panel": panel, "metric": variable, "outcome": outcome,
        "gene_symbol": "NA", "in_literature_panel": "NA", "in_legacy_panel": "NA",
        "value": "NA", "n": len(numeric), "mean": numeric.mean() if len(numeric) else "NA",
        "median": numeric.median() if len(numeric) else "NA",
        "q1": numeric.quantile(0.25) if len(numeric) else "NA",
        "q3": numeric.quantile(0.75) if len(numeric) else "NA",
        "minimum": numeric.min() if len(numeric) else "NA", "maximum": numeric.max() if len(numeric) else "NA",
        "notes": notes,
    }


def build_comparison(
    provenance: pd.DataFrame,
    alpha: pd.DataFrame,
    beta: pd.DataFrame,
    mediation: pd.DataFrame,
    cis_mediation: pd.DataFrame,
    strict_mediation: pd.DataFrame,
    target_mapping: pd.DataFrame,
) -> None:
    new_genes = sorted(provenance.loc[provenance["inclusion_status"] == "primary_high_confidence_panel", "gene_symbol"].unique())
    legacy_selection = pd.read_csv(LEGACY_SELECTION, dtype=str)
    legacy_genes = sorted(legacy_selection["Gene"].unique())
    union = sorted(set(new_genes) | set(legacy_genes))
    overlap = set(new_genes) & set(legacy_genes)
    rows: list[dict[str, object]] = []
    for gene in union:
        rows.append({
            "record_type": "membership", "panel": "both" if gene in overlap else "literature" if gene in new_genes else "legacy",
            "metric": "panel_membership", "outcome": "NA", "gene_symbol": gene,
            "in_literature_panel": gene in new_genes, "in_legacy_panel": gene in legacy_genes,
            "value": "NA", "n": 1, "notes": "Gene-level membership only; assay-level mapping remains separate.",
        })
    union_count = len(set(union))
    summaries = {
        "literature_unique_genes": len(new_genes), "legacy_unique_genes": len(legacy_genes),
        "overlap_genes": len(overlap), "union_genes": union_count,
        "Jaccard_index": len(overlap) / union_count if union_count else math.nan,
        "overlap_coefficient": len(overlap) / min(len(new_genes), len(legacy_genes)),
    }
    for metric, value in summaries.items():
        rows.append({
            "record_type": "panel_summary", "panel": "new_vs_legacy", "metric": metric, "outcome": "NA",
            "gene_symbol": "NA", "in_literature_panel": "NA", "in_legacy_panel": "NA",
            "value": value, "n": "NA", "notes": "No study count is interpreted as independent replication.",
        })

    mapping_eligible = alpha.get("eligible_for_two_step_mapping", "false").astype(str).str.lower().eq("true")
    new_alpha = alpha[(alpha["availability_status"] == "direct_variant_available") & mapping_eligible]
    alpha_note = f"Direct alpha with eligible target-to-assay mapping; {new_alpha['gene_symbol'].nunique()} unique Olink assay genes."
    legacy_alpha = pd.read_csv(LEGACY_ALPHA)
    rows.append(distribution_record("literature", "alpha_rs429358", "NA", new_alpha.loc[new_alpha["variant"] == "rs429358", "beta"],
                                    alpha_note))
    rows.append(distribution_record("literature", "alpha_rs7412", "NA", new_alpha.loc[new_alpha["variant"] == "rs7412", "beta"],
                                    alpha_note))
    rows.append(distribution_record("legacy", "alpha_rs429358", "NA", legacy_alpha.loc[legacy_alpha["Gene"] != "APOE", "BETA"],
                                    "Legacy alpha table; APOE protein excluded."))
    rows.append(distribution_record("legacy", "alpha_rs7412", "NA", pd.Series(dtype=float),
                                    "Legacy analysis did not extract rs7412 alpha."))

    mapped_genes = set(target_mapping.loc[
        target_mapping["eligible_for_alpha_beta_triangulation"].astype(str).str.lower().eq("true"),
        "gene_symbol",
    ])
    current_primary = beta[
        (beta["method_role"] == "primary")
        & (beta["beta_status"] == "reestimated")
        & beta["gene_symbol"].isin(mapped_genes)
    ]
    legacy_beta = pd.read_csv(LEGACY_BETA)
    legacy_beta = legacy_beta[(legacy_beta["method"] == "Inverse variance weighted") & (legacy_beta["exposure"] != "APOE")]
    legacy_beta["outcome_standard"] = legacy_beta["outcome"].map(LEGACY_OUTCOME_MAP)
    for outcome in OUTCOMES:
        rows.append(distribution_record("literature", "beta_reestimated", outcome,
                                        current_primary.loc[current_primary["outcome"] == outcome, "beta"],
                                        f"Current A1 harmonization and clumping; {current_primary.loc[current_primary['outcome'] == outcome, 'gene_symbol'].nunique()} re-estimated assay genes."))
        rows.append(distribution_record("legacy", "beta_legacy_IVW", outcome,
                                        legacy_beta.loc[legacy_beta["outcome_standard"] == outcome, "b"],
                                        "Legacy 29-protein IVW estimates; pipeline limitations are documented in the audit."))

    legacy_med = pd.read_csv(LEGACY_MEDIATION)
    legacy_med["outcome_standard"] = legacy_med["outcome"].map(LEGACY_OUTCOME_MAP)
    totals = pd.read_csv(TOTALS, sep="\t")
    new_protein_med = mediation[mediation["row_type"] == "protein"]
    for variant in ["rs429358", "rs7412"]:
        for outcome in OUTCOMES:
            total_row = totals[(totals["variant"] == variant) & (totals["outcome"] == outcome)]
            denominator = float(total_row["beta"].iloc[0]) if len(total_row) else math.nan
            new_subset = new_protein_med[(new_protein_med["variant"] == variant) & (new_protein_med["outcome"] == outcome)]
            new_sum = pd.to_numeric(new_subset["indirect_effect"], errors="coerce").sum(min_count=1)
            rows.append({
                "record_type": "total_mediation", "panel": "literature", "metric": "total_mediated_proportion",
                "outcome": outcome, "variant": variant, "gene_symbol": "TOTAL",
                "value": new_sum / denominator if pd.notna(new_sum) else "NA",
                "n": new_subset["gene_symbol"].nunique(),
                "notes": "Enriched-panel estimate among proteins with both alpha and re-estimated beta; not total circulating mediation.",
            })
            if variant == "rs429358":
                legacy_subset = legacy_med[legacy_med["outcome_standard"] == outcome]
                legacy_sum = pd.to_numeric(legacy_subset["mediation"], errors="coerce").sum()
                rows.append({
                    "record_type": "total_mediation", "panel": "legacy", "metric": "total_mediated_proportion_rebased",
                    "outcome": outcome, "variant": variant, "gene_symbol": "TOTAL", "value": legacy_sum / denominator,
                    "n": legacy_subset["gene"].nunique(),
                    "notes": "Legacy indirect effects divided by newly verified current-A1 total effect for denominator comparability; not a rerun of legacy beta.",
                })

    primary_evidence = provenance[provenance["inclusion_status"] == "primary_high_confidence_panel"].copy()
    primary_evidence["PP_H4_num"] = pd.to_numeric(primary_evidence["PP_H4"], errors="coerce")
    subset_genes = {
        "all_primary": set(new_genes),
        "literature_cis_evidence_only": set(primary_evidence.loc[primary_evidence["cis_trans_status"].str.contains("cis", case=False, na=False), "gene_symbol"]),
        "coloc_supported_only": set(primary_evidence.loc[primary_evidence["PP_H4_num"] >= 0.8, "gene_symbol"]),
        "independent_outcome_only": set(primary_evidence.loc[primary_evidence["replication_status"].eq("independent_outcome_replication_corrected_same_direction"), "gene_symbol"]),
        "replication_supported_only": set(primary_evidence.loc[primary_evidence["replication_status"].str.startswith("independent_", na=False), "gene_symbol"]),
    }
    for subset_name, genes in subset_genes.items():
        eligible_subset = set(new_protein_med["gene_symbol"]) & genes
        rows.append({
            "record_type": "sensitivity_subset", "panel": "literature", "metric": subset_name,
            "outcome": "all", "gene_symbol": "NA", "value": len(genes), "n": len(eligible_subset),
            "notes": "value=primary genes meeting evidence filter; n=genes also eligible for mediation.",
        })
        for variant in ["rs429358", "rs7412"]:
            for outcome in OUTCOMES:
                subset = new_protein_med[
                    new_protein_med["gene_symbol"].isin(genes)
                    & new_protein_med["variant"].eq(variant)
                    & new_protein_med["outcome"].eq(outcome)
                ]
                total_row = totals[(totals["variant"] == variant) & (totals["outcome"] == outcome)]
                denominator = float(total_row["beta"].iloc[0]) if len(total_row) else math.nan
                indirect_sum = pd.to_numeric(subset["indirect_effect"], errors="coerce").sum(min_count=1)
                rows.append({
                    "record_type": "sensitivity_total_mediation", "panel": "literature",
                    "metric": subset_name, "outcome": outcome, "variant": variant,
                    "gene_symbol": "TOTAL", "value": indirect_sum / denominator if pd.notna(indirect_sum) else "NA",
                    "n": subset["gene_symbol"].nunique(),
                    "notes": "Main-instrument indirect effects filtered by prespecified literature evidence attribute; empty subsets are NA, not zero.",
                })

    cis_total = cis_mediation[cis_mediation["row_type"] == "total"]
    for variant in ["rs429358", "rs7412"]:
        for outcome in OUTCOMES:
            subset = cis_total[(cis_total["variant"] == variant) & (cis_total["outcome"] == outcome)]
            rows.append({
                "record_type": "sensitivity_total_mediation", "panel": "literature",
                "metric": "cis_only_instrument_set", "outcome": outcome, "variant": variant,
                "gene_symbol": "TOTAL",
                "value": subset["mediated_proportion"].iloc[0] if len(subset) else "NA",
                "n": subset["number_of_eligible_proteins"].iloc[0] if len(subset) else 0,
                "notes": "Two-step mediation rerun using cis-only protein instruments; unavailable results are NA, not zero.",
            })

    strict_protein = strict_mediation[strict_mediation["row_type"] == "protein"]
    strict_total = strict_mediation[strict_mediation["row_type"] == "total"]
    strict_genes = set(strict_protein["gene_symbol"])
    rows.append({
        "record_type": "sensitivity_subset", "panel": "literature", "metric": "strict_annotation_sensitivity",
        "outcome": "all", "gene_symbol": "NA", "value": len(strict_genes), "n": len(strict_genes),
        "notes": "High-confidence target annotation only; all retained alpha and beta estimates use the same UKB-PPP assay.",
    })
    for variant in ["rs429358", "rs7412"]:
        for outcome in OUTCOMES:
            subset = strict_total[(strict_total["variant"] == variant) & (strict_total["outcome"] == outcome)]
            rows.append({
                "record_type": "sensitivity_total_mediation", "panel": "literature",
                "metric": "strict_annotation_sensitivity", "outcome": outcome, "variant": variant,
                "gene_symbol": "TOTAL", "value": subset["mediated_proportion"].iloc[0] if len(subset) else "NA",
                "n": subset["number_of_eligible_proteins"].iloc[0] if len(subset) else 0,
                "notes": "High-confidence mapping sensitivity; moderate mappings are excluded uniformly rather than protein-by-protein.",
            })

    comparison = pd.DataFrame(rows).replace({np.nan: "NA"})
    comparison.to_csv(COMPARISON, sep="\t", index=False)


def build_current_exclusions(
    provenance: pd.DataFrame,
    alpha: pd.DataFrame,
    beta: pd.DataFrame,
    mediation: pd.DataFrame,
    target_mapping: pd.DataFrame,
) -> None:
    genes = sorted(provenance.loc[
        provenance["inclusion_status"] == "primary_high_confidence_panel", "gene_symbol"
    ].unique())
    final_genes = set(mediation.loc[mediation["row_type"] == "protein", "gene_symbol"])
    rows = []
    for gene in genes:
        mapped = target_mapping[
            target_mapping["standardized_gene_symbol"].eq(gene)
            & target_mapping["eligible_for_primary"].astype(str).str.lower().eq("true")
        ]
        alpha_gene = alpha[
            alpha["gene_symbol"].eq(gene) & alpha["availability_status"].eq("direct_variant_available")
        ]
        beta_gene = beta[
            beta["gene_symbol"].eq(gene) & beta["method_role"].eq("primary")
            & beta["analysis_set"].eq("genome_wide_instruments_primary")
            & beta["beta_status"].eq("reestimated")
        ]
        if gene in final_genes:
            continue
        if mapped.empty:
            reasons = sorted(set(
                target_mapping.loc[target_mapping["standardized_gene_symbol"].eq(gene), "exclusion_reason"]
            ) - {"NA", ""})
            reason = ";".join(reasons) or "UKB_PPP_assay_mapping_not_eligible"
        elif set(alpha_gene["variant"]) != {"rs429358", "rs7412"}:
            reason = "one_or_both_direct_APOE_alpha_estimates_unavailable"
        elif beta_gene["outcome"].nunique() < 4:
            reason = "one_or_more_current_A1_outcome_beta_estimates_not_reestimable"
        else:
            reason = "eligibility_assertion_requires_review"
        rows.append({
            "gene_symbol": gene,
            "literature_target_entries": target_mapping.loc[
                target_mapping["standardized_gene_symbol"].eq(gene), "literature_target_entry"
            ].nunique(),
            "candidate_UKB_assays": pd.to_numeric(
                target_mapping.loc[target_mapping["standardized_gene_symbol"].eq(gene), "number_of_candidate_UKB_assays"],
                errors="coerce",
            ).max(),
            "selected_Olink_assay": mapped["UKB_PPP_OID"].iloc[0] if len(mapped) else "NA",
            "mapping_confidence": (
                "high" if any(mapped["mapping_confidence"].eq("high"))
                else "moderate" if len(mapped) else "not_mapped"
            ),
            "rs429358_alpha_available": "rs429358" in set(alpha_gene["variant"]),
            "rs7412_alpha_available": "rs7412" in set(alpha_gene["variant"]),
            "estimable_outcome_beta_count": beta_gene["outcome"].nunique(),
            "eligible_for_final_mediation": False,
            "exclusion_reason": reason,
            "missing_values_treated_as_zero": False,
            "notes": "Published beta was not substituted for a missing same-assay UKB-PPP estimate.",
        })
    output = pd.DataFrame(rows)
    if len(output):
        assert output["exclusion_reason"].notna().all() and output["exclusion_reason"].ne("").all()
    output.to_csv(CURRENT_EXCLUSIONS, sep="\t", index=False)


def main() -> None:
    provenance = pd.read_csv(PROVENANCE, sep="\t", dtype=str).fillna("NA")
    alpha = pd.read_csv(ALPHA, sep="\t", dtype=str).fillna("NA")
    beta = pd.read_csv(BETA, sep="\t", dtype=str).fillna("NA")
    mediation = pd.read_csv(MEDIATION, sep="\t", dtype=str).fillna("NA")
    cis_mediation = pd.read_csv(MEDIATION_CIS, sep="\t", dtype=str).fillna("NA")
    strict_mediation = pd.read_csv(MEDIATION_STRICT, sep="\t", dtype=str).fillna("NA")
    target_mapping = pd.read_csv(TARGET_MAPPING, sep="\t", dtype=str).fillna("NA")
    eligible = set(mediation.loc[mediation["row_type"] == "protein", "gene_symbol"])
    provenance = update_provenance(provenance, alpha, eligible, target_mapping)
    primary_targets = provenance.loc[provenance["inclusion_status"] == "primary_high_confidence_panel",
                                     ["gene_symbol", "protein_name", "protein_form_or_isoform"]].drop_duplicates()
    beta = add_beta_placeholders(beta, primary_targets, target_mapping)
    build_comparison(provenance, alpha, beta, mediation, cis_mediation, strict_mediation, target_mapping)
    build_current_exclusions(provenance, alpha, beta, mediation, target_mapping)
    print(f"Primary target entries: {len(primary_targets)}")
    print(f"Two-step eligible genes: {len(eligible)}")
    print(f"Output: {COMPARISON}")


if __name__ == "__main__":
    main()
