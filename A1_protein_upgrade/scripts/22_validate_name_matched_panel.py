#!/usr/bin/env python3
"""Validate dynamic counts, assay identity, corrections, and figure source data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
FIGURE_SOURCE = ROOT / "figures" / "source_data"
OUTPUT = ROOT / "logs" / "name_match_revision_QA.tsv"


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, sep="\t", dtype=str).fillna("NA")


def bh(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def check(results: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    results.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> None:
    mapping = read("APOE_linkable_target_assay_mapping.tsv")
    alpha = read("APOE_variant_to_literature_proteins_alpha.tsv")
    beta = read("literature_panel_beta_results.tsv")
    med = read("APOE_linkable_two_step_mediation.tsv")
    med_cis = read("APOE_linkable_two_step_mediation_cis_sensitivity.tsv")
    med_strict = read("APOE_linkable_two_step_mediation_strict_sensitivity.tsv")
    exclusions = read("name_matched_panel_excluded_genes.tsv")
    flow = read("APOE_linkable_subset_flow.tsv")
    results: list[dict[str, object]] = []

    target_n = len(mapping)
    gene_n = mapping["standardized_gene_symbol"].nunique()
    eligible_map = mapping[mapping["eligible_for_primary"].str.lower().eq("true")]
    matched_assays = eligible_map[["standardized_gene_symbol", "UKB_PPP_OID"]].drop_duplicates()
    strict_assays = mapping[mapping["eligible_for_strict_sensitivity"].str.lower().eq("true")][
        ["standardized_gene_symbol", "UKB_PPP_OID"]
    ].drop_duplicates()
    check(results, "target_entry_count", target_n == 41, f"observed={target_n}; expected=41")
    check(results, "unique_gene_count", gene_n == 33, f"observed={gene_n}; expected=33")
    check(
        results, "included_mapping_confidence",
        set(eligible_map["mapping_confidence"]).issubset({"high", "moderate"}),
        f"values={sorted(set(eligible_map['mapping_confidence']))}",
    )
    excluded_map = mapping[~mapping["eligible_for_primary"].str.lower().eq("true")]
    check(
        results, "mapping_exclusion_reasons",
        excluded_map["exclusion_reason"].notna().all() and excluded_map["exclusion_reason"].ne("NA").all(),
        f"excluded_target_entries={len(excluded_map)}",
    )

    alpha_available = alpha[alpha["availability_status"].eq("direct_variant_available")]
    check(
        results, "alpha_assay_unique_per_gene",
        alpha_available.groupby("gene_symbol")["assay_target_ID"].nunique().max() == 1,
        f"genes={alpha_available['gene_symbol'].nunique()}",
    )

    for label, frame in [("main", med), ("cis", med_cis), ("strict", med_strict)]:
        proteins = frame[frame["row_type"].eq("protein")].copy()
        final_n = proteins["gene_symbol"].nunique()
        expected_paths = final_n * 2 * 4
        check(results, f"{label}_path_count", len(proteins) == expected_paths,
              f"genes={final_n}; observed={len(proteins)}; expected={expected_paths}")
        check(results, f"{label}_alpha_beta_present",
              proteins[["alpha", "SE_alpha", "beta", "SE_beta"]].replace("NA", np.nan).notna().all().all(),
              f"paths={len(proteins)}")
        check(results, f"{label}_same_assay",
              proteins["alpha_assay_target_ID"].eq(proteins["beta_assay_target_ID"]).all(),
              f"paths={len(proteins)}")
        family = pd.to_numeric(proteins["mediation_family_size"], errors="coerce")
        check(results, f"{label}_mediation_family_size",
              family.eq(expected_paths).all(), f"family={sorted(family.dropna().unique().tolist())}")

    final_main = med[med["row_type"].eq("protein")]
    flow_final = int(flow.loc[flow["stage"].eq("final_mediation_genes"), "n"].iloc[0])
    flow_paths = int(flow.loc[flow["stage"].eq("mediation_path_count"), "n"].iloc[0])
    check(results, "flow_matches_main", flow_final == final_main["gene_symbol"].nunique() and flow_paths == len(final_main),
          f"flow_genes={flow_final}; flow_paths={flow_paths}")

    primary_beta = beta[
        beta["method_role"].eq("primary")
        & beta["beta_status"].isin(["reestimated", "reestimated_cis_only"])
    ].copy()
    for (analysis_set, outcome), group in primary_beta.groupby(["analysis_set", "outcome"]):
        observed = pd.to_numeric(group["P_FDR_observed_reestimated"], errors="coerce").to_numpy(float)
        expected = bh(group["P_value"])
        check(results, f"FDR_{analysis_set}_{outcome}", np.allclose(observed, expected, rtol=1e-9, atol=1e-12),
              f"n={len(group)}")
        planned = len(matched_assays)
        bonf_observed = pd.to_numeric(group["P_Bonferroni_planned_family"], errors="coerce").to_numpy(float)
        bonf_expected = np.minimum(1.0, pd.to_numeric(group["P_value"], errors="coerce").to_numpy(float) * planned)
        family_observed = set(pd.to_numeric(group["planned_family_size"], errors="coerce").dropna().astype(int))
        check(results, f"Bonferroni_{analysis_set}_{outcome}",
              np.allclose(bonf_observed, bonf_expected, rtol=1e-9, atol=1e-12) and family_observed == {planned},
              f"planned_assays={planned}; recorded={sorted(family_observed)}")

    check(results, "current_exclusion_reasons",
          exclusions["exclusion_reason"].notna().all() and exclusions["exclusion_reason"].ne("NA").all(),
          f"excluded_genes={len(exclusions)}")
    check(results, "strict_subset_of_primary",
          set(strict_assays["standardized_gene_symbol"]).issubset(set(matched_assays["standardized_gene_symbol"])),
          f"strict={len(strict_assays)}; primary={len(matched_assays)}")

    figure_files = [
        "Figure_5b_APOE_linkage_scatter.csv", "Figure_5c_protein_mediation_forest.csv",
        "Figure_5d_panel_boundary_comparison.csv",
    ]
    for filename in figure_files:
        check(results, f"figure_source_{filename}", (FIGURE_SOURCE / filename).exists(), filename)

    output = pd.DataFrame(results)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)
    failures = output[output["status"].eq("FAIL")]
    print(output.to_string(index=False))
    if len(failures):
        raise SystemExit(f"QA failed: {len(failures)} checks")


if __name__ == "__main__":
    main()
