#!/usr/bin/env python3
"""Hard assertions for the completed deCODE same-platform extension."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "audit" / "decode_extension_validation_report.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    ledger = pd.read_csv(ROOT / "logs" / "decode_download_ledger.tsv", sep="\t")
    require(len(ledger[ledger.normalization.eq("smp")]) == 9, "Expected nine completed SMP objects")
    require(len(ledger[ledger.normalization.eq("raw")]) == 2, "Expected two gated raw objects")
    require(ledger.status.isin(["downloaded_complete", "already_complete", "resumed_complete"]).all(),
            "Every ledger object must be complete")
    require(not list((ROOT / "data_raw" / "decode_somascan").rglob("*.part")), "No .part file may remain")

    alpha = pd.read_csv(ROOT / "tables" / "APOE_variant_to_decode_somascan_alpha.tsv", sep="\t")
    require(len(alpha) == 18 and not alpha.duplicated(["assay_target_ID", "variant"]).any(),
            "Expected 18 unique direct SMP alpha rows")
    require(alpha.direct_variant_not_proxy.eq(True).all(), "All alpha rows must be direct variants")

    beta = pd.read_csv(ROOT / "tables" / "decode_smp_beta_results.tsv", sep="\t")
    primary_beta = beta[(beta.method_role.eq("primary")) & beta.beta_status.eq("reestimated")]
    require(len(primary_beta) == 36, "Expected 36 primary SMP beta estimates")
    require(primary_beta.assay_target_ID.nunique() == 9, "Expected nine SMP assays")
    require(primary_beta.outcome.nunique() == 4, "Expected four frozen outcomes")
    zero_display = primary_beta.loc[primary_beta.P_value.eq(0), "P_value_display"].astype(str)
    require((zero_display == "<1e-300").all(), "Underflowed P values require a non-zero display bound")

    expected_survivors = {
        "decode_smp_two_step_mediation.tsv": 3,
        "decode_smp_two_step_mediation_cis_only.tsv": 1,
        "decode_smp_two_step_mediation_cis_only_PAV_filtered.tsv": 0,
    }
    for filename, expected in expected_survivors.items():
        frame = pd.read_csv(ROOT / "tables" / filename, sep="\t")
        require(len(frame) == 72, f"{filename} must retain the planned 72-row grid")
        require(int(frame.passes_Bonferroni_72.eq(True).sum()) == expected,
                f"Unexpected survivor count in {filename}")
        require(frame.interpretation_boundary.str.contains("not_merged", na=False).all(),
                f"{filename} must preserve the no-merge boundary")

    raw_summary = pd.read_csv(ROOT / "tables" / "decode_raw_gated_sensitivity_summary.tsv", sep="\t")
    require(len(raw_summary) == 3, "Expected three gated raw scopes")
    require(raw_summary.Bonferroni_16_survivors.eq(0).all(), "No gated raw path should survive")
    require(raw_summary.independent_confirmation.eq(False).all(), "Raw sensitivity cannot be independent confirmation")

    gate = pd.read_csv(ROOT / "tables" / "decode_same_platform_feasibility_gate.tsv", sep="\t")
    require(len(gate) == 9, "Expected nine assay-level gate rows")
    require(gate.current_gate_status.eq("same_platform_alpha_beta_mediation_complete").all(),
            "Every gate row must be complete")
    require(gate.reestimated_outcome_beta_count.eq(4).all(), "Each assay must have four beta estimates")

    attrition = pd.read_csv(ROOT / "tables" / "TableS28_V2_candidate_attrition_provenance.tsv", sep="\t")
    require(len(attrition) == 8, "Expected eight excluded unique genes")
    require(attrition.same_platform_alpha_beta_loop_status.eq("complete_separate_non_UKB_sensitivity").all(),
            "Every excluded gene must have a completed deCODE loop")
    require(attrition.deCODE_reestimated_beta_count.sum() == 36, "Gene-level beta counts must sum to 36")

    report = "\n".join([
        "# deCODE extension validation",
        "",
        "- Status: PASS",
        "- Download ledger: 9 SMP + 2 gated raw objects; no partial file remains.",
        "- Direct alpha: 18 unique assay-variant rows.",
        "- Primary beta: 36 assay-outcome estimates.",
        "- SMP mediation survivors: genome-wide 3; cis-only 1; target-gene PAV-filtered 0.",
        "- Gated raw survivors: 0 in all three scopes; not treated as independent confirmation.",
        "- Final gate: 9/9 assays complete; Table S28 contains 8/8 completed separate sensitivity loops.",
    ]) + "\n"
    REPORT.write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
