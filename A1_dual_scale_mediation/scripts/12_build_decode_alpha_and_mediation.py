#!/usr/bin/env python3
"""Build auditable deCODE APOE alpha and separate same-platform mediation tables."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
ALPHA_SOURCE = ROOT / "data_processed" / "decode_smp_direct_APOE_alpha.tsv"
LEDGER = ROOT / "logs" / "decode_download_ledger.tsv"
OUTCOMES_CONFIG = PROJECT / "A1_protein_upgrade" / "config" / "outcomes.yml"

ALPHA_OUTPUT = ROOT / "tables" / "APOE_variant_to_decode_somascan_alpha.tsv"
MEDIATION_OUTPUT = ROOT / "tables" / "decode_smp_two_step_mediation.tsv"
MEDIATION_CIS_OUTPUT = ROOT / "tables" / "decode_smp_two_step_mediation_cis_only.tsv"
MEDIATION_CIS_PAV_OUTPUT = ROOT / "tables" / "decode_smp_two_step_mediation_cis_only_PAV_filtered.tsv"
SUMMARY_OUTPUT = ROOT / "tables" / "decode_smp_mediation_summary.tsv"
LOG_OUTPUT = ROOT / "logs" / "decode_smp_mediation.log"

EXPECTED_ALLELES = {
    "rs429358": ("C", "T"),
    "rs7412": ("T", "C"),
}
VARIANT_CONFIG_KEYS = {
    "rs429358": "rs429358_C",
    "rs7412": "rs7412_T",
}
OUTCOME_ORDER = ["AD", "any_AMD", "dry_AMD", "wet_AMD"]
PLANNED_ASSAYS = 9
PLANNED_VARIANTS = 2
PLANNED_OUTCOMES = 4
PLANNED_MEDIATION_FAMILY = PLANNED_ASSAYS * PLANNED_VARIANTS * PLANNED_OUTCOMES
BOOTSTRAP_DRAWS = 100_000
SEED = 20260719


def bh_adjust(values: pd.Series, planned_n: int) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().sort_values()
    if valid.empty:
        return result
    ranks = np.arange(1, len(valid) + 1)
    adjusted = valid.to_numpy() * planned_n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[valid.index] = np.minimum(1.0, adjusted)
    return result


def p_display(p: float, minus_log10_p: float) -> str:
    if not np.isfinite(minus_log10_p):
        return "NA"
    if minus_log10_p > 300:
        return "<1e-300"
    return f"{p:.3e}"


def build_alpha() -> pd.DataFrame:
    alpha = pd.read_csv(ALPHA_SOURCE, sep="\t", dtype={"requested_variant": str})
    ledger = pd.read_csv(LEDGER, sep="\t", dtype=str)
    ledger = ledger[ledger["normalization"].eq("smp")].rename(columns={
        "SomaScan_SeqId": "assay_target_ID",
        "object_key": "source_object_key",
        "sha256": "source_file_sha256",
        "etag": "source_object_etag",
    })
    alpha = alpha.rename(columns={
        "SomaScan_SeqId": "assay_target_ID",
        "requested_variant": "variant",
        "Beta": "alpha",
        "SE": "alpha_SE",
        "Pval": "alpha_P",
        "N": "alpha_N",
        "ImpMAF": "imputation_MAF",
        "effectAllele": "effect_allele",
        "otherAllele": "other_allele",
        "Chrom": "chromosome_hg38",
        "Pos": "position_hg38",
    })
    alpha = alpha.merge(
        ledger[["assay_target_ID", "source_object_key", "source_file_sha256", "source_object_etag"]],
        on="assay_target_ID", how="left", validate="many_to_one",
    )
    if len(alpha) != PLANNED_ASSAYS * PLANNED_VARIANTS:
        raise ValueError(f"Expected 18 direct alpha rows, found {len(alpha)}")
    if alpha.duplicated(["assay_target_ID", "variant"]).any():
        raise ValueError("Duplicate assay-variant alpha rows")
    for variant, (effect, other) in EXPECTED_ALLELES.items():
        subset = alpha[alpha["variant"].eq(variant)]
        if not (subset["effect_allele"].eq(effect) & subset["other_allele"].eq(other)).all():
            raise ValueError(f"Allele orientation mismatch for {variant}")

    alpha["original_effect_allele"] = alpha["effect_allele"]
    alpha["original_other_allele"] = alpha["other_allele"]
    alpha["harmonized_effect_allele"] = alpha["effect_allele"]
    alpha["harmonized_other_allele"] = alpha["other_allele"]
    alpha["allele_flipped"] = False
    alpha["direct_variant_not_proxy"] = True
    alpha["effect_unit"] = "SomaScan_SMP_normalized_protein_per_effect_allele"
    alpha["imputation_MAF_not_EAF"] = True
    alpha["alpha_F_statistic"] = (alpha["alpha"] / alpha["alpha_SE"]) ** 2
    alpha["alpha_P_FDR_18"] = bh_adjust(alpha["alpha_P"], len(alpha))
    alpha["alpha_P_Bonferroni_18"] = np.minimum(1.0, alpha["alpha_P"] * len(alpha))
    alpha["alpha_significance_not_used_for_eligibility"] = True
    alpha["manual_verification_status"] = "direct_variant_alleles_and_full_file_row_verified"
    alpha["alpha_source"] = "deCODE_SomaScan_SMP_full_summary_statistics"

    columns = [
        "gene_symbol", "assay_target_ID", "UniProt_ID", "variant",
        "original_effect_allele", "original_other_allele", "harmonized_effect_allele",
        "harmonized_other_allele", "allele_flipped", "direct_variant_not_proxy",
        "chromosome_hg38", "position_hg38", "alpha", "alpha_SE", "alpha_P", "alpha_N",
        "imputation_MAF", "imputation_MAF_not_EAF", "effect_unit", "alpha_F_statistic",
        "alpha_P_FDR_18", "alpha_P_Bonferroni_18", "alpha_significance_not_used_for_eligibility",
        "alpha_source", "source_object_key", "source_object_etag", "source_file_sha256",
        "manual_verification_status",
    ]
    return alpha[columns].sort_values(["assay_target_ID", "variant"]).reset_index(drop=True)


def run_mediation(alpha: pd.DataFrame, beta_path: Path, scope: str, output: Path, total_effects: dict) -> pd.DataFrame:
    beta = pd.read_csv(beta_path, sep="\t")
    beta = beta[beta["method_role"].eq("primary")].copy()
    if beta.duplicated(["assay_target_ID", "outcome"]).any():
        raise ValueError(f"Duplicate primary beta rows in {beta_path}")

    grid = pd.MultiIndex.from_product(
        [sorted(alpha["assay_target_ID"].unique()), EXPECTED_ALLELES.keys(), OUTCOME_ORDER],
        names=["assay_target_ID", "variant", "outcome"],
    ).to_frame(index=False)
    alpha_join = alpha.drop(columns=["manual_verification_status"]).copy()
    merged = grid.merge(alpha_join, on=["assay_target_ID", "variant"], how="left", validate="many_to_one")
    beta_columns = [
        "gene_symbol", "assay_target_ID", "UniProt_ID", "outcome", "method", "nsnp", "n_clumped",
        "beta", "SE", "P_value", "P_value_display", "P_Bonferroni_within_outcome",
        "P_Bonferroni_global_36", "Cochran_Q", "Q_df", "Q_P", "Egger_intercept",
        "Egger_intercept_P", "mean_F", "min_F", "beta_status", "exclusion_reason", "instrument_scope",
    ]
    beta = beta[[column for column in beta_columns if column in beta.columns]].rename(columns={
        "gene_symbol": "beta_gene_symbol",
        "UniProt_ID": "beta_UniProt_ID",
        "beta": "protein_outcome_beta",
        "SE": "protein_outcome_SE",
        "P_value": "protein_outcome_P",
        "P_value_display": "protein_outcome_P_display",
        "beta_status": "protein_outcome_beta_status",
        "exclusion_reason": "protein_outcome_exclusion_reason",
        "mean_F": "protein_outcome_mean_F",
        "min_F": "protein_outcome_min_F",
    })
    merged = merged.merge(beta, on=["assay_target_ID", "outcome"], how="left", validate="many_to_one")
    if not (merged["gene_symbol"].fillna(merged["beta_gene_symbol"]) == merged["beta_gene_symbol"].fillna(merged["gene_symbol"])).all():
        raise ValueError("Gene mismatch between alpha and beta tables")
    merged["analysis_scope"] = scope
    merged["planned_mediation_family"] = PLANNED_MEDIATION_FAMILY
    merged["total_effect"] = [
        total_effects[VARIANT_CONFIG_KEYS[variant]][outcome]["beta"]
        for variant, outcome in zip(merged["variant"], merged["outcome"])
    ]
    merged["total_effect_source"] = "frozen_current_A1_harmonized_GWAS"
    merged["mediation_status"] = np.where(
        merged["protein_outcome_beta_status"].eq("reestimated")
        & merged["alpha"].notna() & merged["protein_outcome_beta"].notna(),
        "estimable", "not_estimable",
    )

    scope_seed = {"genome_wide_APOE_excluded": SEED, "cis_only": SEED + 1,
                  "cis_only_target_gene_PAV_filtered": SEED + 2}[scope]
    rng = np.random.default_rng(scope_seed)
    rows = []
    for _, row in merged.iterrows():
        record = row.to_dict()
        if record["mediation_status"] != "estimable":
            record.update({
                "indirect_effect": np.nan, "indirect_effect_delta_SE": np.nan,
                "indirect_effect_delta_CI_lower": np.nan, "indirect_effect_delta_CI_upper": np.nan,
                "indirect_effect_delta_P": np.nan, "indirect_effect_delta_minus_log10_P": np.nan,
                "indirect_effect_delta_P_display": "NA", "indirect_effect_bootstrap_CI_lower": np.nan,
                "indirect_effect_bootstrap_CI_upper": np.nan, "mediated_proportion": np.nan,
                "mediated_proportion_bootstrap_CI_lower": np.nan,
                "mediated_proportion_bootstrap_CI_upper": np.nan,
                "direction_classification": "not_estimable",
            })
            rows.append(record)
            continue
        alpha_value = float(record["alpha"])
        alpha_se = float(record["alpha_SE"])
        beta_value = float(record["protein_outcome_beta"])
        beta_se = float(record["protein_outcome_SE"])
        total = float(record["total_effect"])
        indirect = alpha_value * beta_value
        delta_se = math.sqrt((beta_value * alpha_se) ** 2 + (alpha_value * beta_se) ** 2)
        z_value = indirect / delta_se if delta_se > 0 else np.nan
        log_p = math.log(2) + norm.logsf(abs(z_value)) if np.isfinite(z_value) else np.nan
        p_value = math.exp(log_p) if np.isfinite(log_p) and log_p > math.log(np.finfo(float).tiny) else 0.0
        minus_log10 = -log_p / math.log(10) if np.isfinite(log_p) else np.nan
        alpha_draws = rng.normal(alpha_value, alpha_se, BOOTSTRAP_DRAWS)
        beta_draws = rng.normal(beta_value, beta_se, BOOTSTRAP_DRAWS)
        indirect_draws = alpha_draws * beta_draws
        proportion_draws = indirect_draws / total
        direction = "concordant_mediation" if np.sign(indirect) == np.sign(total) else "opposing_or_suppressing"
        record.update({
            "indirect_effect": indirect,
            "indirect_effect_delta_SE": delta_se,
            "indirect_effect_delta_CI_lower": indirect - 1.96 * delta_se,
            "indirect_effect_delta_CI_upper": indirect + 1.96 * delta_se,
            "indirect_effect_delta_P": p_value,
            "indirect_effect_delta_minus_log10_P": minus_log10,
            "indirect_effect_delta_P_display": p_display(p_value, minus_log10),
            "indirect_effect_bootstrap_CI_lower": np.quantile(indirect_draws, 0.025),
            "indirect_effect_bootstrap_CI_upper": np.quantile(indirect_draws, 0.975),
            "mediated_proportion": indirect / total,
            "mediated_proportion_bootstrap_CI_lower": np.quantile(proportion_draws, 0.025),
            "mediated_proportion_bootstrap_CI_upper": np.quantile(proportion_draws, 0.975),
            "direction_classification": direction,
        })
        rows.append(record)

    result = pd.DataFrame(rows)
    estimable = result["mediation_status"].eq("estimable")
    result["mediation_P_FDR_72"] = np.nan
    result.loc[estimable, "mediation_P_FDR_72"] = bh_adjust(
        result.loc[estimable, "indirect_effect_delta_P"], PLANNED_MEDIATION_FAMILY
    )
    result["mediation_P_Bonferroni_72"] = np.where(
        estimable, np.minimum(1.0, result["indirect_effect_delta_P"] * PLANNED_MEDIATION_FAMILY), np.nan
    )
    result["passes_Bonferroni_72"] = result["mediation_P_Bonferroni_72"].lt(0.05)
    result["bootstrap_draws"] = BOOTSTRAP_DRAWS
    result["bootstrap_seed"] = scope_seed
    result["bootstrap_covariance_model"] = "independent_normal_alpha_beta_total_effect_fixed"
    result["bootstrap_limitation"] = "does_not_model_alpha_beta_covariance_mapping_uncertainty_or_total_effect_uncertainty"
    result["interpretation_boundary"] = "same_platform_non_UKB_sensitivity_not_merged_into_frozen_25_assay_aggregate"
    result.to_csv(output, sep="\t", index=False, na_rep="NA")
    return result


def main() -> None:
    for required in (ALPHA_SOURCE, LEDGER, OUTCOMES_CONFIG,
                     ROOT / "tables" / "decode_smp_beta_results.tsv",
                     ROOT / "tables" / "decode_smp_beta_results_cis_only.tsv",
                     ROOT / "tables" / "decode_smp_beta_results_cis_only_PAV_filtered.tsv"):
        if not required.exists():
            raise FileNotFoundError(required)
    with OUTCOMES_CONFIG.open("rt", encoding="utf-8") as handle:
        total_effects = yaml.safe_load(handle)["apoe_total_effects"]

    alpha = build_alpha()
    alpha.to_csv(ALPHA_OUTPUT, sep="\t", index=False, na_rep="NA")
    genome = run_mediation(alpha, ROOT / "tables" / "decode_smp_beta_results.tsv",
                           "genome_wide_APOE_excluded", MEDIATION_OUTPUT, total_effects)
    cis = run_mediation(alpha, ROOT / "tables" / "decode_smp_beta_results_cis_only.tsv",
                        "cis_only", MEDIATION_CIS_OUTPUT, total_effects)
    cis_pav = run_mediation(alpha, ROOT / "tables" / "decode_smp_beta_results_cis_only_PAV_filtered.tsv",
                            "cis_only_target_gene_PAV_filtered", MEDIATION_CIS_PAV_OUTPUT, total_effects)

    summaries = []
    for label, frame in [("genome_wide_APOE_excluded", genome), ("cis_only", cis),
                         ("cis_only_target_gene_PAV_filtered", cis_pav)]:
        estimable = frame[frame["mediation_status"].eq("estimable")]
        survivors = estimable[estimable["passes_Bonferroni_72"]]
        summaries.append({
            "analysis_scope": label,
            "planned_paths": PLANNED_MEDIATION_FAMILY,
            "estimable_paths": len(estimable),
            "not_estimable_paths": PLANNED_MEDIATION_FAMILY - len(estimable),
            "bonferroni_surviving_paths": len(survivors),
            "surviving_path_ids": ";".join(
                f"{r.assay_target_ID}:{r.variant}:{r.outcome}:{r.direction_classification}"
                for r in survivors.itertuples()
            ) or "none",
            "aggregate_computed": False,
            "aggregate_reason": "independent_platform_extension_not_merged_or_summed",
        })
    pd.DataFrame(summaries).to_csv(SUMMARY_OUTPUT, sep="\t", index=False)
    LOG_OUTPUT.write_text(
        "\n".join([
            f"alpha_rows={len(alpha)}",
            f"planned_mediation_paths_per_scope={PLANNED_MEDIATION_FAMILY}",
            f"bootstrap_draws={BOOTSTRAP_DRAWS}",
            f"seed={SEED}",
            "alpha_significance_not_used_for_eligibility=true",
            "aggregate_computed=false",
            "bootstrap_covariance_and_mapping_uncertainty_not_modelled=true",
        ]) + "\n", encoding="utf-8"
    )
    print(pd.DataFrame(summaries).to_string(index=False))
    print(f"alpha_output={ALPHA_OUTPUT}")


if __name__ == "__main__":
    main()
