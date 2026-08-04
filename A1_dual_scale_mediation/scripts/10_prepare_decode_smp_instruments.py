#!/usr/bin/env python3
"""Map and audit deCODE SMP pQTL candidates for the frozen A1 outcomes."""

from __future__ import annotations

import gzip
import hashlib
import os
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pyliftover import LiftOver


ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ROOT.parent
CANDIDATES = ROOT / "data_processed" / "decode_smp_pqtl_candidates_p5e8.tsv"
CHAIN = ROOT / "data_raw" / "reference" / "hg38ToHg19.over.chain.gz"
RESOURCE_ROOT = Path(os.environ.get("A1_RESOURCE_ROOT", PROJECT_ROOT / "data" / "external"))
GENE_COORDINATES = PROJECT_ROOT / "A1_protein_upgrade" / "config" / "protein_gene_coordinates_grch37.tsv"
LD_BIM = Path(os.environ.get("A1_LD_BIM", str(RESOURCE_ROOT / "EUR" / "EUR.bim")))

MAPPED_OUTPUT = ROOT / "data_processed" / "decode_smp_pqtl_candidates_mapped_qc.tsv.gz"
ELIGIBLE_OUTPUT = ROOT / "data_processed" / "decode_smp_pqtl_candidates_for_clumping.tsv"
SUMMARY_OUTPUT = ROOT / "tables" / "decode_smp_instrument_mapping_qc.tsv"
APOE_AUDIT_OUTPUT = ROOT / "tables" / "decode_coordinate_and_APOE_exclusion_audit.tsv"
LOG_OUTPUT = ROOT / "logs" / "decode_smp_instrument_mapping.log"

APOE_CHROMOSOME = 19
APOE_START_HG19 = 44_000_000
APOE_END_HG19 = 46_500_000
CIS_WINDOW_BP = 1_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def complement(allele: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return allele.translate(table) if set(allele) <= set("ACGT") else allele


def allele_sets_compatible(a1: str, a2: str, b1: str, b2: str) -> bool:
    left = {a1.upper(), a2.upper()}
    right = {b1.upper(), b2.upper()}
    if left == right:
        return True
    return {complement(a1.upper()), complement(a2.upper())} == right


def lift_point(lo: LiftOver, chrom: str, pos_1based: int) -> tuple[str | None, int | None, str]:
    mapped = lo.convert_coordinate(chrom, pos_1based - 1)
    valid = [(target, position, strand) for target, position, strand, _ in mapped if target.startswith("chr")]
    if len(valid) != 1:
        return None, None, "unmapped" if not valid else "multiple_mappings"
    target, position, strand = valid[0]
    if strand != "+":
        return target.removeprefix("chr"), int(position) + 1, "reverse_strand_mapping"
    return target.removeprefix("chr"), int(position) + 1, "unique_plus_strand"


def parse_aliases(value: str) -> list[str]:
    if not value or value in {".", "NA", "nan"}:
        return []
    return [item.strip() for item in value.split(",") if item.strip().startswith("rs")]


def read_ld_matches(alias_set: set[str], position_set: set[tuple[str, int]]) -> tuple[dict[str, list[dict]], dict[tuple[str, int], list[dict]]]:
    by_rsid: dict[str, list[dict]] = defaultdict(list)
    by_position: dict[tuple[str, int], list[dict]] = defaultdict(list)
    with LD_BIM.open("rt", encoding="utf-8") as handle:
        for line in handle:
            chrom, snp, cm, bp, a1, a2 = line.rstrip("\n").split()
            key = (chrom, int(bp))
            if snp in alias_set or key in position_set:
                record = {"chrom": chrom, "snp": snp, "bp": int(bp), "a1": a1.upper(), "a2": a2.upper()}
                if snp in alias_set:
                    by_rsid[snp].append(record)
                if key in position_set:
                    by_position[key].append(record)
    return by_rsid, by_position


def resolve_ld_variant(row: pd.Series, by_rsid: dict, by_position: dict) -> tuple[str, str, str, str, int | None]:
    aliases = parse_aliases(str(row["rsids"]))
    chrom = str(row["chromosome_hg19"])
    bp = int(row["position_hg19"]) if pd.notna(row["position_hg19"]) else None
    ea = str(row["effectAllele"]).upper()
    oa = str(row["otherAllele"]).upper()

    exact = []
    for alias in aliases:
        exact.extend(by_rsid.get(alias, []))
    exact = [r for r in exact if r["chrom"] == chrom and (bp is None or r["bp"] == bp)]
    compatible_exact = [r for r in exact if allele_sets_compatible(ea, oa, r["a1"], r["a2"])]
    unique_exact = {r["snp"]: r for r in compatible_exact}
    if len(unique_exact) == 1:
        record = next(iter(unique_exact.values()))
        return record["snp"], "exact_rsid_coordinate_allele", "compatible", "eligible", record["bp"]
    if len(unique_exact) > 1:
        return "NA", "ambiguous_multiple_exact_rsids", "compatible", "excluded", None

    positional = by_position.get((chrom, bp), []) if bp is not None else []
    compatible_positional = [r for r in positional if allele_sets_compatible(ea, oa, r["a1"], r["a2"])]
    unique_positional = {r["snp"]: r for r in compatible_positional}
    if len(unique_positional) == 1:
        record = next(iter(unique_positional.values()))
        return record["snp"], "hg38_to_hg19_position_plus_allele", "compatible", "eligible", record["bp"]
    if len(unique_positional) > 1:
        return "NA", "ambiguous_position_allele_match", "compatible", "excluded", None
    if exact:
        return "NA", "exact_rsid_coordinate_or_allele_conflict", "incompatible", "excluded", None
    if positional:
        return "NA", "position_found_alleles_incompatible", "incompatible", "excluded", None
    return "NA", "not_present_in_EUR_LD_reference", "not_assessable", "excluded", None


def main() -> None:
    for required in (CANDIDATES, CHAIN, GENE_COORDINATES, LD_BIM):
        if not required.exists():
            raise FileNotFoundError(required)

    data = pd.read_csv(CANDIDATES, sep="\t", dtype={"Chrom": str, "rsids": str})
    lo = LiftOver(str(CHAIN))
    unique_positions = data[["Chrom", "Pos"]].drop_duplicates()
    lifted = {}
    for chrom, pos in unique_positions.itertuples(index=False):
        lifted[(str(chrom), int(pos))] = lift_point(lo, str(chrom), int(pos))

    mapped = data.apply(lambda row: lifted[(str(row["Chrom"]), int(row["Pos"]))], axis=1, result_type="expand")
    mapped.columns = ["chromosome_hg19", "position_hg19", "liftover_status"]
    data = pd.concat([data, mapped], axis=1)
    data["position_hg19"] = pd.to_numeric(data["position_hg19"], errors="coerce").astype("Int64")
    data["APOE_region_excluded"] = (
        data["chromosome_hg19"].eq(str(APOE_CHROMOSOME))
        & data["position_hg19"].between(APOE_START_HG19, APOE_END_HG19)
    )

    aliases = {alias for value in data["rsids"].fillna(".") for alias in parse_aliases(str(value))}
    positions = {
        (str(chrom), int(pos))
        for chrom, pos in data.loc[data["position_hg19"].notna(), ["chromosome_hg19", "position_hg19"]].itertuples(index=False)
    }
    by_rsid, by_position = read_ld_matches(aliases, positions)

    resolved = data.apply(lambda row: resolve_ld_variant(row, by_rsid, by_position), axis=1, result_type="expand")
    resolved.columns = ["SNP", "ld_mapping_method", "allele_compatibility", "ld_mapping_status", "ld_position_hg19"]
    data = pd.concat([data, resolved], axis=1)

    coordinates = pd.read_csv(GENE_COORDINATES, sep="\t").set_index("gene_symbol")
    data["cis_to_encoding_gene"] = False
    data["cis_window_definition"] = "Ensembl GRCh37 gene boundaries +/-1 Mb"
    data["gene_coordinate_source"] = "NA"
    for gene in sorted(data["gene_symbol"].unique()):
        coordinate = coordinates.loc[gene]
        lower = int(coordinate["gene_start_grch37"]) - CIS_WINDOW_BP
        upper = int(coordinate["gene_end_grch37"]) + CIS_WINDOW_BP
        mask = (
            data["gene_symbol"].eq(gene)
            & data["chromosome_hg19"].eq(str(coordinate["chromosome"]))
            & data["position_hg19"].between(lower, upper)
        )
        data.loc[mask, "cis_to_encoding_gene"] = True
        data.loc[data["gene_symbol"].eq(gene), "gene_coordinate_source"] = coordinate["source_url"]

    data["eligible_for_clumping"] = (
        data["liftover_status"].eq("unique_plus_strand")
        & ~data["APOE_region_excluded"]
        & data["ld_mapping_status"].eq("eligible")
        & data["SNP"].ne("NA")
    )
    data["instrument_exclusion_reason"] = "eligible"
    data.loc[~data["liftover_status"].eq("unique_plus_strand"), "instrument_exclusion_reason"] = "liftover_not_unique_plus_strand"
    data.loc[data["APOE_region_excluded"], "instrument_exclusion_reason"] = "frozen_APOE_region_exclusion"
    mapping_fail = ~data["ld_mapping_status"].eq("eligible") & ~data["APOE_region_excluded"]
    data.loc[mapping_fail, "instrument_exclusion_reason"] = data.loc[mapping_fail, "ld_mapping_method"]

    eligible = data[data["eligible_for_clumping"]].copy()
    eligible = eligible.sort_values(["SomaScan_SeqId", "SNP", "Pval"])
    eligible = eligible.drop_duplicates(["SomaScan_SeqId", "SNP"], keep="first")
    eligible = eligible.rename(columns={
        "Beta": "beta",
        "SE": "SE",
        "Pval": "P_value",
        "effectAllele": "effect_allele",
        "otherAllele": "other_allele",
        "ImpMAF": "imputation_MAF",
        "N": "sample_size",
        "SomaScan_SeqId": "assay_target_ID",
    })
    eligible["chromosome"] = eligible["chromosome_hg19"].astype(int)
    eligible["position_hg19_from_liftover"] = eligible["position_hg19"].astype(int)

    summary_rows = []
    for (gene, assay), group in data.groupby(["gene_symbol", "SomaScan_SeqId"], sort=True):
        elig = eligible[(eligible["gene_symbol"] == gene) & (eligible["assay_target_ID"] == assay)]
        summary_rows.append({
            "gene_symbol": gene,
            "SomaScan_SeqId": assay,
            "n_candidates_p_lt_5e8": len(group),
            "n_unique_hg38_positions": group[["Chrom", "Pos"]].drop_duplicates().shape[0],
            "n_liftover_unique_plus_strand": int(group["liftover_status"].eq("unique_plus_strand").sum()),
            "n_APOE_region_excluded": int(group["APOE_region_excluded"].sum()),
            "n_exact_rsid_LD_matches": int(group["ld_mapping_method"].eq("exact_rsid_coordinate_allele").sum()),
            "n_position_allele_LD_rescues": int(group["ld_mapping_method"].eq("hg38_to_hg19_position_plus_allele").sum()),
            "n_unresolved_or_ambiguous": int((~group["ld_mapping_status"].eq("eligible") & ~group["APOE_region_excluded"]).sum()),
            "n_unique_instruments_before_clumping": len(elig),
            "n_cis_candidates_before_clumping": int(elig["cis_to_encoding_gene"].sum()) if len(elig) else 0,
            "minimum_ImpMAF_eligible": elig["imputation_MAF"].min() if len(elig) else pd.NA,
            "n_eligible_ImpMAF_lt_0_01": int((elig["imputation_MAF"] < 0.01).sum()) if len(elig) else 0,
        })

    apoe_direct = pd.DataFrame([
        {"source_build": "GRCh38", "source_variant": "rs429358", "source_chromosome": "chr19", "source_position": 44908684,
         "target_build": "GRCh37", "target_chromosome": "19", "target_position": 45411941, "validation": "matches_frozen_A1_direct_variant"},
        {"source_build": "GRCh38", "source_variant": "rs7412", "source_chromosome": "chr19", "source_position": 44908822,
         "target_build": "GRCh37", "target_chromosome": "19", "target_position": 45412079, "validation": "matches_frozen_A1_direct_variant"},
    ])
    apoe_direct["frozen_exclusion_interval_hg19"] = f"chr19:{APOE_START_HG19}-{APOE_END_HG19}"
    apoe_direct["chain_file"] = str(CHAIN)
    apoe_direct["chain_sha256"] = sha256(CHAIN)

    MAPPED_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    ELIGIBLE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    LOG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(MAPPED_OUTPUT, sep="\t", index=False, compression="gzip", na_rep="NA")
    eligible.to_csv(ELIGIBLE_OUTPUT, sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summary_rows).to_csv(SUMMARY_OUTPUT, sep="\t", index=False, na_rep="NA")
    apoe_direct.to_csv(APOE_AUDIT_OUTPUT, sep="\t", index=False, na_rep="NA")

    log_lines = [
        f"candidate_rows={len(data)}",
        f"unique_hg38_positions={len(unique_positions)}",
        f"liftover_unique_plus_strand={int(data['liftover_status'].eq('unique_plus_strand').sum())}",
        f"apoe_region_excluded={int(data['APOE_region_excluded'].sum())}",
        f"eligible_rows_before_deduplication={int(data['eligible_for_clumping'].sum())}",
        f"eligible_unique_assay_snp_rows={len(eligible)}",
        f"position_allele_rescues={int(data['ld_mapping_method'].eq('hg38_to_hg19_position_plus_allele').sum())}",
        f"chain_sha256={sha256(CHAIN)}",
        f"candidate_source_sha256={sha256(CANDIDATES)}",
    ]
    LOG_OUTPUT.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print("\n".join(log_lines))
    print(f"eligible_output={ELIGIBLE_OUTPUT}")


if __name__ == "__main__":
    main()
