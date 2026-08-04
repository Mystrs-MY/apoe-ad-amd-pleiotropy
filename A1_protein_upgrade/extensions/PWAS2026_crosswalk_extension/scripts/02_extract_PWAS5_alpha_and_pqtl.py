#!/usr/bin/env python3
"""Extract APOE alpha rows and genome-wide significant pQTLs for the frozen PWAS5 set."""

from __future__ import annotations

import csv
import gzip
import io
import math
import re
import tarfile
from pathlib import Path

import pandas as pd
import yaml


EXT_ROOT = Path(__file__).resolve().parents[1]
UPGRADE_ROOT = EXT_ROOT.parents[1]
RESOURCE_CONFIG = yaml.safe_load((UPGRADE_ROOT / "config" / "resources.yml").read_text(encoding="utf-8"))
RESOURCE_DIR = next(
    Path(path) for path in RESOURCE_CONFIG["paths"]["ukb_ppp_tar_dirs"]
    if Path(path).name == "syn51365303_European_discovery"
)
RSID_DIR = Path(RESOURCE_CONFIG["paths"]["rsid_map_dir"])
MEMBERS = EXT_ROOT / "config" / "PWAS5_frozen_members.tsv"
COORDINATES = EXT_ROOT / "config" / "PWAS5_gene_coordinates_grch37.tsv"
INTEGRITY_QA = EXT_ROOT / "tables" / "PWAS5_integrity_QA.tsv"
ALPHA_OUT = EXT_ROOT / "tables" / "PWAS5_APOE_alpha.tsv"
CANDIDATE_OUT = EXT_ROOT / "data_processed" / "PWAS5_pqtl_candidates.tsv"
EXTRACTION_QA_OUT = EXT_ROOT / "tables" / "PWAS5_extraction_QA.tsv"
CROSSWALK_OUT = EXT_ROOT / "tables" / "PWAS5_crosswalk_mapping.tsv"

P_THRESHOLD = 5e-8
LOG10P_THRESHOLD = -math.log10(P_THRESHOLD)
MAF_THRESHOLD = 0.01
INFO_THRESHOLD = 0.8
F_THRESHOLD = 10.0
APOE_START = 44_000_000
APOE_END = 46_500_000
APOE_VARIANTS = {
    "rs429358": {"position": 45_411_941, "effect_allele": "C", "other_allele": "T"},
    "rs7412": {"position": 45_412_079, "effect_allele": "T", "other_allele": "C"},
}


def chromosome_from_member(name: str) -> str:
    match = re.search(r"discovery_chr([0-9XY]+)_", name, re.IGNORECASE)
    return match.group(1).upper() if match else "unresolved"


def parse_variant_id(value: str) -> tuple[int, int]:
    parts = value.split(":")
    if len(parts) < 2:
        raise ValueError(f"Unexpected UKB-PPP variant ID: {value}")
    return int(parts[0]), int(parts[1])


def p_value(log10p: str) -> float:
    value = float(log10p)
    return 0.0 if value > 323 else 10 ** (-value)


def add_rsids(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["SNP"] = "mapping_unresolved"
    for chromosome, group in frame.groupby("chromosome", sort=True):
        wanted = set(group["source_variant_ID"].astype(str))
        map_path = RSID_DIR / f"olink_rsid_map_mac5_info03_b0_7_chr{chromosome}_patched_v2.tsv.gz"
        if not map_path.exists():
            raise FileNotFoundError(map_path)
        matches: dict[str, str] = {}
        for chunk in pd.read_csv(map_path, sep="\t", usecols=["ID", "rsid"], dtype=str, chunksize=500_000):
            selected = chunk[chunk["ID"].isin(wanted)]
            if len(selected):
                matches.update(dict(zip(selected["ID"], selected["rsid"])))
            if len(matches) == len(wanted):
                break
        index = group.index
        frame.loc[index, "SNP"] = frame.loc[index, "source_variant_ID"].map(matches).fillna("mapping_unresolved")
    frame["rsid_mapping_status"] = frame["SNP"].ne("mapping_unresolved").map({True: "mapped", False: "mapping_unresolved"})
    return frame


def orient_alpha(row: dict[str, str], variant: str) -> dict[str, object]:
    specification = APOE_VARIANTS[variant]
    original_effect = row["ALLELE1"].upper()
    original_other = row["ALLELE0"].upper()
    requested_effect = specification["effect_allele"]
    requested_other = specification["other_allele"]
    beta = float(row["BETA"])
    eaf = float(row["A1FREQ"])
    if (original_effect, original_other) == (requested_effect, requested_other):
        flipped = False
    elif (original_effect, original_other) == (requested_other, requested_effect):
        beta = -beta
        eaf = 1 - eaf
        flipped = True
    else:
        raise ValueError(
            f"APOE allele mismatch for {variant}: source={original_effect}/{original_other}, "
            f"requested={requested_effect}/{requested_other}"
        )
    chromosome, position = parse_variant_id(row["ID"])
    return {
        "variant": variant,
        "chromosome": chromosome,
        "position_hg19_from_ID": position,
        "position_hg38_from_GENPOS": int(row["GENPOS"]),
        "requested_effect_allele": requested_effect,
        "requested_other_allele": requested_other,
        "original_effect_allele": original_effect,
        "original_other_allele": original_other,
        "allele_flipped": flipped,
        "beta": beta,
        "SE": float(row["SE"]),
        "P_value": p_value(row["LOG10P"]),
        "LOG10P": float(row["LOG10P"]),
        "effect_allele_frequency": eaf,
        "effect_allele_frequency_original": float(row["A1FREQ"]),
        "sample_size": int(row["N"]),
        "INFO": float(row["INFO"]),
        "source_variant_ID": row["ID"],
        "availability_status": "direct_variant_available",
        "proxy_used": False,
    }


def scan_archive(spec: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    path = RESOURCE_DIR / spec["archive_name"]
    alpha_hits: dict[str, list[dict[str, str]]] = {variant: [] for variant in APOE_VARIANTS}
    candidates: list[dict[str, object]] = []
    scanned_rows = 0
    with tarfile.open(path, "r") as archive:
        members = [m for m in archive.getmembers() if m.isfile() and m.name.endswith(".gz")]
        for member in members:
            chromosome_label = chromosome_from_member(member.name)
            if chromosome_label in {"X", "Y", "XY", "UNRESOLVED"}:
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise OSError(f"Cannot read {member.name}")
            with gzip.GzipFile(fileobj=stream) as compressed:
                text = io.TextIOWrapper(compressed, encoding="utf-8")
                header = text.readline().strip().split()
                for line in text:
                    if not line.strip():
                        continue
                    scanned_rows += 1
                    values = dict(zip(header, line.strip().split()))
                    chromosome, position = parse_variant_id(values["ID"])
                    if chromosome == 19:
                        for variant, requested in APOE_VARIANTS.items():
                            if position == requested["position"]:
                                alpha_hits[variant].append(values)
                    if values.get("TEST") != "ADD":
                        continue
                    log10p = float(values["LOG10P"])
                    if log10p < LOG10P_THRESHOLD:
                        continue
                    eaf = float(values["A1FREQ"])
                    maf = min(eaf, 1 - eaf)
                    info = float(values["INFO"])
                    beta = float(values["BETA"])
                    se = float(values["SE"])
                    f_statistic = (beta / se) ** 2 if math.isfinite(beta) and math.isfinite(se) and se > 0 else math.nan
                    if maf < MAF_THRESHOLD or info < INFO_THRESHOLD or not math.isfinite(f_statistic) or f_statistic <= F_THRESHOLD:
                        continue
                    if chromosome == 19 and APOE_START <= position <= APOE_END:
                        continue
                    candidates.append({
                        "gene_symbol": spec["gene_symbol"],
                        "UniProt_ID": spec["UniProt_ID"],
                        "assay_target_ID": spec["UKB_PPP_OID"],
                        "Olink_panel": spec["Olink_panel"],
                        "source_tar": str(path),
                        "source_member": member.name,
                        "source_variant_ID": values["ID"],
                        "chromosome": chromosome,
                        "position_hg19_from_ID": position,
                        "position_hg38_from_GENPOS": int(values["GENPOS"]),
                        "effect_allele": values["ALLELE1"].upper(),
                        "other_allele": values["ALLELE0"].upper(),
                        "effect_allele_frequency": eaf,
                        "MAF": maf,
                        "INFO": info,
                        "TEST": values["TEST"],
                        "sample_size": int(values["N"]),
                        "beta": beta,
                        "SE": se,
                        "F_statistic": f_statistic,
                        "P_value": p_value(values["LOG10P"]),
                        "LOG10P": log10p,
                        "APOE_region_excluded": True,
                        "QC_rule": "ADD; P<5e-8; MAF>=0.01; INFO>=0.8; F>10; APOE region excluded",
                    })

    alpha_rows: list[dict[str, object]] = []
    for variant in APOE_VARIANTS:
        hits = alpha_hits[variant]
        if len(hits) == 1:
            row = orient_alpha(hits[0], variant)
        else:
            row = {
                "variant": variant,
                "chromosome": 19,
                "position_hg19_from_ID": APOE_VARIANTS[variant]["position"],
                "requested_effect_allele": APOE_VARIANTS[variant]["effect_allele"],
                "requested_other_allele": APOE_VARIANTS[variant]["other_allele"],
                "availability_status": "missing" if len(hits) == 0 else "duplicate_unresolved",
                "proxy_used": False,
            }
        row.update({
            "gene_symbol": spec["gene_symbol"],
            "protein_name": spec["literature_protein_name"],
            "UniProt_ID": spec["UniProt_ID"],
            "assay_target_name": spec["UKB_PPP_target_name"],
            "assay_target_ID": spec["UKB_PPP_OID"],
            "Olink_panel": spec["Olink_panel"],
            "source_tar": str(path),
            "same_assay_extension": True,
            "variant_match_count": len(hits),
        })
        alpha_rows.append(row)
    qa = {
        "gene_symbol": spec["gene_symbol"],
        "UKB_PPP_OID": spec["UKB_PPP_OID"],
        "rows_scanned_autosomes": scanned_rows,
        "rs429358_match_count": len(alpha_hits["rs429358"]),
        "rs7412_match_count": len(alpha_hits["rs7412"]),
        "n_pqtl_candidates_after_APOE_exclusion": len(candidates),
        "candidate_QC_rule": "ADD; P<5e-8; MAF>=0.01; INFO>=0.8; F>10; APOE region excluded",
        "alpha_unique_pass": len(alpha_hits["rs429358"]) == 1 and len(alpha_hits["rs7412"]) == 1,
    }
    return alpha_rows, candidates, qa


def main() -> None:
    integrity = pd.read_csv(INTEGRITY_QA, sep="\t")
    if len(integrity) != 5 or not integrity["overall_pass"].astype(str).str.lower().eq("true").all():
        raise RuntimeError("PWAS5 archive integrity must pass before extraction")
    members = pd.read_csv(MEMBERS, sep="\t", dtype=str).sort_values("analysis_order")
    coordinates = pd.read_csv(COORDINATES, sep="\t").set_index("gene_symbol")
    alpha_rows: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    qa_rows: list[dict[str, object]] = []
    for spec in members.to_dict("records"):
        print(f"Extracting {spec['gene_symbol']}", flush=True)
        alpha, candidate, qa = scan_archive(spec)
        alpha_rows.extend(alpha)
        candidates.extend(candidate)
        qa_rows.append(qa)

    alpha = pd.DataFrame(alpha_rows)
    candidate_frame = pd.DataFrame(candidates)
    if not candidate_frame.empty:
        candidate_frame = add_rsids(candidate_frame)
        candidate_frame["cis_to_encoding_gene"] = False
        candidate_frame["cis_window_definition"] = "encoding gene GRCh37 boundaries +/-1 Mb"
        candidate_frame["gene_coordinate_source"] = "NA"
        for gene in members["gene_symbol"]:
            coord = coordinates.loc[gene]
            lower = int(coord["gene_start_grch37"]) - int(coord["cis_window_bp"])
            upper = int(coord["gene_end_grch37"]) + int(coord["cis_window_bp"])
            gene_mask = candidate_frame["gene_symbol"].eq(gene)
            cis_mask = (
                gene_mask
                & candidate_frame["chromosome"].eq(int(coord["chromosome"]))
                & candidate_frame["position_hg19_from_ID"].between(lower, upper)
            )
            candidate_frame.loc[cis_mask, "cis_to_encoding_gene"] = True
            candidate_frame.loc[gene_mask, "gene_coordinate_source"] = coord["source_url"]
        candidate_frame = candidate_frame.sort_values(["gene_symbol", "chromosome", "position_hg19_from_ID"])

    crosswalk = members.rename(columns={
        "UKB_PPP_OID": "assay_target_ID",
        "UKB_PPP_target_name": "assay_target_name",
    }).copy()
    crosswalk["candidate_assay_count"] = 1
    crosswalk["protein_form"] = "circulating_plasma_protein; platform-specific binding assay"
    crosswalk["cross_platform_correlation"] = "not_reported_in_prespecified_crosswalk"
    crosswalk["exact_assay_replication"] = False
    crosswalk["interpretation_boundary"] = (
        "Gene/protein-level cross-platform crosswalk; not exact SomaScan aptamer replication by Olink."
    )

    ALPHA_OUT.parent.mkdir(parents=True, exist_ok=True)
    CANDIDATE_OUT.parent.mkdir(parents=True, exist_ok=True)
    alpha.to_csv(ALPHA_OUT, sep="\t", index=False, na_rep="NA")
    candidate_frame.to_csv(CANDIDATE_OUT, sep="\t", index=False, na_rep="NA")
    pd.DataFrame(qa_rows).to_csv(EXTRACTION_QA_OUT, sep="\t", index=False)
    crosswalk.to_csv(CROSSWALK_OUT, sep="\t", index=False)

    if len(alpha) != 10:
        raise RuntimeError(f"Expected 10 alpha rows; observed {len(alpha)}")
    if not alpha["variant_match_count"].eq(1).all():
        raise RuntimeError("At least one APOE alpha variant was missing or duplicated")
    print(f"Alpha rows: {len(alpha)}; pQTL candidates: {len(candidate_frame)}")
    print(f"Mapped rsIDs: {candidate_frame['rsid_mapping_status'].eq('mapped').sum()}")


if __name__ == "__main__":
    main()
