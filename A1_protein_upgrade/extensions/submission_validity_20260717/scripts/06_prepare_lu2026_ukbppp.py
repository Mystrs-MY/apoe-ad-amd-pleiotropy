#!/usr/bin/env python3
"""Prepare exact-assay UKB-PPP inputs for the bounded Lu 2026 extension."""

from __future__ import annotations

import csv
import gzip
import io
import math
import re
import tarfile
import time
from datetime import date
from pathlib import Path

import pandas as pd
import requests
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
EXTENSION_ROOT = SCRIPT_DIR.parent
UPGRADE_ROOT = EXTENSION_ROOT.parents[1]

SELECTION = EXTENSION_ROOT / "config" / "lu2026_exact_assay_download_selection.tsv"
TAR_DIR = EXTENSION_ROOT / "data_raw" / "lu2026_ukbppp_exact_assays"
PROCESSED_DIR = EXTENSION_ROOT / "data_processed"
TABLE_DIR = EXTENSION_ROOT / "tables"
LOG_DIR = EXTENSION_ROOT / "logs"
RESOURCES = yaml.safe_load((UPGRADE_ROOT / "config" / "resources.yml").read_text(encoding="utf-8"))
RSID_DIR = Path(RESOURCES["paths"]["rsid_map_dir"])

P_THRESHOLD = 5e-8
LOG10P_THRESHOLD = -math.log10(P_THRESHOLD)
APOE_START = 44_000_000
APOE_END = 46_500_000
ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens"

VARIANTS = {
    "rs429358": {"chromosome": 19, "position_hg19": 45_411_941, "effect_allele": "C", "other_allele": "T"},
    "rs7412": {"chromosome": 19, "position_hg19": 45_412_079, "effect_allele": "T", "other_allele": "C"},
}


def extension_relative(path: Path) -> str:
    return path.resolve().relative_to(EXTENSION_ROOT.resolve()).as_posix()


def p_from_log10(value: str) -> float:
    log10p = float(value)
    return 0.0 if log10p > 323 else 10 ** (-log10p)


def parse_tar_name(path: Path) -> dict[str, str]:
    match = re.match(
        r"^(?P<gene>[^_]+)_(?P<uniprot>[^_]+)_(?P<oid>OID\d+)_v(?P<version>\d+)_(?P<panel>.+)\.tar$",
        path.name,
        re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"UKB-PPP filename cannot be parsed at assay level: {path.name}")
    return {key: value.upper() if key in {"gene", "uniprot", "oid"} else value for key, value in match.groupdict().items()}


def locate_selected_tars(selection: pd.DataFrame) -> dict[str, tuple[Path, dict[str, str]]]:
    found: dict[str, tuple[Path, dict[str, str]]] = {}
    for path in sorted(TAR_DIR.glob("*.tar")):
        meta = parse_tar_name(path)
        selected = selection[
            selection["gene_symbol"].str.upper().eq(meta["gene"])
            & selection["UniProt_ID"].str.upper().eq(meta["uniprot"])
            & selection["Olink_target_ID"].str.upper().eq(meta["oid"])
        ]
        if len(selected) == 1:
            if meta["gene"] in found:
                raise ValueError(f"Multiple exact tar files found for {meta['gene']}")
            found[meta["gene"]] = (path, meta)
    missing = sorted(set(selection["gene_symbol"].str.upper()) - set(found))
    if missing:
        raise FileNotFoundError("Verified Synapse tar files are missing: " + ", ".join(missing))
    return found


def fetch_coordinates(genes: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for gene in genes:
        url = f"{ENSEMBL_GRCH37}/{gene}?content-type=application/json"
        last_error = ""
        for attempt in range(1, 6):
            try:
                response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
                if response.status_code != 200:
                    raise RuntimeError(f"HTTP {response.status_code}")
                record = response.json()
                rows.append({
                    "gene_symbol": gene,
                    "Ensembl_gene_ID": record.get("id", "mapping_unresolved"),
                    "chromosome": record.get("seq_region_name", "mapping_unresolved"),
                    "gene_start_grch37": record.get("start", "mapping_unresolved"),
                    "gene_end_grch37": record.get("end", "mapping_unresolved"),
                    "strand": record.get("strand", "mapping_unresolved"),
                    "cis_window_bp": 1_000_000,
                    "source_url": url,
                    "verification_date": date.today().isoformat(),
                    "mapping_status": "verified_ensembl_grch37",
                })
                break
            except Exception as exc:  # network failures are retained, not guessed
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt < 5:
                    time.sleep(3 * attempt)
        else:
            rows.append({
                "gene_symbol": gene,
                "Ensembl_gene_ID": "mapping_unresolved",
                "chromosome": "mapping_unresolved",
                "gene_start_grch37": "mapping_unresolved",
                "gene_end_grch37": "mapping_unresolved",
                "strand": "mapping_unresolved",
                "cis_window_bp": 1_000_000,
                "source_url": url,
                "verification_date": date.today().isoformat(),
                "mapping_status": "requires_manual_verification",
                "error": last_error,
            })
    frame = pd.DataFrame(rows)
    if not frame["mapping_status"].eq("verified_ensembl_grch37").all():
        unresolved = ", ".join(frame.loc[frame["mapping_status"] != "verified_ensembl_grch37", "gene_symbol"])
        raise RuntimeError(f"GRCh37 coordinates unresolved; cis analysis stopped for: {unresolved}")
    return frame


def chromosome_from_member(name: str) -> str:
    match = re.search(r"discovery_chr([0-9XY]+)_", name, re.IGNORECASE)
    return match.group(1).upper() if match else "UNRESOLVED"


def scan_tar(gene: str, tar_path: Path, meta: dict[str, str]) -> tuple[list[dict[str, object]], dict[str, dict[str, str]]]:
    candidates: list[dict[str, object]] = []
    alpha_rows: dict[str, dict[str, str]] = {}
    with tarfile.open(tar_path, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile() and member.name.endswith(".gz")]
        if not members:
            raise ValueError(f"No chromosome gzip members found in {tar_path.name}")
        for member in members:
            chromosome_label = chromosome_from_member(member.name)
            if chromosome_label in {"X", "Y", "XY", "UNRESOLVED"}:
                continue
            member_stream = archive.extractfile(member)
            if member_stream is None:
                raise OSError(f"Cannot read {member.name}")
            with gzip.GzipFile(fileobj=member_stream) as compressed:
                text = io.TextIOWrapper(compressed, encoding="utf-8")
                header = text.readline().strip().split()
                required = {"CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "LOG10P"}
                if not required.issubset(header):
                    raise ValueError(f"Unexpected UKB-PPP columns in {member.name}")
                for line in text:
                    values = dict(zip(header, line.strip().split()))
                    parts = values["ID"].split(":")
                    if len(parts) < 2:
                        continue
                    chromosome = int(parts[0])
                    position_hg19 = int(parts[1])
                    if chromosome == 19:
                        for rsid, specification in VARIANTS.items():
                            if position_hg19 == specification["position_hg19"]:
                                alpha_rows[rsid] = values.copy()
                    if float(values["LOG10P"]) < LOG10P_THRESHOLD:
                        continue
                    if chromosome == 19 and APOE_START <= position_hg19 <= APOE_END:
                        continue
                    candidates.append({
                        "gene_symbol": gene,
                        "UniProt_ID": meta["uniprot"],
                        "assay_target_ID": meta["oid"],
                        "Olink_panel": meta["panel"],
                        "source_tar": extension_relative(tar_path),
                        "source_member": member.name,
                        "source_variant_ID": values["ID"],
                        "chromosome": chromosome,
                        "position_hg19_from_ID": position_hg19,
                        "position_hg38_from_GENPOS": int(values["GENPOS"]),
                        "effect_allele": values["ALLELE1"].upper(),
                        "other_allele": values["ALLELE0"].upper(),
                        "effect_allele_frequency": float(values["A1FREQ"]),
                        "INFO": float(values["INFO"]),
                        "sample_size": int(values["N"]),
                        "beta": float(values["BETA"]),
                        "SE": float(values["SE"]),
                        "P_value": p_from_log10(values["LOG10P"]),
                        "LOG10P": float(values["LOG10P"]),
                        "APOE_region_excluded": "true",
                    })
    return candidates, alpha_rows


def add_rsids(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["SNP"] = "mapping_unresolved"
    for chromosome, group in frame.groupby("chromosome"):
        wanted = set(group["source_variant_ID"])
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
        frame.loc[group.index, "SNP"] = frame.loc[group.index, "source_variant_ID"].map(matches).fillna("mapping_unresolved")
    frame["rsid_mapping_status"] = frame["SNP"].ne("mapping_unresolved").map({True: "mapped", False: "mapping_unresolved"})
    return frame


def align_alpha(gene: str, meta: dict[str, str], tar_path: Path, found: dict[str, dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rsid, specification in VARIANTS.items():
        source = found.get(rsid)
        if source is None:
            rows.append({
                "gene_symbol": gene,
                "variant": rsid,
                "requested_effect_allele": specification["effect_allele"],
                "requested_other_allele": specification["other_allele"],
                "UniProt_ID": meta["uniprot"],
                "assay_target_ID": meta["oid"],
                "availability_status": "direct_variant_absent_in_assay_GWAS",
                "exclusion_reason": f"{rsid} was not found at the verified GRCh37 position.",
                "alpha_source": extension_relative(tar_path),
            })
            continue
        original_effect = source["ALLELE1"].upper()
        original_other = source["ALLELE0"].upper()
        requested_effect = specification["effect_allele"]
        requested_other = specification["other_allele"]
        if original_effect == requested_effect and original_other == requested_other:
            beta = float(source["BETA"])
            allele_flipped = False
        elif original_effect == requested_other and original_other == requested_effect:
            beta = -float(source["BETA"])
            allele_flipped = True
        else:
            rows.append({
                "gene_symbol": gene,
                "variant": rsid,
                "requested_effect_allele": requested_effect,
                "requested_other_allele": requested_other,
                "original_effect_allele": original_effect,
                "original_other_allele": original_other,
                "UniProt_ID": meta["uniprot"],
                "assay_target_ID": meta["oid"],
                "availability_status": "alleles_not_harmonizable",
                "exclusion_reason": "Observed alleles do not match the prespecified APOE alleles.",
                "alpha_source": extension_relative(tar_path),
            })
            continue
        rows.append({
            "gene_symbol": gene,
            "protein_name": gene,
            "protein_form_or_isoform": "Olink assay-level target",
            "variant": rsid,
            "chromosome": 19,
            "position_hg19": specification["position_hg19"],
            "position_hg38_from_GENPOS": source["GENPOS"],
            "requested_effect_allele": requested_effect,
            "requested_other_allele": requested_other,
            "original_effect_allele": original_effect,
            "original_other_allele": original_other,
            "allele_flipped": allele_flipped,
            "beta": beta,
            "SE": float(source["SE"]),
            "P_value": p_from_log10(source["LOG10P"]),
            "effect_allele_frequency_original": float(source["A1FREQ"]),
            "sample_size": int(source["N"]),
            "imputation_INFO": float(source["INFO"]),
            "source_variant_ID": source["ID"],
            "proxy_used": False,
            "proteomic_platform": "UKB-PPP Olink Explore 3072",
            "UniProt_ID": meta["uniprot"],
            "assay_target_ID": meta["oid"],
            "assay_version": f"v{meta['version']}",
            "Olink_panel": meta["panel"],
            "effect_unit": "per genetically predicted Olink NPX unit as reported by UKB-PPP",
            "alpha_source": extension_relative(tar_path),
            "availability_status": "direct_variant_available",
            "exclusion_reason": "NA",
            "analysis_label": "Lu_2026_exact_assay_external_candidate_sensitivity",
        })
    return rows


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    selection = pd.read_csv(SELECTION, sep="\t", dtype=str).fillna("NA")
    selection = selection[selection["selection_status"] == "selected_for_download"].copy()
    if len(selection) != 5:
        raise ValueError(f"Expected exactly five gate-passing Lu candidates; observed {len(selection)}")

    selected_tars = locate_selected_tars(selection)
    coordinates = fetch_coordinates(sorted(selected_tars))
    coordinates.to_csv(EXTENSION_ROOT / "config" / "lu2026_gene_coordinates_grch37.tsv", sep="\t", index=False)
    coordinate_map = coordinates.set_index("gene_symbol")

    candidate_rows: list[dict[str, object]] = []
    alpha_output: list[dict[str, object]] = []
    log_rows: list[dict[str, object]] = []
    for gene in sorted(selected_tars):
        tar_path, meta = selected_tars[gene]
        print(f"Scanning {gene}: {tar_path.name}", flush=True)
        candidates, alpha_found = scan_tar(gene, tar_path, meta)
        candidate_rows.extend(candidates)
        alpha_output.extend(align_alpha(gene, meta, tar_path, alpha_found))
        log_rows.append({
            "gene_symbol": gene,
            "source_tar": extension_relative(tar_path),
            "n_p_lt_5e8_after_APOE_exclusion": len(candidates),
            "direct_APOE_variants_found": len(alpha_found),
            "status": "complete",
        })
        print(f"{gene}: {len(candidates)} pQTL candidates; {len(alpha_found)}/2 APOE variants", flush=True)

    candidates = pd.DataFrame(candidate_rows)
    if candidates.empty:
        raise RuntimeError("No genome-wide significant pQTL candidates were extracted for Lu 2026 exact assays.")
    candidates = add_rsids(candidates)
    candidates["cis_to_encoding_gene"] = False
    candidates["cis_window_definition"] = "encoding gene GRCh37 boundaries +/-1 Mb"
    candidates["gene_coordinate_source"] = "NA"
    for gene in sorted(selected_tars):
        coordinate = coordinate_map.loc[gene]
        lower = int(coordinate["gene_start_grch37"]) - int(coordinate["cis_window_bp"])
        upper = int(coordinate["gene_end_grch37"]) + int(coordinate["cis_window_bp"])
        mask = (
            candidates["gene_symbol"].eq(gene)
            & candidates["chromosome"].eq(int(coordinate["chromosome"]))
            & candidates["position_hg19_from_ID"].between(lower, upper)
        )
        candidates.loc[mask, "cis_to_encoding_gene"] = True
        candidates.loc[candidates["gene_symbol"].eq(gene), "gene_coordinate_source"] = coordinate["source_url"]
    candidates = candidates.sort_values(["gene_symbol", "chromosome", "position_hg19_from_ID"])
    candidates.to_csv(PROCESSED_DIR / "Lu2026_exact_assay_pqtl_candidates.tsv", sep="\t", index=False)

    alpha = pd.DataFrame(alpha_output)
    alpha.to_csv(TABLE_DIR / "Lu2026_APOE_variant_to_protein_alpha.tsv", sep="\t", index=False)
    pd.DataFrame(log_rows).to_csv(LOG_DIR / "Lu2026_UKB_PPP_extraction_log.tsv", sep="\t", index=False)

    mapped = int(candidates["rsid_mapping_status"].eq("mapped").sum())
    direct_alpha = int(alpha["availability_status"].eq("direct_variant_available").sum())
    print(f"Total pQTL candidates: {len(candidates)}; mapped rsIDs: {mapped}")
    print(f"Direct alpha rows: {direct_alpha}/10")


if __name__ == "__main__":
    main()
