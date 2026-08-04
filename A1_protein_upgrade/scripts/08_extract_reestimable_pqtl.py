#!/usr/bin/env python3
"""Extract genome-wide significant pQTL candidates for re-estimable primary proteins."""

from __future__ import annotations

import gzip
import io
import math
import re
import tarfile
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "config" / "resources.yml").read_text(encoding="utf-8"))
TAR_DIR = Path(CONFIG["paths"]["ukb_ppp_tar_dir"])
RSID_DIR = Path(CONFIG["paths"]["rsid_map_dir"])
ALPHA = ROOT / "tables" / "APOE_variant_to_literature_proteins_alpha.tsv"
MANIFEST = ROOT / "data_processed" / "ukbppp_local_assay_manifest.tsv"
OUTPUT = ROOT / "data_processed" / "literature_panel_pqtl_candidates.tsv"
LOG = ROOT / "logs" / "pqtl_candidate_extraction.tsv"
COORDINATES = ROOT / "config" / "protein_gene_coordinates_grch37.tsv"

P_THRESHOLD = 5e-8
LOG10P_THRESHOLD = -math.log10(P_THRESHOLD)
APOE_START = 44_000_000
APOE_END = 46_500_000


def chromosome_from_member(name: str) -> str:
    match = re.search(r"discovery_chr([0-9XY]+)_", name, re.IGNORECASE)
    return match.group(1).upper() if match else "unresolved"


def scan_tar(gene: str, tar_path: Path, assay: pd.Series) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with tarfile.open(tar_path, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile() and member.name.endswith(".gz")]
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
                    values = dict(zip(header, line.strip().split()))
                    if float(values["LOG10P"]) < LOG10P_THRESHOLD:
                        continue
                    parts = values["ID"].split(":")
                    if len(parts) < 2:
                        continue
                    chromosome = int(parts[0])
                    position_hg19 = int(parts[1])
                    if chromosome == 19 and APOE_START <= position_hg19 <= APOE_END:
                        continue
                    rows.append({
                        "gene_symbol": gene,
                        "UniProt_ID": assay["UniProt_ID"],
                        "assay_target_ID": assay["Olink_target_ID"],
                        "Olink_panel": assay["Olink_panel"],
                        "source_tar": str(tar_path),
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
                        "P_value": 0.0 if float(values["LOG10P"]) > 323 else 10 ** (-float(values["LOG10P"])),
                        "LOG10P": float(values["LOG10P"]),
                        "APOE_region_excluded": "true",
                    })
    return rows


def add_rsids(frame: pd.DataFrame) -> pd.DataFrame:
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
        indexes = group.index
        frame.loc[indexes, "SNP"] = frame.loc[indexes, "source_variant_ID"].map(matches).fillna("mapping_unresolved")
    return frame


def main() -> None:
    alpha = pd.read_csv(ALPHA, sep="\t", dtype=str).fillna("NA")
    mapping_eligible = alpha["eligible_for_two_step_mapping"].astype(str).str.lower().eq("true")
    genes = sorted(alpha.loc[
        (alpha["availability_status"] == "direct_variant_available") & mapping_eligible,
        "gene_symbol",
    ].unique())
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).set_index("gene_symbol")
    rows: list[dict[str, object]] = []
    logs = []
    for gene in genes:
        assay = manifest.loc[gene]
        tar_path = Path(assay["tar_path"])
        extracted = scan_tar(gene, tar_path, assay)
        rows.extend(extracted)
        logs.append({"gene_symbol": gene, "source_tar": str(tar_path), "n_p_lt_5e8_after_APOE_exclusion": len(extracted), "status": "complete"})
        print(f"{gene}: {len(extracted)} candidates")
    frame = pd.DataFrame(rows)
    if len(frame):
        frame = add_rsids(frame)
        coordinates = pd.read_csv(COORDINATES, sep="\t").set_index("gene_symbol")
        frame["cis_to_encoding_gene"] = False
        frame["cis_window_definition"] = "encoding gene GRCh37 boundaries +/-1 Mb"
        frame["gene_coordinate_source"] = "NA"
        for gene in genes:
            coordinate = coordinates.loc[gene]
            lower = int(coordinate["gene_start_grch37"]) - int(coordinate["cis_window_bp"])
            upper = int(coordinate["gene_end_grch37"]) + int(coordinate["cis_window_bp"])
            mask = (
                (frame["gene_symbol"] == gene)
                & (frame["chromosome"] == int(coordinate["chromosome"]))
                & frame["position_hg19_from_ID"].between(lower, upper)
            )
            frame.loc[mask, "cis_to_encoding_gene"] = True
            frame.loc[frame["gene_symbol"] == gene, "gene_coordinate_source"] = coordinate["source_url"]
        frame["rsid_mapping_status"] = frame["SNP"].where(frame["SNP"] != "mapping_unresolved", "mapping_unresolved")
        frame.loc[frame["SNP"] != "mapping_unresolved", "rsid_mapping_status"] = "mapped"
        frame = frame.sort_values(["gene_symbol", "chromosome", "position_hg19_from_ID"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUTPUT, sep="\t", index=False)
    pd.DataFrame(logs).to_csv(LOG, sep="\t", index=False)
    print(f"Total candidates: {len(frame)}")
    print(f"Mapped rsIDs: {(frame['rsid_mapping_status'] == 'mapped').sum() if len(frame) else 0}")
    if len(frame):
        print("Cis candidates: " + ", ".join(f"{gene}={count}" for gene, count in frame[frame['cis_to_encoding_gene']].groupby('gene_symbol').size().items()))
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
