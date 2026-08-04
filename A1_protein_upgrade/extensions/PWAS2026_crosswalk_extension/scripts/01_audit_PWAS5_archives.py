#!/usr/bin/env python3
"""Audit the five frozen 2026 PWAS crosswalk UKB-PPP archives without extraction."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import tarfile
import os
from datetime import datetime
from pathlib import Path


EXT_ROOT = Path(__file__).resolve().parents[1]
UPGRADE_ROOT = EXT_ROOT.parents[1]
PROJECT_ROOT = EXT_ROOT.parents[2]
RESOURCE_ROOT = Path(os.environ.get("A1_RESOURCE_ROOT", PROJECT_ROOT / "data" / "external"))
RESOURCE_DIR = Path(os.environ.get(
    "UKB_PPP_DOWNLOAD_DIR",
    RESOURCE_ROOT / "UKB-PPP" / "syn51365303_European_discovery",
))
CONFIG = EXT_ROOT / "config" / "PWAS5_frozen_members.tsv"
TABLE_DIR = EXT_ROOT / "tables"
LOG_DIR = EXT_ROOT / "logs"
MANIFEST_OUT = TABLE_DIR / "PWAS5_file_manifest.tsv"
QA_OUT = TABLE_DIR / "PWAS5_integrity_QA.tsv"
MEMBER_OUT = TABLE_DIR / "PWAS5_member_integrity_detail.tsv"
EXPECTED_HEADER = [
    "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO",
    "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA",
]


def digest(path: Path, algorithm: str) -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest().upper()


def chromosome_from_member(name: str) -> str:
    match = re.search(r"discovery_chr([0-9XY]+)_", name, re.IGNORECASE)
    return match.group(1).upper() if match else "unresolved"


def inspect_member(archive: tarfile.TarFile, member: tarfile.TarInfo) -> dict[str, object]:
    stream = archive.extractfile(member)
    if stream is None:
        raise OSError(f"Cannot open {member.name}")
    row_count = 0
    first_data = ""
    last_data = ""
    with gzip.GzipFile(fileobj=stream) as compressed:
        text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
        header_line = text.readline().strip()
        header = header_line.split()
        for line in text:
            stripped = line.strip()
            if not stripped:
                continue
            row_count += 1
            if not first_data:
                first_data = stripped
            last_data = stripped
    first_fields = first_data.split()
    last_fields = last_data.split()
    return {
        "member_name": member.name,
        "chromosome_label": chromosome_from_member(member.name),
        "member_size_bytes": member.size,
        "header": header_line,
        "header_exact_match": header == EXPECTED_HEADER,
        "row_count": row_count,
        "first_variant_ID": first_fields[0] if first_fields else "NA",
        "last_variant_ID": last_fields[0] if last_fields else "NA",
        "first_row_field_count": len(first_fields),
        "last_row_field_count": len(last_fields),
        "first_and_last_complete": len(first_fields) == len(header) and len(last_fields) == len(header),
    }


def read_config() -> list[dict[str, str]]:
    with CONFIG.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"Refusing to write empty audit table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    qa_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []

    for spec in read_config():
        path = RESOURCE_DIR / spec["archive_name"]
        if not path.exists():
            raise FileNotFoundError(path)
        print(f"Auditing {path.name}", flush=True)
        stat = path.stat()
        md5 = digest(path, "md5")
        sha256 = digest(path, "sha256")
        tar_readable = False
        gzip_members_readable = False
        detail: list[dict[str, object]] = []
        try:
            with tarfile.open(path, "r") as archive:
                tar_readable = True
                members = [m for m in archive.getmembers() if m.isfile() and m.name.endswith(".gz")]
                for member in sorted(members, key=lambda m: (chromosome_from_member(m.name).zfill(2), m.name)):
                    row = inspect_member(archive, member)
                    row.update({
                        "gene_symbol": spec["gene_symbol"],
                        "UKB_PPP_OID": spec["UKB_PPP_OID"],
                        "archive_name": path.name,
                    })
                    detail.append(row)
                gzip_members_readable = len(detail) == len(members) and all(x["row_count"] > 0 for x in detail)
        except Exception as exc:
            print(f"ERROR {path.name}: {exc}", flush=True)
            members = []

        chromosomes = {str(x["chromosome_label"]) for x in detail}
        expected_chromosomes = {str(i) for i in range(1, 23)} | {"X"}
        all_headers_match = bool(detail) and all(bool(x["header_exact_match"]) for x in detail)
        all_rows_complete = bool(detail) and all(bool(x["first_and_last_complete"]) for x in detail)
        total_rows = sum(int(x["row_count"]) for x in detail)
        manifest_rows.append({
            "analysis_order": spec["analysis_order"],
            "gene_symbol": spec["gene_symbol"],
            "UKB_PPP_OID": spec["UKB_PPP_OID"],
            "UniProt_ID": spec["UniProt_ID"],
            "Olink_panel": spec["Olink_panel"],
            "synapse_id": spec["synapse_id"],
            "download_source": f"Synapse {spec['synapse_id']}; provider-authorized archive verified",
            "archive_name": path.name,
            "archive_path": str(path),
            "file_size_bytes": stat.st_size,
            "modified_time": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
            "md5": md5,
            "sha256": sha256,
            "expected_md5": spec["expected_md5"],
            "md5_matches_frozen_config": md5 == spec["expected_md5"].upper(),
            "expected_sha256": spec["observed_sha256"],
            "sha256_matches_frozen_config": sha256 == spec["observed_sha256"].upper(),
            "member_count": len(members),
            "total_data_rows": total_rows,
            "genome_build_variant_ID": "GRCh37/hg19",
            "GENPOS_build": "GRCh38/hg38",
            "variant_ID_format": "CHR:GRCh37_POS:ALLELE0:ALLELE1",
            "tar_readable": tar_readable,
            "all_gzip_members_readable": gzip_members_readable,
        })
        qa_rows.append({
            "analysis_order": spec["analysis_order"],
            "gene_symbol": spec["gene_symbol"],
            "UKB_PPP_OID": spec["UKB_PPP_OID"],
            "archive_name": path.name,
            "size_matches_frozen_config": stat.st_size == int(spec["expected_size_bytes"]),
            "md5_matches_frozen_config": md5 == spec["expected_md5"].upper(),
            "sha256_matches_frozen_config": sha256 == spec["observed_sha256"].upper(),
            "tar_readable": tar_readable,
            "member_count_is_23": len(members) == 23,
            "chromosomes_1_to_22_and_X_complete": chromosomes == expected_chromosomes,
            "all_headers_exact_match": all_headers_match,
            "all_first_last_rows_complete": all_rows_complete,
            "all_members_have_data": gzip_members_readable,
            "OID_matches_archive_name": spec["UKB_PPP_OID"] in path.name,
            "target_matches_archive_name": spec["gene_symbol"] in path.name,
            "UniProt_matches_archive_name": spec["UniProt_ID"] in path.name,
            "overall_pass": all([
                stat.st_size == int(spec["expected_size_bytes"]),
                md5 == spec["expected_md5"].upper(),
                sha256 == spec["observed_sha256"].upper(),
                tar_readable,
                len(members) == 23,
                chromosomes == expected_chromosomes,
                all_headers_match,
                all_rows_complete,
                gzip_members_readable,
            ]),
            "stale_partial_file_present": (RESOURCE_DIR / f"{path.name}.part").exists(),
            "stale_partial_action": "quarantine_after_complete_archive_pass; never use as input",
        })
        member_rows.extend(detail)

    write_tsv(MANIFEST_OUT, manifest_rows)
    write_tsv(QA_OUT, qa_rows)
    write_tsv(MEMBER_OUT, member_rows)
    if len(manifest_rows) != 5 or not all(bool(row["overall_pass"]) for row in qa_rows):
        raise RuntimeError("PWAS5 archive integrity audit failed; inspect PWAS5_integrity_QA.tsv")
    print(f"PASS: five frozen archives; manifest={MANIFEST_OUT}")


if __name__ == "__main__":
    main()
