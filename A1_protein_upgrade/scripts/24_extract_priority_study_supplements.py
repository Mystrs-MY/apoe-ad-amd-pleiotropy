#!/usr/bin/env python3
"""Inventory and search the four priority-study supplementary packages."""

from __future__ import annotations

import hashlib
import re
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "literature" / "fulltext_supplements"
OUT = ROOT / "data_processed" / "priority_study_supplements"
LOG = ROOT / "logs"

TARGETS = {
    "CLN5": [r"\bCLN5\b", r"ceroid[- ]lipofuscinosis.*5"],
    "COL10A1": [r"\bCOL10A1\b", r"collagen.*type\s*X.*alpha\s*1"],
    "PLOD2": [r"\bPLOD2\b", r"lysyl hydroxylase\s*2"],
    "SDF2": [r"\bSDF2\b", r"stromal cell[- ]derived factor\s*2"],
    "TMEM106B": [r"\bTMEM106B\b", r"transmembrane protein\s*106B"],
    "VTN": [r"\bVTN\b", r"\bvitronectin\b", r"serum spreading factor"],
    "TREM2": [r"\bTREM2\b", r"triggering receptor expressed on myeloid cells\s*2"],
    "ACE": [r"\bACE\b", r"angiotensin[- ]converting enzyme"],
    "LBP": [r"\bLBP\b", r"lipopolysaccharide[- ]binding protein"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def study_id(path: Path) -> str:
    match = re.search(r"PMID(\d+)", str(path))
    return match.group(1) if match else "unresolved"


def unpack_science_zip() -> Path:
    archive = SOURCE / "PMID42384774_tables_s1_to_s30.zip"
    destination = OUT / "PMID42384774"
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as package:
        for member in package.infolist():
            target = (destination / member.filename).resolve()
            if destination.resolve() not in target.parents and target != destination.resolve():
                raise RuntimeError(f"Unsafe archive member: {member.filename}")
            package.extract(member, destination)
    files = sorted(destination.glob("*.xlsx"))
    if len(files) != 1:
        raise RuntimeError("Expected exactly one Science supplementary workbook")
    return files[0]


def read_sheet(path: Path, sheet: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, dtype=object, engine="openpyxl")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    science_book = unpack_science_zip()

    workbooks = sorted(SOURCE.glob("PMID*.xlsx")) + [science_book]
    inventory_rows: list[dict[str, object]] = []
    hit_rows: list[dict[str, object]] = []

    for workbook in workbooks:
        sid = study_id(workbook)
        excel = pd.ExcelFile(workbook, engine="openpyxl")
        for sheet in excel.sheet_names:
            try:
                frame = read_sheet(workbook, sheet)
            except Exception as exc:  # preserve unreadable sheets in the audit trail
                inventory_rows.append(
                    {
                        "PMID": sid,
                        "file": str(workbook),
                        "sheet": sheet,
                        "rows": None,
                        "columns": None,
                        "column_names": None,
                        "read_status": f"error: {type(exc).__name__}: {exc}",
                    }
                )
                continue

            inventory_rows.append(
                {
                    "PMID": sid,
                    "file": str(workbook),
                    "sheet": sheet,
                    "rows": len(frame),
                    "columns": len(frame.columns),
                    "column_names": " | ".join(map(str, frame.columns)),
                    "read_status": "ok",
                }
            )
            if frame.empty:
                continue

            joined = frame.fillna("").astype(str).agg(" | ".join, axis=1)
            for gene, patterns in TARGETS.items():
                mask = pd.Series(False, index=frame.index)
                for pattern in patterns:
                    mask |= joined.str.contains(pattern, case=False, regex=True, na=False)
                for row_index in frame.index[mask]:
                    row = frame.loc[row_index]
                    nonempty = {
                        str(column): str(value)
                        for column, value in row.items()
                        if pd.notna(value) and str(value).strip() not in {"", "nan"}
                    }
                    hit_rows.append(
                        {
                            "PMID": sid,
                            "target_gene": gene,
                            "file": str(workbook),
                            "sheet": sheet,
                            "excel_row": int(row_index) + 2,
                            "matched_row": " | ".join(
                                f"{key}={value}" for key, value in nonempty.items()
                            )[:12000],
                        }
                    )

    inventory = pd.DataFrame(inventory_rows)
    hits = pd.DataFrame(hit_rows).drop_duplicates()
    inventory.to_csv(OUT / "priority_study_sheet_inventory.tsv", sep="\t", index=False)
    hits.to_csv(OUT / "priority_study_target_hits.tsv", sep="\t", index=False)

    manifest_rows = []
    for path in sorted(SOURCE.glob("PMID*")):
        if path.is_file() and path.suffix.lower() in {".xlsx", ".docx", ".pdf", ".zip"}:
            manifest_rows.append(
                {
                    "PMID": study_id(path),
                    "file": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(manifest_rows).to_csv(
        LOG / "priority_study_download_manifest.tsv", sep="\t", index=False
    )

    expected = {"34381170", "40397384", "40452368", "42384774"}
    observed = set(pd.DataFrame(manifest_rows)["PMID"])
    if not expected.issubset(observed):
        raise RuntimeError(f"Missing priority-study packages: {sorted(expected - observed)}")
    print(f"Inventoried {len(inventory):,} sheets; retained {len(hits):,} target-hit rows")


if __name__ == "__main__":
    main()
