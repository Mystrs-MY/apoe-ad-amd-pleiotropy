#!/usr/bin/env python3
"""Inventory tabular literature supplements without modifying source files."""

from __future__ import annotations

import csv
import json
import re
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / "data_raw" / "literature_fulltext"
OUTPUT_TSV = ROOT / "audit" / "literature_supplement_table_inventory.tsv"
OUTPUT_JSON = ROOT / "audit" / "literature_supplement_table_inventory.json"


def pmid_from_path(path: Path) -> str:
    for part in path.parts:
        match = re.match(r"PMID(\d+)", part)
        if match:
            return f"PMID{match.group(1)}"
    return "unresolved"


def clean_value(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()[:300]


def workbook_rows(path: Path, source_container: str = "") -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    try:
        workbook = pd.ExcelFile(path)
        for sheet in workbook.sheet_names:
            try:
                frame = pd.read_excel(path, sheet_name=sheet, header=None, nrows=8)
                preview = [[clean_value(value) for value in row] for row in frame.values.tolist()]
                nonempty = frame.dropna(how="all").dropna(axis=1, how="all")
                rows.append(
                    {
                        "PMID": pmid_from_path(path),
                        "source_file": str(path.relative_to(ROOT)),
                        "source_container": source_container,
                        "file_type": path.suffix.lower(),
                        "sheet_or_member": sheet,
                        "reported_rows": workbook.book[sheet].max_row,
                        "reported_columns": workbook.book[sheet].max_column,
                        "preview_nonempty_rows": len(nonempty.index),
                        "preview_nonempty_columns": len(nonempty.columns),
                        "preview_json": json.dumps(preview, ensure_ascii=False),
                        "read_status": "readable",
                        "error": "",
                    }
                )
            except Exception as exc:  # preserve sheet-level failures for audit
                rows.append(error_row(path, sheet, source_container, exc))
    except Exception as exc:
        rows.append(error_row(path, "", source_container, exc))
    return rows


def delimited_row(path: Path) -> dict[str, object]:
    try:
        frame = pd.read_csv(path, sep=None, engine="python", nrows=8)
        preview = [[clean_value(value) for value in frame.columns.tolist()]]
        preview.extend([[clean_value(value) for value in row] for row in frame.values.tolist()])
        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            line_count = sum(1 for _ in handle)
        return {
            "PMID": pmid_from_path(path),
            "source_file": str(path.relative_to(ROOT)),
            "source_container": "",
            "file_type": path.suffix.lower(),
            "sheet_or_member": "",
            "reported_rows": max(line_count - 1, 0),
            "reported_columns": len(frame.columns),
            "preview_nonempty_rows": len(frame.index),
            "preview_nonempty_columns": len(frame.columns),
            "preview_json": json.dumps(preview, ensure_ascii=False),
            "read_status": "readable",
            "error": "",
        }
    except Exception as exc:
        return error_row(path, "", "", exc)


def error_row(path: Path, member: str, container: str, exc: Exception) -> dict[str, object]:
    return {
        "PMID": pmid_from_path(path),
        "source_file": str(path.relative_to(ROOT)),
        "source_container": container,
        "file_type": path.suffix.lower(),
        "sheet_or_member": member,
        "reported_rows": "",
        "reported_columns": "",
        "preview_nonempty_rows": "",
        "preview_nonempty_columns": "",
        "preview_json": "[]",
        "read_status": "unreadable",
        "error": f"{type(exc).__name__}: {exc}"[:500],
    }


def extract_zip_tables(path: Path) -> list[Path]:
    destination = ROOT / "data_processed" / "literature_supplements_extracted" / pmid_from_path(path) / path.stem
    destination.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            member_path = Path(info.filename)
            if info.is_dir() or member_path.suffix.lower() not in {".xlsx", ".xls", ".csv", ".tsv"}:
                continue
            safe_name = "__".join(member_path.parts)
            target = destination / safe_name
            resolved = target.resolve()
            if destination.resolve() not in resolved.parents:
                raise ValueError(f"Unsafe ZIP member: {info.filename}")
            if not target.exists() or target.stat().st_size != info.file_size:
                with archive.open(info) as source, target.open("wb") as output:
                    output.write(source.read())
            extracted.append(target)
    return extracted


def main() -> None:
    rows: list[dict[str, object]] = []
    source_tables = sorted(
        path
        for path in SOURCE_ROOT.rglob("*")
        if path.is_file() and path.suffix.lower() in {".xlsx", ".xls", ".csv", ".tsv"}
    )
    for path in source_tables:
        if path.suffix.lower() in {".xlsx", ".xls"}:
            rows.extend(workbook_rows(path))
        else:
            rows.append(delimited_row(path))

    for archive in sorted(SOURCE_ROOT.rglob("*.zip")):
        try:
            for extracted in extract_zip_tables(archive):
                container = str(archive.relative_to(ROOT))
                if extracted.suffix.lower() in {".xlsx", ".xls"}:
                    rows.extend(workbook_rows(extracted, container))
                else:
                    row = delimited_row(extracted)
                    row["source_container"] = container
                    rows.append(row)
        except Exception as exc:
            rows.append(error_row(archive, "", "", exc))

    deduplicated: list[dict[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for row in rows:
        key = (
            row["PMID"],
            Path(str(row["source_file"])).name.lower(),
            row["sheet_or_member"],
            row["reported_rows"],
            row["reported_columns"],
            row["preview_json"],
        )
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(row)
    rows = deduplicated

    columns = [
        "PMID",
        "source_file",
        "source_container",
        "file_type",
        "sheet_or_member",
        "reported_rows",
        "reported_columns",
        "preview_nonempty_rows",
        "preview_nonempty_columns",
        "preview_json",
        "read_status",
        "error",
    ]
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    OUTPUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Inventoried {len(source_tables)} source tables into {len(rows)} unique sheet/table records")
    print(f"Output: {OUTPUT_TSV}")


if __name__ == "__main__":
    main()
