#!/usr/bin/env python3
"""Stream the deCODE public XLSX package for eight prespecified target genes."""

from __future__ import annotations

import csv
import io
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data_raw" / "decode_public" / "41586_2023_6563_MOESM3_ESM.xlsx"
JSON_OUT = ROOT / "audit" / "decode_2023_supplement_stream_audit.json"
TSV_OUT = ROOT / "tables" / "decode_2023_target_assay_matches.tsv"
TARGETS = ("CLN5", "COL10A1", "PLOD2", "SDF2", "TMEM106B", "VTN", "BRD2", "IL20RB")

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def text_of(element: ET.Element) -> str:
    return "".join(node.text or "" for node in element.iter() if node.tag.endswith("}t"))


def column_number(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    if not letters:
        return 0
    value = 0
    for char in letters.group(0):
        value = value * 26 + ord(char) - 64
    return value


with zipfile.ZipFile(INPUT) as archive:
    shared_strings: list[str] = []
    with archive.open("xl/sharedStrings.xml") as stream:
        for _, element in ET.iterparse(stream, events=("end",)):
            if element.tag == f"{{{NS_MAIN}}}si":
                shared_strings.append(text_of(element))
                element.clear()

    target_indices: dict[int, list[str]] = {}
    for index, value in enumerate(shared_strings):
        normalized = value.upper()
        hits = [target for target in TARGETS if re.search(rf"(?<![A-Z0-9]){re.escape(target)}(?![A-Z0-9])", normalized)]
        if hits:
            target_indices[index] = hits

    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    rel_map = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in relationships.findall(f"{{{NS_PKG_REL}}}Relationship")
    }
    sheet_map: dict[str, str] = {}
    for sheet in workbook.findall(f".//{{{NS_MAIN}}}sheet"):
        relationship_id = sheet.attrib[f"{{{NS_REL}}}id"]
        target = rel_map[relationship_id].lstrip("/")
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        sheet_map[target] = sheet.attrib["name"]

    index_patterns = [f"<v>{index}</v>".encode() for index in target_indices]
    candidate_sheets: dict[str, bytes] = {}
    for member in archive.infolist():
        if not member.filename.startswith("xl/worksheets/sheet") or not member.filename.endswith(".xml"):
            continue
        data = archive.read(member.filename)
        if any(pattern in data for pattern in index_patterns):
            candidate_sheets[member.filename] = data

    matches: list[dict[str, object]] = []
    sheet_headers: dict[str, list[dict[str, str]]] = {}
    for member_name, data in candidate_sheets.items():
        sheet_name = sheet_map.get(member_name, member_name)
        header_rows: list[dict[str, str]] = []
        for _, row in ET.iterparse(io.BytesIO(data), events=("end",)):
            if row.tag != f"{{{NS_MAIN}}}row":
                continue
            row_number = int(row.attrib.get("r", "0"))
            values: dict[str, str] = {}
            row_targets: set[str] = set()
            for cell in row.findall(f"{{{NS_MAIN}}}c"):
                cell_ref = cell.attrib.get("r", "")
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{NS_MAIN}}}v")
                value = ""
                if cell_type == "s" and value_node is not None and value_node.text is not None:
                    string_index = int(value_node.text)
                    if 0 <= string_index < len(shared_strings):
                        value = shared_strings[string_index]
                    row_targets.update(target_indices.get(string_index, []))
                elif cell_type == "inlineStr":
                    value = text_of(cell)
                elif value_node is not None and value_node.text is not None:
                    value = value_node.text
                values[cell_ref] = value

            ordered_values = {
                ref: values[ref]
                for ref in sorted(values, key=column_number)
            }
            if 0 < row_number <= 12:
                header_rows.append(ordered_values)
            if row_targets:
                matches.append(
                    {
                        "sheet_name": sheet_name,
                        "sheet_xml": member_name,
                        "row_number": row_number,
                        "targets": sorted(row_targets),
                        "cells": ordered_values,
                    }
                )
            row.clear()
        sheet_headers[sheet_name] = header_rows

JSON_OUT.write_text(
    json.dumps(
        {
            "source_file": str(INPUT),
            "source_role": "deCODE 2023 public supplementary workbook; assay mapping audit only",
            "shared_string_count": len(shared_strings),
            "target_shared_string_indices": {str(key): value for key, value in target_indices.items()},
            "candidate_sheets": {name: sheet_map.get(name, name) for name in candidate_sheets},
            "headers": sheet_headers,
            "matches": matches,
        },
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)

with TSV_OUT.open("w", encoding="utf-8", newline="") as handle:
    fields = ["target", "sheet_name", "row_number", "row_cells_json", "source_file", "interpretive_role"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for match in matches:
        for target in match["targets"]:
            writer.writerow(
                {
                    "target": target,
                    "sheet_name": match["sheet_name"],
                    "row_number": match["row_number"],
                    "row_cells_json": json.dumps(match["cells"], ensure_ascii=False, separators=(",", ":")),
                    "source_file": INPUT.name,
                    "interpretive_role": "mapping_audit_only_not_full_GWAS",
                }
            )

print(f"Shared strings: {len(shared_strings)}")
print(f"Target string indices: {len(target_indices)}")
print(f"Candidate sheets: {len(candidate_sheets)}")
print(f"Matched rows: {len(matches)}")
print(f"JSON: {JSON_OUT}")
print(f"TSV: {TSV_OUT}")
