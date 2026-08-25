#!/usr/bin/env python3
"""Fail closed when a candidate public release contains unsafe or stale assets."""

from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 50 * 1024 * 1024
RELEASE_VERSION = "1.0.0"
RELEASE_DATE = "2026-08-26"
SKIP_DIRECTORIES = {
    ".git", "__pycache__", "main", "previews", "supplementary_figures",
    "_Figure_4_panels", "grouped_supplementary_panel_assets",
}

FORBIDDEN_SUFFIXES = {
    ".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam", ".rds", ".rdata",
    ".docx", ".pptx", ".xlsx", ".xls", ".tar", ".gz", ".zip",
}
ALLOWED_BINARY: set[Path] = set()
FORBIDDEN_RELATIVE_PATHS = {
    Path("05_mediation/09_figures_main.R"),
    Path("01_meta/PRISMA/build_prisma_documents.ps1"),
    Path("01_meta/PRISMA/template_contract.md"),
    Path("figures_submission/code/Figure_1_study_design.R"),
    Path("figures_submission/code/FigS8_decode_extension.R"),
    Path("figures_submission/code/Figure_3_vector.R"),
    Path("figures_submission/code/Figure_4_vector.R"),
    Path("figures_submission/code/Figure_5_name_matched_assembled.R"),
    Path("figures_submission/code/compose_Figure_4_pdf.py"),
    Path("figures_submission/code/compose_supplementary_FigS3_FigS7.R"),
    Path("figures_submission/code/figure_utils.R"),
    Path("02_genetic_arch/MiXeR/plot_mixer_venn.R"),
    Path("02_genetic_arch/LAVA/LAVAFigure_Updated.R"),
    Path("02_genetic_arch/HyPrColoc_v1/LocusZoom_Article1.R"),
    Path("03_causal_lock/12_fig_apoe_scatter.R"),
    Path("P0_finemap/05_figS6_gcta_cojo.R"),
    Path("A1_protein_upgrade/scripts/18_plot_independent_protein_figures.R"),
    Path("A1_protein_upgrade/scripts/19_qa_independent_protein_figures.py"),
    Path("figures_submission/assets/Fig4b_LAVA_Dual_Panel.pdf"),
}
FORBIDDEN_PATH_PATTERNS = {
    "excluded diagram or graphical-abstract asset": re.compile(
        r"^figures_submission/(?:code|assets|source_data)/.*(?:graphical[_ -]?abstract|figure[_ -]?5a)",
        re.I,
    ),
    "excluded multi-panel assembly code": re.compile(
        r"^figures_submission/code/(?:compose_|.*assembled|compose_supplementary)",
        re.I,
    ),
}
TEXT_SUFFIXES = {
    ".r", ".py", ".ps1", ".sh", ".awk", ".md", ".txt", ".yml", ".yaml",
    ".json", ".csv", ".tsv", ".cff", ".out", ".cojo", "",
}
SECRET_PATTERNS = {
    "JWT": re.compile("eyJ0eXAi" + "OiJKV1Qi"),
    "private key": re.compile(r"BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY"),
    "AWS secret assignment": re.compile(r"AWS_SECRET_ACCESS_KEY\s*[=:]\s*[^<$\s]", re.I),
}
MACHINE_PATH_PATTERNS = {
    "local Article_1 path": re.compile(r"[A-Za-z]:[\\/]Article_1", re.I),
    "local user path": re.compile(r"[A-Za-z]:[\\/]Users[\\/]", re.I),
    "deprecated resource path": re.compile(r"[A-Za-z]:[\\/]AD_AMD[\\/]Resource", re.I),
}


def main() -> int:
    failures: list[str] = []
    files = []
    listed = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    candidates = sorted({entry.decode("utf-8") for entry in listed.split(b"\0") if entry})
    for relative_name in candidates:
        path = ROOT / relative_name
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if relative.parts and relative.parts[0] == "figures":
            continue
        if any(part in SKIP_DIRECTORIES for part in relative.parts):
            continue
        if relative.parts[:2] == ("figures_submission", "code") and (
            path.name.endswith("_export_QA.tsv") or path.name.endswith("_sessionInfo.txt")
        ):
            continue
        if relative.parts[:1] == ("figures_submission",) and path.name.endswith("_vector_QA.csv"):
            continue
        files.append(path)

    for path in files:
        relative = path.relative_to(ROOT)
        if relative in FORBIDDEN_RELATIVE_PATHS:
            failures.append(f"excluded diagram-layout asset present: {relative}")
        relative_text = relative.as_posix()
        for label, pattern in FORBIDDEN_PATH_PATTERNS.items():
            if pattern.search(relative_text):
                failures.append(f"{label}: {relative}")
        if path.stat().st_size > MAX_FILE_BYTES:
            failures.append(f"oversized file: {relative} ({path.stat().st_size} bytes)")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES and relative not in ALLOWED_BINARY:
            failures.append(f"forbidden raw/binary suffix: {relative}")
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8-sig")
        except UnicodeDecodeError:
            failures.append(f"non-UTF-8 text candidate: {relative}")
            continue
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                failures.append(f"{label} pattern: {relative}")
        for label, pattern in MACHINE_PATH_PATTERNS.items():
            if pattern.search(text):
                failures.append(f"{label}: {relative}")

    required = [
        "README.md", "LICENSE", "CITATION.cff", "NOTICE.md",
        "environment/release_environment_1.0.0.txt",
        "A1_protein_upgrade/tables/Table_Literature_Prioritized_Protein_Provenance.tsv",
        "A1_protein_upgrade/tables/APOE_linkable_two_step_mediation.tsv",
        "A1_dual_scale_mediation/tables/decode_smp_two_step_mediation.tsv",
        "figures_submission/source_data/Figure_2_source_data.csv",
        "figures_submission/source_data/FigS6_GCTA_COJO_source_data.tsv",
    ]
    for relative in required:
        if not (ROOT / relative).exists():
            failures.append(f"required release asset missing: {relative}")

    citation_path = ROOT / "CITATION.cff"
    if citation_path.exists():
        citation = citation_path.read_text(encoding="utf-8-sig")
        expected_citation_fields = {
            "version": RELEASE_VERSION,
            "date-released": RELEASE_DATE,
            "repository-code": "https://github.com/Mystrs-MY/apoe-ad-amd-pleiotropy",
            "repository-artifact": (
                "https://github.com/Mystrs-MY/apoe-ad-amd-pleiotropy/releases/tag/v1.0.0"
            ),
            "license": "MIT",
        }
        for field, expected in expected_citation_fields.items():
            pattern = re.compile(
                rf"(?m)^{re.escape(field)}:\s*[\"']?{re.escape(expected)}[\"']?\s*$"
            )
            if not pattern.search(citation):
                failures.append(
                    f"CITATION.cff {field} does not match release value {expected}"
                )

    readme_path = ROOT / "README.md"
    if readme_path.exists():
        readme = readme_path.read_text(encoding="utf-8-sig")
        if "public version 1.0.0 reproducibility release" not in readme:
            failures.append("README.md does not declare the public v1.0.0 release")
        if "pre-submission private repository" in readme:
            failures.append("README.md retains the obsolete private-repository status")

    decode_boundary_files = [
        ROOT / "docs/data_access.md",
        ROOT / "A1_dual_scale_mediation/config/decode_extension_parameters.tsv",
    ]
    for path in decode_boundary_files:
        if path.exists():
            boundary_text = path.read_text(encoding="utf-8-sig")
            if re.search(r"Independent same-platform", boundary_text, re.I):
                failures.append(
                    f"deCODE extension is incorrectly labelled independent: {path.relative_to(ROOT)}"
                )

    provenance_path = ROOT / "A1_protein_upgrade/tables/Table_Literature_Prioritized_Protein_Provenance.tsv"
    if provenance_path.exists():
        with provenance_path.open(encoding="utf-8-sig", newline="") as handle:
            provenance = list(csv.DictReader(handle, delimiter="\t"))
        tier1 = sum(row.get("evidence_tier") == "Tier1" for row in provenance)
        tier2 = sum(row.get("evidence_tier") == "Tier2" for row in provenance)
        priority_pmids = {"34381170", "40397384", "40452368", "42384774"}
        priority_rows = sum(row.get("PMID") in priority_pmids for row in provenance)
        record_ids = [row.get("record_id", "") for row in provenance]
        if len(provenance) != 345:
            failures.append(f"provenance row count is {len(provenance)}; expected 345")
        if (tier1, tier2) != (52, 293):
            failures.append(f"provenance tiers are Tier1={tier1}, Tier2={tier2}; expected 52 and 293")
        if priority_rows != 61:
            failures.append(f"priority-study provenance subset is {priority_rows}; expected 61")
        if len(set(record_ids)) != len(record_ids) or not all(record_ids):
            failures.append("provenance record_id values must be non-empty and unique")

    if failures:
        print("Release validation FAILED")
        for failure in sorted(set(failures)):
            print(f"- {failure}")
        return 1

    total_bytes = sum(path.stat().st_size for path in files)
    print(f"Release validation passed: {len(files)} files, {total_bytes / 1024 / 1024:.2f} MiB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
