#!/usr/bin/env python3
"""Fail closed when a candidate public release contains unsafe or stale assets."""

from __future__ import annotations

import csv
import math
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 50 * 1024 * 1024
RELEASE_VERSION = "1.0.0"
RELEASE_DATE = "2026-08-26"
CORRECTED_RUN_TAG = "correctedN_20260903"
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
    "excluded plotting or figure-export script": re.compile(
        r"(?:^|/)[^/]*(?:plot|figure|forest|visual|redraw|export.*asset)[^/]*\.(?:r|py|ps1|sh|bat)$",
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
    "local WSL mount path": re.compile(r"/mnt/[a-z]/(?:Article_1|AD_AMD|Users)/", re.I),
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
        "README.md", "CHANGELOG.md", "LICENSE", "CITATION.cff", "NOTICE.md",
        "environment/release_environment_1.0.0.txt",
        "02_genetic_arch/LAVA/0_prepare_LAVA_totalN_inputs.R",
        "02_genetic_arch/LAVA/1.LAVA.R",
        "02_genetic_arch/MiXeR/prepare_mixer_inputs.R",
        "02_genetic_arch/MiXeR/run_MiXeR_corrected.sh",
        "02_genetic_arch/LDSC/LDSC_Results_Formatted.csv",
        "02_genetic_arch/HDL/Figure_3_HDL_results.csv",
        "02_genetic_arch/MiXeR/MiXeR_bivariate_results.csv",
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
            "doi": "10.5281/zenodo.22102006",
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
        if "10.5281/zenodo.22102006" not in readme:
            failures.append("README.md does not report the version-specific Zenodo DOI")
        if "current `main` branch contains a post-v1.0 technical correction" not in readme:
            failures.append("README.md does not distinguish corrected main from immutable v1.0.0")

    def read_csv_rows(relative: str) -> list[dict[str, str]]:
        with (ROOT / relative).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))

    def assert_close(label: str, observed: float, expected: float, tolerance: float = 1e-10) -> None:
        if not math.isclose(observed, expected, rel_tol=tolerance, abs_tol=tolerance):
            failures.append(f"{label} is {observed}; expected {expected}")

    expected_ldsc = {
        "Dry_AMD": (-0.111368922314042, 0.160119852732122),
        "Wet_AMD": (0.0302085075153177, 0.140511950225366),
        "Any_AMD": (0.00411709827365306, 0.155753217266933),
    }
    ldsc_rows = read_csv_rows("02_genetic_arch/LDSC/LDSC_Results_Formatted.csv")
    ldsc_ad = {row["p2"]: row for row in ldsc_rows if row["p1"] == "AD"}
    for outcome, (expected_rg, expected_se) in expected_ldsc.items():
        row = ldsc_ad.get(outcome)
        if row is None:
            failures.append(f"corrected LDSC AD-{outcome} row is missing")
            continue
        if row.get("analysis_run_tag") != CORRECTED_RUN_TAG:
            failures.append(f"LDSC AD-{outcome} run tag is stale")
        assert_close(f"LDSC AD-{outcome} rg", float(row["rg"]), expected_rg)
        assert_close(f"LDSC AD-{outcome} SE", float(row["se"]), expected_se)

    expected_hdl = {
        "Dry_AMD": (0.0593806219831374, 0.100275317783007),
        "Wet_AMD": (-0.0321337347871692, 0.0642618005248076),
        "Any_AMD": (0.0209624102170959, 0.0736883169478447),
    }
    hdl_rows = read_csv_rows("02_genetic_arch/HDL/Figure_3_HDL_results.csv")
    hdl_ad = {row["Trait2"]: row for row in hdl_rows if row["Trait1"] == "AD"}
    for outcome, (expected_rg, expected_se) in expected_hdl.items():
        row = hdl_ad.get(outcome)
        if row is None:
            failures.append(f"corrected HDL AD-{outcome} row is missing")
            continue
        if row.get("analysis_run_tag") != CORRECTED_RUN_TAG:
            failures.append(f"HDL AD-{outcome} run tag is stale")
        assert_close(f"HDL AD-{outcome} rg", float(row["rg"]), expected_rg)
        assert_close(f"HDL AD-{outcome} SE", float(row["se"]), expected_se)

    expected_mixer = {
        "AMD_Dry": (52.5924987238301, -0.651768119418754),
        "AMD_Wet": (2.47775918736052, -0.321797019048213),
        "AMD_Any": (20.9145657945658, -0.852595482028776),
    }
    mixer_rows = read_csv_rows("02_genetic_arch/MiXeR/MiXeR_bivariate_results.csv")
    mixer = {row["Trait2"]: row for row in mixer_rows}
    for outcome, (expected_overlap, expected_rho) in expected_mixer.items():
        row = mixer.get(outcome)
        if row is None:
            failures.append(f"corrected MiXeR {outcome} row is missing")
            continue
        if row.get("analysis_run_tag") != CORRECTED_RUN_TAG:
            failures.append(f"MiXeR {outcome} run tag is stale")
        assert_close(f"MiXeR {outcome} overlap", float(row["Overlap_pct"]), expected_overlap)
        assert_close(f"MiXeR {outcome} rho_beta", float(row["rho_beta"]), expected_rho)

    lava_files = {
        "02_genetic_arch/LAVA/LAVA_FullScan_AD_vs_DryAMD_Final.csv": (
            2029, -0.319401, 8.02075e-14
        ),
        "02_genetic_arch/LAVA/LAVA_FullScan_AD_vs_WetAMD_Final.csv": (
            2036, -0.192170, 8.68393e-7
        ),
        "02_genetic_arch/LAVA/LAVA_FullScan_AD_vs_AnyAMD_Final.csv": (
            2086, -0.259934, 9.23234e-13
        ),
    }
    lava_supported = []
    lava_threshold = 0.05 / 2495
    for relative, (expected_estimable_rows, expected_rho, expected_p) in lava_files.items():
        rows = read_csv_rows(relative)
        if len(rows) != expected_estimable_rows:
            failures.append(
                f"{relative} has {len(rows)} estimable rows; expected {expected_estimable_rows}"
            )
        supported_here = []
        for row in rows:
            if row.get("p") not in (None, "", "NA") and float(row["p"]) < lava_threshold:
                lava_supported.append((relative, int(float(row["locus"]))))
                supported_here.append(row)
        if len(supported_here) == 1:
            assert_close(
                f"LAVA {relative} APOE rho", float(supported_here[0]["rho"]),
                expected_rho, tolerance=1e-6
            )
            assert_close(
                f"LAVA {relative} APOE P", float(supported_here[0]["p"]),
                expected_p, tolerance=1e-18
            )
    if len(lava_supported) != 3 or any(locus != 2351 for _, locus in lava_supported):
        failures.append(
            "corrected LAVA support must contain exactly three rows, all at locus 2351"
        )

    magma_path = ROOT / "02_genetic_arch/MAGMA/source_data/FigS2_MAGMA_source_data.tsv"
    with magma_path.open(encoding="utf-8-sig", newline="") as handle:
        magma_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(magma_rows) != 4 or {row.get("locus") for row in magma_rows} != {"APOE"}:
        failures.append("corrected MAGMA locus summary must contain four APOE rows")
    expected_best_genes = {
        "AD": "BCL3", "Dry AMD": "APOC1", "Wet AMD": "APOC1", "Any AMD": "APOC1"
    }
    for row in magma_rows:
        if row.get("best_gene") != expected_best_genes.get(row.get("trait")):
            failures.append(f"unexpected corrected MAGMA best gene for {row.get('trait')}")

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
