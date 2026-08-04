#!/usr/bin/env python3
"""Extract and harmonize direct APOE variant effects from current A1 outcome GWAS files."""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
OUTCOMES_CONFIG = yaml.safe_load((ROOT / "config" / "outcomes.yml").read_text(encoding="utf-8"))
OUTPUT = ROOT / "tables" / "APOE_variant_total_effects_current_A1.tsv"

VARIANTS = {
    "rs429358": {"position_hg19": 45411941, "effect_allele": "C", "other_allele": "T"},
    "rs7412": {"position_hg19": 45412079, "effect_allele": "T", "other_allele": "C"},
}


def read_rows(path: Path) -> dict[str, dict[str, str]]:
    found: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chromosome = str(row["CHR"]).replace("chr", "")
            if chromosome != "19":
                continue
            position = int(float(row["BP"]))
            for rsid, specification in VARIANTS.items():
                if row["SNP"] == rsid or position == specification["position_hg19"]:
                    found[rsid] = row
            if len(found) == len(VARIANTS):
                break
    return found


def main() -> None:
    rows = []
    for outcome, specification in OUTCOMES_CONFIG["outcomes"].items():
        path = Path(specification["local_file"])
        found = read_rows(path)
        for rsid, variant in VARIANTS.items():
            source = found.get(rsid)
            if source is None:
                rows.append({
                    "variant": rsid, "outcome": outcome, "position_hg19": variant["position_hg19"],
                    "requested_effect_allele": variant["effect_allele"], "requested_other_allele": variant["other_allele"],
                    "source_SNP": "NA", "original_effect_allele": "NA", "original_other_allele": "NA",
                    "allele_flipped": "NA", "beta": "NA", "SE": "NA", "P_value": "NA", "EAF_original": "NA",
                    "sample_size": "NA", "outcome_GWAS": specification["current_A1_source"],
                    "source_file": str(path), "availability_status": "direct_variant_not_found",
                    "exclusion_reason": "Variant absent by both rsID and verified hg19 position.",
                })
                continue
            effect = source["A1"].upper()
            other = source["A2"].upper()
            source_beta = float(source["BETA"])
            source_eaf = float(source["FREQ"])
            if effect == variant["effect_allele"] and other == variant["other_allele"]:
                beta = source_beta
                eaf = source_eaf
                flipped = "false"
            elif other == variant["effect_allele"] and effect == variant["other_allele"]:
                beta = -source_beta
                eaf = 1 - source_eaf
                flipped = "true"
            else:
                beta = eaf = math.nan
                flipped = "not_harmonizable"
            rows.append({
                "variant": rsid, "outcome": outcome, "position_hg19": variant["position_hg19"],
                "requested_effect_allele": variant["effect_allele"], "requested_other_allele": variant["other_allele"],
                "source_SNP": source["SNP"], "original_effect_allele": effect, "original_other_allele": other,
                "allele_flipped": flipped, "beta": beta, "SE": source["SE"], "P_value": source["P"],
                "EAF_harmonized": eaf, "EAF_original": source_eaf, "sample_size": source["N"],
                "outcome_GWAS": specification["current_A1_source"], "source_file": str(path),
                "availability_status": "direct_variant_available" if flipped != "not_harmonizable" else "alleles_not_harmonizable",
                "exclusion_reason": "NA" if flipped != "not_harmonizable" else f"Observed {effect}/{other}.",
            })
    output = pd.DataFrame(rows).replace({math.nan: "NA"})
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)
    print(f"Direct total-effect rows: {(output['availability_status'] == 'direct_variant_available').sum()}/{len(output)}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
