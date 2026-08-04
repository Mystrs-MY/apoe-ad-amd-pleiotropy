"""High-recall, auditable title/abstract pre-screen for literature records."""

from __future__ import annotations

import csv
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "literature" / "search_results_deduplicated.tsv"
OUTPUT = ROOT / "literature" / "study_screening.tsv"
KEYWORDS = ROOT / "config" / "screening_keywords.yml"


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def any_term(text: str, terms: list[str]) -> bool:
    return any(term.lower() in text for term in terms)


def main() -> int:
    records = read_tsv(INPUT)
    keywords = yaml.safe_load(KEYWORDS.read_text(encoding="utf-8"))
    fields = [
        "record_id",
        "title",
        "first_author",
        "publication_year",
        "journal",
        "PMID",
        "DOI",
        "publication_status",
        "source_databases",
        "disease_scope",
        "circulating_protein_scope",
        "proteome_scale",
        "genetic_causal_method",
        "title_abstract_decision",
        "title_abstract_exclusion_reason",
        "full_text_status",
        "full_text_decision",
        "full_text_exclusion_reason",
        "supplement_status",
        "duplicate_group",
        "evidence_independence",
        "manual_verification_status",
        "notes",
    ]

    counts: dict[str, int] = {}
    output_rows = []
    for record in records:
        text = f"{record.get('title', '')} {record.get('abstract', '')}".lower()
        disease_hits = [
            name for name, terms in keywords["disease"].items() if any_term(text, terms)
        ]
        protein_hit = any_term(text, keywords["circulating_protein"])
        causal_hit = any_term(text, keywords["genetic_causal"])
        background_hit = any_term(record.get("title", "").lower(), keywords["background_only"])

        if disease_hits and protein_hit and causal_hit and not background_hit:
            decision = "candidate_full_text"
            exclusion = ""
        elif disease_hits and (protein_hit or causal_hit) and not record.get("abstract", ""):
            decision = "manual_review_missing_abstract"
            exclusion = ""
        else:
            decision = "exclude_title_abstract"
            reasons = []
            if not disease_hits:
                reasons.append("disease_scope_not_confirmed")
            if not protein_hit:
                reasons.append("circulating_or_large_scale_protein_scope_not_confirmed")
            if not causal_hit:
                reasons.append("genetic_causal_method_not_confirmed")
            if background_hit:
                reasons.append("review_protocol_or_commentary")
            exclusion = ";".join(reasons)

        counts[decision] = counts.get(decision, 0) + 1
        output_rows.append(
            {
                "record_id": record.get("record_id", ""),
                "title": record.get("title", ""),
                "first_author": record.get("first_author", ""),
                "publication_year": record.get("publication_year", ""),
                "journal": record.get("journal", ""),
                "PMID": record.get("PMID", ""),
                "DOI": record.get("DOI", ""),
                "publication_status": record.get("publication_status", ""),
                "source_databases": record.get("source_databases", ""),
                "disease_scope": ";".join(disease_hits) if disease_hits else "not_confirmed",
                "circulating_protein_scope": "candidate" if protein_hit else "not_confirmed",
                "proteome_scale": "requires_full_text_verification" if protein_hit else "not_confirmed",
                "genetic_causal_method": "candidate" if causal_hit else "not_confirmed",
                "title_abstract_decision": decision,
                "title_abstract_exclusion_reason": exclusion,
                "full_text_status": "not_assessed",
                "full_text_decision": "",
                "full_text_exclusion_reason": "",
                "supplement_status": "not_assessed",
                "duplicate_group": "",
                "evidence_independence": "not_assessed",
                "manual_verification_status": "pending",
                "notes": "Automated high-recall pre-screen; no evidence tier assigned.",
            }
        )

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    for key in sorted(counts):
        print(f"{key}: {counts[key]}")
    print(f"Output: {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
