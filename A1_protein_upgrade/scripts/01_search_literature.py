"""Reproducible PubMed and CrossRef discovery search for the A1 protein upgrade."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit("PyYAML is required: python -m pip install pyyaml") from exc


ROOT = Path(__file__).resolve().parents[1]
QUERY_FILE = ROOT / "literature" / "search_queries.yml"
RAW_DIR = ROOT / "data_raw" / "literature_search"
OUT_FILE = ROOT / "literature" / "search_results_deduplicated.tsv"
LOG_FILE = ROOT / "logs" / "literature_search_run.json"


def request_text(url: str, user_agent: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8")


def normalize_doi(value: str) -> str:
    value = html.unescape(value or "").strip().lower()
    value = re.sub(r"^https?://(dx\.)?doi\.org/", "", value)
    return value.rstrip(" .;,)")


def normalize_title(value: str) -> list[str]:
    stop = {"a", "an", "the", "in", "of", "for", "on", "to", "and", "with", "by", "et", "al"}
    cleaned = re.sub(r"[^a-z0-9 ]+", " ", html.unescape(value or "").lower())
    return [token for token in cleaned.split() if token and token not in stop]


def title_similarity(a: str, b: str) -> float:
    aa, bb = set(normalize_title(a)), set(normalize_title(b))
    union = aa | bb
    return len(aa & bb) / len(union) if union else 0.0


def first_author_surname(value: str) -> str:
    value = (value or "").strip().lower()
    return re.split(r"[\s,]+", value)[0] if value else ""


def pubmed_search(name: str, query: str, config: dict, user_agent: str) -> tuple[list[dict], dict]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": "10000",
        "tool": config.get("ncbi_tool", "A1_protein_upgrade"),
    }
    if config.get("contact_email"):
        params["email"] = config["contact_email"]
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    search_raw = request_text(search_url, user_agent)
    search_json = json.loads(search_raw)
    ids = search_json["esearchresult"].get("idlist", [])
    count = int(search_json["esearchresult"].get("count", 0))

    records: list[dict] = []
    for start in range(0, len(ids), 200):
        chunk = ids[start : start + 200]
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(chunk),
            "retmode": "xml",
            "tool": config.get("ncbi_tool", "A1_protein_upgrade"),
        }
        if config.get("contact_email"):
            fetch_params["email"] = config["contact_email"]
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(fetch_params)
        xml_text = request_text(fetch_url, user_agent)
        root = ET.fromstring(xml_text)
        for article in root.findall(".//PubmedArticle"):
            medline = article.find("MedlineCitation")
            citation = medline.find("Article") if medline is not None else None
            if medline is None or citation is None:
                continue
            pmid = medline.findtext("PMID", default="")
            title = "".join(citation.find("ArticleTitle").itertext()) if citation.find("ArticleTitle") is not None else ""
            abstract = " ".join(
                "".join(node.itertext()) for node in citation.findall("Abstract/AbstractText")
            )
            authors = []
            for author in citation.findall("AuthorList/Author"):
                collective = author.findtext("CollectiveName", default="")
                if collective:
                    authors.append(collective)
                else:
                    surname = author.findtext("LastName", default="")
                    initials = author.findtext("Initials", default="")
                    authors.append(" ".join(x for x in (surname, initials) if x))
            journal = citation.findtext("Journal/Title", default="")
            year = (
                citation.findtext("Journal/JournalIssue/PubDate/Year", default="")
                or citation.findtext("Journal/JournalIssue/PubDate/MedlineDate", default="")[:4]
            )
            doi = ""
            for aid in article.findall("PubmedData/ArticleIdList/ArticleId"):
                if aid.attrib.get("IdType") == "doi":
                    doi = aid.text or ""
                    break
            pub_types = [node.text or "" for node in citation.findall("PublicationTypeList/PublicationType")]
            records.append(
                {
                    "title": title,
                    "first_author": authors[0] if authors else "",
                    "authors": "; ".join(authors),
                    "publication_year": year,
                    "journal": journal,
                    "PMID": pmid,
                    "DOI": normalize_doi(doi),
                    "publication_status": "preprint" if any("Preprint" in x for x in pub_types) else "peer_reviewed_or_indexed",
                    "abstract": abstract,
                    "source_databases": "PubMed",
                    "source_queries": name,
                }
            )
        time.sleep(0.36)
    return records, {"query_name": name, "query": query, "count": count, "retrieved": len(records)}


def crossref_search(query: str, email: str, user_agent: str) -> tuple[list[dict], dict]:
    params = {"query.bibliographic": query, "rows": "100", "select": "DOI,title,author,published,container-title,type"}
    if email:
        params["mailto"] = email
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    raw = request_text(url, user_agent)
    payload = json.loads(raw)
    records = []
    for item in payload.get("message", {}).get("items", []):
        title = (item.get("title") or [""])[0]
        authors = []
        for author in item.get("author", []):
            authors.append(" ".join(x for x in (author.get("family", ""), author.get("given", "")) if x))
        date_parts = ((item.get("published") or {}).get("date-parts") or [[""]])[0]
        year = str(date_parts[0]) if date_parts else ""
        record_type = item.get("subtype") or item.get("type") or ""
        records.append(
            {
                "title": title,
                "first_author": authors[0] if authors else "",
                "authors": "; ".join(authors),
                "publication_year": year,
                "journal": (item.get("container-title") or [""])[0],
                "PMID": "",
                "DOI": normalize_doi(item.get("DOI", "")),
                "publication_status": "preprint" if "preprint" in record_type.lower() else "peer_reviewed_or_indexed",
                "abstract": "",
                "source_databases": "CrossRef",
                "source_queries": query,
            }
        )
    return records, {"query": query, "retrieved": len(records)}


def merge_records(records: list[dict]) -> list[dict]:
    merged: list[dict] = []
    for record in records:
        match = None
        doi = normalize_doi(record.get("DOI", ""))
        for existing in merged:
            existing_doi = normalize_doi(existing.get("DOI", ""))
            if doi and existing_doi and doi == existing_doi:
                match = existing
                break
            if (
                first_author_surname(record.get("first_author", ""))
                == first_author_surname(existing.get("first_author", ""))
                and title_similarity(record.get("title", ""), existing.get("title", "")) >= 0.90
            ):
                match = existing
                break
        if match is None:
            merged.append(record.copy())
            continue
        for field in ("title", "first_author", "authors", "publication_year", "journal", "PMID", "DOI", "abstract"):
            if not match.get(field) and record.get(field):
                match[field] = record[field]
        match["source_databases"] = ";".join(sorted(set(match["source_databases"].split(";") + record["source_databases"].split(";"))))
        match["source_queries"] = ";".join(sorted(set(match["source_queries"].split(";") + record["source_queries"].split(";"))))
        if match.get("publication_status") == "preprint" and record.get("publication_status") != "preprint":
            match["publication_status"] = record["publication_status"]
    return merged


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = [
        "record_id",
        "title",
        "first_author",
        "authors",
        "publication_year",
        "journal",
        "PMID",
        "DOI",
        "publication_status",
        "source_databases",
        "source_queries",
        "abstract",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = row.copy()
            if row.get("PMID"):
                stable_id = f"PMID{row['PMID']}"
            elif normalize_doi(row.get("DOI", "")):
                digest = hashlib.sha1(normalize_doi(row["DOI"]).encode("utf-8")).hexdigest()[:12]
                stable_id = f"DOI_{digest}"
            else:
                key = f"{row.get('title', '')}|{row.get('first_author', '')}".lower()
                stable_id = f"REC_{hashlib.sha1(key.encode('utf-8')).hexdigest()[:12]}"
            out["record_id"] = stable_id
            writer.writerow(out)


def main() -> int:
    config = yaml.safe_load(QUERY_FILE.read_text(encoding="utf-8"))
    email = config.get("contact_email", "")
    user_agent = f"{config.get('ncbi_tool', 'A1_protein_upgrade')}/1.0"
    if email:
        user_agent += f" (mailto:{email})"
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    all_records: list[dict] = []
    run_log = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "declared_search_date": str(config.get("search_date", "")),
        "pubmed": [],
        "crossref": [],
        "errors": [],
    }

    for name, query in config.get("pubmed", {}).items():
        raw_file = RAW_DIR / f"pubmed_{name}.json"
        if raw_file.exists():
            records = json.loads(raw_file.read_text(encoding="utf-8"))
            all_records.extend(records)
            run_log["pubmed"].append(
                {
                    "query_name": name,
                    "query": query,
                    "count": len(records),
                    "retrieved": len(records),
                    "cached": True,
                }
            )
            print(f"PubMed {name}: {len(records)} cached records", flush=True)
            LOG_FILE.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")
            continue
        try:
            print(f"PubMed {name}: querying", flush=True)
            records, summary = pubmed_search(name, query, config, user_agent)
            all_records.extend(records)
            run_log["pubmed"].append(summary)
            raw_file.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
            print(f"PubMed {name}: {len(records)} records", flush=True)
        except Exception as exc:
            run_log["errors"].append({"source": "PubMed", "query": name, "error": repr(exc)})
            print(f"PubMed {name}: failed: {exc!r}", file=sys.stderr, flush=True)
        LOG_FILE.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")

    for index, query in enumerate(config.get("crossref", []), 1):
        raw_file = RAW_DIR / f"crossref_{index:02d}.json"
        if raw_file.exists():
            records = json.loads(raw_file.read_text(encoding="utf-8"))
            all_records.extend(records)
            run_log["crossref"].append({"query": query, "retrieved": len(records), "cached": True})
            print(f"CrossRef {index}: {len(records)} cached records", flush=True)
            LOG_FILE.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")
            continue
        try:
            print(f"CrossRef {index}: querying", flush=True)
            records, summary = crossref_search(query, email, user_agent)
            all_records.extend(records)
            run_log["crossref"].append(summary)
            raw_file.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")
            print(f"CrossRef {index}: {len(records)} records", flush=True)
            time.sleep(0.2)
        except Exception as exc:
            run_log["errors"].append({"source": "CrossRef", "query": query, "error": repr(exc)})
            print(f"CrossRef {index}: failed: {exc!r}", file=sys.stderr, flush=True)
        LOG_FILE.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")

    deduplicated = merge_records(all_records)
    deduplicated.sort(key=lambda x: (str(x.get("publication_year", "")), x.get("title", "")), reverse=True)
    write_tsv(OUT_FILE, deduplicated)
    run_log["raw_record_count"] = len(all_records)
    run_log["deduplicated_record_count"] = len(deduplicated)
    run_log["finished_utc"] = datetime.now(timezone.utc).isoformat()
    LOG_FILE.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Raw records: {len(all_records)}")
    print(f"Deduplicated records: {len(deduplicated)}")
    print(f"Output: {OUT_FILE}")
    if run_log["errors"]:
        print(f"Source errors: {len(run_log['errors'])}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
