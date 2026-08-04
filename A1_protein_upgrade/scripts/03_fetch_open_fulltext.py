"""Fetch open PubMed Central full text and supplementary files for review records."""

from __future__ import annotations

import csv
import io
import json
import re
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SEARCH_RESULTS = ROOT / "literature" / "search_results_deduplicated.tsv"
QUEUE_CONFIG = ROOT / "config" / "full_text_review_records.yml"
OUT_DIR = ROOT / "data_raw" / "literature_fulltext"
AVAILABILITY = ROOT / "literature" / "full_text_availability.tsv"
ID_CACHE = ROOT / "data_raw" / "literature_fulltext" / "pmid_pmcid_conversion.json"
USER_AGENT = "A1_protein_upgrade/1.0"
XLINK = "{http://www.w3.org/1999/xlink}href"


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def get_bytes(url: str, timeout: int = 60, attempts: int = 3) -> bytes:
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except Exception as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(1.5 * (attempt + 1))
    assert last_error is not None
    raise last_error


def id_convert(pmids: list[str]) -> dict[str, dict]:
    output: dict[str, dict] = {}
    if ID_CACHE.exists():
        output.update(json.loads(ID_CACHE.read_text(encoding="utf-8")))
    missing = [pmid for pmid in pmids if pmid not in output]
    for start in range(0, len(missing), 50):
        chunk = missing[start : start + 50]
        params = {"format": "json", "ids": ",".join(chunk)}
        url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?" + urllib.parse.urlencode(params)
        payload = json.loads(get_bytes(url).decode("utf-8"))
        for record in payload.get("records", []):
            if record.get("pmid"):
                output[str(record["pmid"])] = record
        ID_CACHE.parent.mkdir(parents=True, exist_ok=True)
        ID_CACHE.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    return output


def clean_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())


def supplementary_links(root: ET.Element) -> list[str]:
    links = []
    for element in root.findall(".//supplementary-material"):
        for node in element.iter():
            href = node.attrib.get(XLINK) or node.attrib.get("href")
            if href:
                links.append(href)
    return sorted(set(links))


def fetch_oa_supplements(pmcid: str, hrefs: list[str], target_dir: Path) -> tuple[list[str], str]:
    supplement_dir = target_dir / "supplements"
    supplement_dir.mkdir(parents=True, exist_ok=True)
    cached = [path for path in supplement_dir.iterdir() if path.is_file() and path.stat().st_size > 0]
    if cached:
        return [str(path.relative_to(ROOT)).replace("\\", "/") for path in cached], "cached_OA_package"

    try:
        europe_pmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles"
        payload = get_bytes(europe_pmc_url, timeout=120)
        if payload.startswith(b"PK"):
            zip_file = target_dir / f"{pmcid}_europepmc_supplementary.zip"
            zip_file.write_bytes(payload)
            extracted = []
            allowed_extensions = {".csv", ".tsv", ".txt", ".xls", ".xlsx", ".doc", ".docx", ".pdf", ".zip"}
            with zipfile.ZipFile(io.BytesIO(payload)) as archive:
                for name in archive.namelist():
                    basename = Path(name).name
                    if not basename or Path(basename.lower()).suffix not in allowed_extensions:
                        continue
                    target = supplement_dir / basename
                    target.write_bytes(archive.read(name))
                    extracted.append(target)
            if extracted:
                return [str(path.relative_to(ROOT)).replace("\\", "/") for path in extracted], "downloaded_Europe_PMC_supplementary"
            return [], "Europe_PMC_zip_no_supported_files"
    except Exception:
        pass

    try:
        oa_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?" + urllib.parse.urlencode({"id": pmcid})
        oa_root = ET.fromstring(get_bytes(oa_url, timeout=30))
        package_link = ""
        for link in oa_root.findall(".//link"):
            if link.attrib.get("format") == "tgz":
                package_link = link.attrib.get("href", "")
                break
        if not package_link:
            return [], "OA_package_not_available"
        if package_link.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
            package_link = package_link.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/", 1)

        payload = get_bytes(package_link, timeout=120)
        if len(payload) > 250 * 1024 * 1024:
            return [], "OA_package_skipped_over_250MB"
        package_file = target_dir / f"{pmcid}_oa_package.tar.gz"
        package_file.write_bytes(payload)

        href_names = {Path(urllib.parse.urlparse(href).path).name.lower() for href in hrefs}
        allowed_extensions = {".csv", ".tsv", ".txt", ".xls", ".xlsx", ".doc", ".docx", ".pdf", ".zip"}
        extracted = []
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
            for member in archive.getmembers():
                if not member.isfile():
                    continue
                basename = Path(member.name).name
                lower = basename.lower()
                extension = Path(lower).suffix
                looks_supplementary = (
                    lower in href_names
                    or any(token in lower for token in ("supp", "moesm", "esm", "mmc", "s001", "s002", "s003"))
                )
                if extension not in allowed_extensions or not looks_supplementary:
                    continue
                source = archive.extractfile(member)
                if source is None:
                    continue
                target = supplement_dir / basename
                target.write_bytes(source.read())
                extracted.append(target)
        status = "downloaded_OA_package" if extracted else "OA_package_downloaded_no_matching_supplement"
        return [str(path.relative_to(ROOT)).replace("\\", "/") for path in extracted], status
    except Exception as exc:
        return [], f"OA_package_failed:{type(exc).__name__}"


def main() -> int:
    records = {row["record_id"]: row for row in read_tsv(SEARCH_RESULTS)}
    queue = yaml.safe_load(QUEUE_CONFIG.read_text(encoding="utf-8"))["records"]
    selected = []
    for item in queue:
        record = records.get(item["record_id"])
        if record is None:
            selected.append({"record_id": item["record_id"], "priority_reason": item["priority_reason"], "missing": True})
        else:
            selected.append({**record, "priority_reason": item["priority_reason"], "missing": False})

    pmids = [item.get("PMID", "") for item in selected if item.get("PMID")]
    conversions = id_convert(pmids)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []

    def write_availability(current_rows: list[dict]) -> None:
        if not current_rows:
            return
        fields = list(current_rows[0].keys())
        with AVAILABILITY.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(current_rows)

    for item in selected:
        pmid = item.get("PMID", "")
        conversion = conversions.get(pmid, {})
        pmcid = conversion.get("pmcid", "")
        full_text_status = "not_available_in_PMC"
        xml_path = ""
        text_path = ""
        supplement_status = "not_assessed"
        supplement_files: list[str] = []
        errors: list[str] = []

        if item.get("missing"):
            full_text_status = "record_missing_after_dedup"
        elif pmcid:
            article_dir = OUT_DIR / f"{item['record_id']}_{pmcid}"
            article_dir.mkdir(parents=True, exist_ok=True)
            try:
                xml_file = article_dir / f"{pmcid}.xml"
                if xml_file.exists() and xml_file.stat().st_size > 0:
                    xml_bytes = xml_file.read_bytes()
                    full_text_status = "PMC_full_text_cached"
                else:
                    params = {"db": "pmc", "id": pmcid, "retmode": "xml"}
                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
                    xml_bytes = get_bytes(url, timeout=45)
                    xml_file.write_bytes(xml_bytes)
                    full_text_status = "PMC_full_text_downloaded"
                xml_path = str(xml_file.relative_to(ROOT)).replace("\\", "/")
                root = ET.fromstring(xml_bytes)

                body_text = "\n\n".join(
                    text for text in (clean_text(node) for node in root.findall(".//body//p")) if text
                )
                table_text = "\n\n".join(
                    clean_text(node) for node in root.findall(".//table-wrap") if clean_text(node)
                )
                text_file = article_dir / f"{pmcid}_fulltext.txt"
                if not text_file.exists() or text_file.stat().st_size == 0:
                    text_file.write_text(body_text + "\n\nTABLES\n\n" + table_text, encoding="utf-8")
                text_path = str(text_file.relative_to(ROOT)).replace("\\", "/")

                links = supplementary_links(root)
                if links:
                    supplement_files, supplement_status = fetch_oa_supplements(pmcid, links, article_dir)
                else:
                    supplement_status = "no_supplement_link_in_JATS"
            except Exception as exc:
                full_text_status = "PMC_fetch_failed"
                errors.append(repr(exc))

        rows.append(
            {
                "record_id": item.get("record_id", ""),
                "priority_reason": item.get("priority_reason", ""),
                "title": item.get("title", ""),
                "PMID": pmid,
                "DOI": item.get("DOI", ""),
                "PMCID": pmcid,
                "full_text_status": full_text_status,
                "xml_path": xml_path,
                "text_path": text_path,
                "supplement_status": supplement_status,
                "supplement_files": ";".join(supplement_files),
                "manual_verification_status": "pending",
                "errors": ";".join(errors),
            }
        )
        write_availability(rows)
        print(f"{item.get('record_id')}: {full_text_status}; {supplement_status}", flush=True)

    write_availability(rows)
    print(f"Output: {AVAILABILITY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
