#!/usr/bin/env python3
"""Inventory and selectively download UKB-PPP files from Synapse.

Authentication is read only from SYNAPSE_AUTH_TOKEN. The token is never written
to manifests or logs. Downloads use bounded HTTP Range requests, retain verified
partial bytes across restarts, and validate completed files against Synapse MD5.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
import time
from pathlib import Path

import requests
import synapseclient


DEFAULT_ENTITY = "syn51365303"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--entity", default=DEFAULT_ENTITY)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--select", type=Path, help="TSV with a required synapse_id column")
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--reuse-inventory", action="store_true")
    parser.add_argument("--retries", type=int, default=12)
    parser.add_argument("--retry-wait", type=int, default=30)
    return parser.parse_args()


def connect(retries: int = 12, retry_wait: int = 15) -> synapseclient.Synapse:
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise SystemExit("SYNAPSE_AUTH_TOKEN is required; do not place it in a file.")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            syn = synapseclient.Synapse()
            syn.multi_threaded = False
            syn.login(authToken=token, silent=True)
            return syn
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(retry_wait * min(attempt, 6))
    raise RuntimeError(f"Synapse login failed after {retries} attempts: {last_error}")


def inventory_tree(syn: synapseclient.Synapse, root_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    queue = [(root_id, "")]
    while queue:
        parent_id, relative_parent = queue.pop(0)
        for child in syn.getChildren(parent_id):
            child_id = child["id"]
            child_name = child.get("name", child_id)
            child_type = child.get("type", "unknown")
            relative_path = str(Path(relative_parent, child_name))
            row = {
                "synapse_id": child_id,
                "parent_id": parent_id,
                "entity_type": child_type,
                "name": child_name,
                "relative_path": relative_path,
                "content_size": "",
                "content_md5": "",
                "version_number": "",
            }
            if child_type == "org.sagebionetworks.repo.model.Folder":
                queue.append((child_id, relative_path))
            rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "synapse_id", "parent_id", "entity_type", "name", "relative_path",
        "content_size", "content_md5", "version_number", "data_file_handle_id",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_inventory(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row.setdefault("data_file_handle_id", "")
    return rows


def retry_call(function, retries: int, retry_wait: int, label: str):
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            return function()
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(retry_wait * min(attempt, 6))
    raise RuntimeError(f"{label} failed after {retries} attempts: {last_error}")


def read_selection(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "synapse_id" not in reader.fieldnames:
            raise SystemExit(f"Selection file lacks synapse_id: {path}")
        return [
            row["synapse_id"].strip()
            for row in reader
            if row["synapse_id"].strip()
            and ("selection_status" not in row or row["selection_status"] == "selected_for_download")
        ]


def resolve_selected_metadata(
    syn: synapseclient.Synapse,
    selected: list[str],
    inventory: dict[str, dict[str, str]],
) -> int:
    total_bytes = 0
    for synapse_id in selected:
        meta = inventory.get(synapse_id)
        if not meta:
            continue
        entity = syn.restGET(f"/entity/{synapse_id}")
        handles = syn.restGET(f"/entity/{synapse_id}/filehandles").get("list", [])
        handle_id = str(entity.get("dataFileHandleId", ""))
        handle = next((item for item in handles if str(item.get("id", "")) == handle_id), {})
        meta["content_size"] = str(handle.get("contentSize", ""))
        meta["content_md5"] = str(handle.get("contentMd5", ""))
        meta["version_number"] = str(entity.get("versionNumber", ""))
        meta["data_file_handle_id"] = handle_id
        total_bytes += int(meta["content_size"] or 0)
    return total_bytes


def md5(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_in_ranges(
    syn: synapseclient.Synapse,
    synapse_id: str,
    handle_id: str,
    partial_path: Path,
    total_size: int,
    retries: int,
    retry_wait: int,
    range_size: int = 8 * 1024 * 1024,
) -> int:
    """Resume a file in bounded HTTP ranges, rolling back only a failed range."""
    current_size = partial_path.stat().st_size if partial_path.exists() else 0
    if current_size > total_size:
        raise RuntimeError(f"Partial file is larger than expected: {current_size} > {total_size}")

    attempts_used = 0
    while current_size < total_size:
        range_start = current_size
        range_end = min(range_start + range_size - 1, total_size - 1)
        last_error = ""
        for attempt in range(1, retries + 1):
            attempts_used += 1
            try:
                signed = syn._getFileHandleDownload(
                    fileHandleId=handle_id,
                    objectId=synapse_id,
                    objectType="FileEntity",
                )
                url = signed.get("preSignedURL", "")
                if not url:
                    raise RuntimeError("Synapse did not return a pre-signed URL")
                with requests.get(
                    url,
                    headers={"Range": f"bytes={range_start}-{range_end}"},
                    stream=True,
                    timeout=(30, 60),
                ) as response:
                    if response.status_code != 206:
                        raise RuntimeError(f"Expected HTTP 206, received {response.status_code}")
                    expected_range = f"bytes {range_start}-{range_end}/{total_size}"
                    if response.headers.get("Content-Range") != expected_range:
                        raise RuntimeError(
                            f"Unexpected Content-Range: {response.headers.get('Content-Range')}"
                        )
                    with partial_path.open("ab") as handle:
                        for chunk in response.iter_content(chunk_size=4 * 1024 * 1024):
                            if chunk:
                                handle.write(chunk)
                        handle.flush()
                        os.fsync(handle.fileno())
                observed_size = partial_path.stat().st_size
                if observed_size != range_end + 1:
                    raise RuntimeError(
                        f"Range length mismatch: observed {observed_size}, expected {range_end + 1}"
                    )
                current_size = observed_size
                print(
                    f"{synapse_id}: {current_size}/{total_size} bytes "
                    f"({100 * current_size / total_size:.1f}%)",
                    flush=True,
                )
                break
            except Exception as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                if partial_path.exists():
                    with partial_path.open("r+b") as handle:
                        handle.truncate(range_start)
                if attempt < retries:
                    time.sleep(retry_wait * min(attempt, 6))
        else:
            raise RuntimeError(
                f"Range {range_start}-{range_end} failed after {retries} attempts: {last_error}"
            )
    return attempts_used


def download_selected(
    syn: synapseclient.Synapse,
    selected: list[str],
    inventory: dict[str, dict[str, str]],
    output: Path,
    retries: int,
    retry_wait: int,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    state_path = output / "download_state.tsv"
    state_fields = ["synapse_id", "path", "status", "bytes", "md5", "attempts", "message"]
    states: list[dict[str, str]] = []

    def persist_state() -> None:
        with state_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=state_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(states)

    for synapse_id in selected:
        meta = inventory.get(synapse_id)
        if not meta:
            states.append({
                "synapse_id": synapse_id, "path": "", "status": "not_in_inventory",
                "bytes": "", "md5": "", "attempts": "0", "message": "selection ID not found",
            })
            persist_state()
            continue
        entity_meta = retry_call(
            lambda: syn.restGET(f"/entity/{synapse_id}"), retries, retry_wait,
            f"entity metadata {synapse_id}",
        )
        handle_id = str(entity_meta.get("dataFileHandleId", ""))
        if not meta.get("content_size") or not meta.get("content_md5"):
            handles = retry_call(
                lambda: syn.restGET(f"/entity/{synapse_id}/filehandles"), retries, retry_wait,
                f"file-handle metadata {synapse_id}",
            ).get("list", [])
            handle = next((item for item in handles if str(item.get("id", "")) == handle_id), {})
            meta["content_size"] = str(handle.get("contentSize", ""))
            meta["content_md5"] = str(handle.get("contentMd5", ""))
        meta["version_number"] = str(entity_meta.get("versionNumber", ""))
        meta["data_file_handle_id"] = handle_id
        relative_parent = Path(meta["relative_path"]).parent
        destination = output / relative_parent
        destination.mkdir(parents=True, exist_ok=True)
        expected_path = destination / meta["name"]
        expected_md5 = meta["content_md5"].lower()

        if expected_path.exists() and expected_md5 and md5(expected_path) == expected_md5:
            states.append({
                "synapse_id": synapse_id, "path": str(expected_path), "status": "verified_existing",
                "bytes": str(expected_path.stat().st_size), "md5": expected_md5,
                "attempts": "0", "message": "",
            })
            persist_state()
            continue

        partial_path = expected_path.with_name(expected_path.name + ".part")
        last_error = ""
        try:
            attempts_used = download_in_ranges(
                syn=syn,
                synapse_id=synapse_id,
                handle_id=handle_id,
                partial_path=partial_path,
                total_size=int(meta["content_size"]),
                retries=retries,
                retry_wait=retry_wait,
            )
            observed_md5 = md5(partial_path)
            if expected_md5 and observed_md5 != expected_md5:
                raise RuntimeError(f"MD5 mismatch: {observed_md5} != {expected_md5}")
            os.replace(partial_path, expected_path)
            states.append({
                "synapse_id": synapse_id, "path": str(expected_path), "status": "downloaded_verified",
                "bytes": str(expected_path.stat().st_size), "md5": observed_md5,
                "attempts": str(attempts_used), "message": "",
            })
            print(f"[{len(states)}/{len(selected)}] verified {meta['name']}")
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            states.append({
                "synapse_id": synapse_id, "path": str(expected_path), "status": "failed_resumable",
                "bytes": str(partial_path.stat().st_size) if partial_path.exists() else "0",
                "md5": "", "attempts": str(retries), "message": last_error,
            })

        persist_state()

    persist_state()
    failures = [row for row in states if row["status"] not in {"verified_existing", "downloaded_verified"}]
    if failures:
        failed_ids = ", ".join(row["synapse_id"] for row in failures)
        raise RuntimeError(f"Selected downloads failed integrity or availability checks: {failed_ids}")


def main() -> int:
    args = parse_args()
    syn = connect(retries=args.retries, retry_wait=args.retry_wait)
    if args.reuse_inventory and args.inventory.exists():
        rows = read_inventory(args.inventory)
        print(f"Reused inventory: {args.inventory}")
    else:
        rows = inventory_tree(syn, args.entity)
        write_tsv(args.inventory, rows)
    file_rows = {
        row["synapse_id"]: row
        for row in rows
        if row["entity_type"] == "org.sagebionetworks.repo.model.FileEntity"
    }
    print(f"Inventoried {len(file_rows)} files; size and MD5 are resolved only for selected files.")
    selected = read_selection(args.select) if args.select else []
    if selected and not args.reuse_inventory:
        total_bytes = resolve_selected_metadata(syn, selected, file_rows)
        write_tsv(args.inventory, rows)
        print(f"Selected {len(selected)} files ({total_bytes / 1024**3:.2f} GiB).")
    if args.inventory_only:
        return 0
    if not args.select:
        raise SystemExit("--select is required unless --inventory-only is used")
    download_selected(
        syn, selected, file_rows, args.output,
        retries=args.retries, retry_wait=args.retry_wait,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
