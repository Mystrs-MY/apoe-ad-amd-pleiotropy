#!/usr/bin/env python3
import argparse
import csv
import hashlib
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import botocore.session
from botocore.config import Config


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "tables" / "decode_target_object_inventory.tsv"
DEFAULT_OUTPUT = ROOT / "data_raw" / "decode_somascan"
DEFAULT_LEDGER = ROOT / "logs" / "decode_download_ledger.tsv"
PROFILE = "decode-download"
ENDPOINT = "https://s3-ext.decode.is"
CHUNK_SIZE = 8 * 1024 * 1024
REPORT_INTERVAL = 256 * 1024 * 1024
RESERVE_BYTES = 2 * 1024 ** 3


def client_from_profile():
    session = botocore.session.Session(profile=PROFILE)
    return session.create_client(
        "s3",
        endpoint_url=ENDPOINT,
        region_name="us-east-1",
        config=Config(
            signature_version="s3v4",
            retries={"max_attempts": 10, "mode": "standard"},
            connect_timeout=30,
            read_timeout=180,
            max_pool_connections=2,
            s3={"addressing_style": "path"},
        ),
    )


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_ledger(path):
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {(row["normalization"], row["SomaScan_SeqId"]): row for row in rows}


def write_ledger(path, records):
    fields = [
        "normalization",
        "gene_symbol",
        "SomaScan_SeqId",
        "object_key",
        "local_path",
        "content_length_bytes",
        "etag",
        "sha256",
        "completed_utc",
        "status",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    ordered = [records[key] for key in sorted(records)]
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ordered)
    temporary.replace(path)


def verify_remote(client, row):
    head = client.head_object(Bucket=row["bucket"], Key=row["object_key"])
    size = int(row["content_length_bytes"])
    etag = row["etag"]
    if head["ContentLength"] != size:
        raise RuntimeError(f"Remote size changed for {row['object_key']}")
    if head["ETag"].strip('"') != etag:
        raise RuntimeError(f"Remote ETag changed for {row['object_key']}")
    return size


def download_with_resume(client, row, output_root, max_retries):
    expected = verify_remote(client, row)
    destination = output_root / row["normalization"] / Path(row["object_key"]).name
    partial = destination.with_suffix(destination.suffix + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists():
        if destination.stat().st_size != expected:
            raise RuntimeError(f"Completed file has unexpected size: {destination}")
        return destination, "already_complete"

    if partial.exists() and partial.stat().st_size > expected:
        raise RuntimeError(f"Partial file exceeds remote size: {partial}")

    free = shutil.disk_usage(destination.parent).free
    remaining = expected - (partial.stat().st_size if partial.exists() else 0)
    if free < remaining + RESERVE_BYTES:
        raise RuntimeError(
            f"Insufficient free space for {destination.name}: need {remaining + RESERVE_BYTES}, have {free}"
        )

    retries = 0
    while True:
        offset = partial.stat().st_size if partial.exists() else 0
        if offset == expected:
            os.replace(partial, destination)
            return destination, "resumed_complete"
        if retries > max_retries:
            raise RuntimeError(f"Retry limit exceeded for {row['object_key']}")

        body = None
        try:
            response = client.get_object(
                Bucket=row["bucket"],
                Key=row["object_key"],
                Range=f"bytes={offset}-",
                IfMatch=row["etag"],
            )
            content_range = response.get("ContentRange", "")
            if not content_range.startswith(f"bytes {offset}-"):
                raise RuntimeError(
                    f"Unexpected Content-Range for {row['object_key']}: {content_range}"
                )
            body = response["Body"]
            checkpoint = offset
            started = time.monotonic()
            with partial.open("ab") as handle:
                while True:
                    chunk = body.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    handle.write(chunk)
                    offset += len(chunk)
                    if offset - checkpoint >= REPORT_INTERVAL:
                        handle.flush()
                        os.fsync(handle.fileno())
                        elapsed = max(time.monotonic() - started, 0.001)
                        rate = (offset - checkpoint) / elapsed / (1024 ** 2)
                        print(
                            f"{row['normalization']} {row['gene_symbol']} {row['SomaScan_SeqId']}: "
                            f"{offset / expected:.1%} ({rate:.1f} MiB/s)",
                            flush=True,
                        )
                        checkpoint = offset
                        started = time.monotonic()
                handle.flush()
                os.fsync(handle.fileno())
            body.close()
            if partial.stat().st_size != expected:
                raise IOError(
                    f"Stream ended at {partial.stat().st_size}, expected {expected}"
                )
            os.replace(partial, destination)
            return destination, "downloaded_complete"
        except Exception as error:
            if body is not None:
                body.close()
            retries += 1
            delay = min(60, 2 ** min(retries, 6))
            current = partial.stat().st_size if partial.exists() else 0
            print(
                f"Retry {retries}/{max_retries} for {row['SomaScan_SeqId']} from byte "
                f"{current} after {type(error).__name__}: {error}",
                flush=True,
            )
            time.sleep(delay)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--normalization", choices=["smp", "raw", "both"], default="smp")
    parser.add_argument("--assays", nargs="+", default=None,
                        help="Optional exact SomaScan SeqIds for a gated subset download")
    parser.add_argument("--max-retries", type=int, default=50)
    return parser.parse_args()


def main():
    args = parse_args()
    with args.manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = rows if args.normalization == "both" else [
        row for row in rows if row["normalization"] == args.normalization
    ]
    if args.assays:
        requested = set(args.assays)
        selected = [row for row in selected if row["SomaScan_SeqId"] in requested]
        observed = {row["SomaScan_SeqId"] for row in selected}
        if observed != requested:
            raise RuntimeError(f"Requested assays absent from manifest: {sorted(requested - observed)}")
    expected_count = (18 if args.normalization == "both" else 9) if not args.assays else (
        len(set(args.assays)) * (2 if args.normalization == "both" else 1)
    )
    if len(selected) != expected_count:
        raise RuntimeError(f"Expected {expected_count} manifest rows, found {len(selected)}")

    client = client_from_profile()
    ledger = load_ledger(args.ledger)
    for row in selected:
        print(f"Starting {row['normalization']} {row['gene_symbol']} {row['SomaScan_SeqId']}")
        path, transfer_status = download_with_resume(
            client, row, args.output_root, args.max_retries
        )
        digest = sha256_file(path)
        key = (row["normalization"], row["SomaScan_SeqId"])
        ledger[key] = {
            "normalization": row["normalization"],
            "gene_symbol": row["gene_symbol"],
            "SomaScan_SeqId": row["SomaScan_SeqId"],
            "object_key": row["object_key"],
            "local_path": str(path.resolve()),
            "content_length_bytes": row["content_length_bytes"],
            "etag": row["etag"],
            "sha256": digest,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "status": transfer_status,
        }
        write_ledger(args.ledger, ledger)
        print(f"Verified SHA-256 {digest}  {path.name}")

    print(f"Completed {len(selected)} target objects; ledger: {args.ledger}")


if __name__ == "__main__":
    main()
