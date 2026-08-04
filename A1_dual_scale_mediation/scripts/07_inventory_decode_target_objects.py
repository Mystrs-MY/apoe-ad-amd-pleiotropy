#!/usr/bin/env python3
import csv
from datetime import timezone
from pathlib import Path

import botocore.session
from botocore.config import Config


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "tables" / "decode_same_platform_feasibility_gate.tsv"
OUTPUT = ROOT / "tables" / "decode_target_object_inventory.tsv"
PROFILE = "decode-download"
ENDPOINT = "https://s3-ext.decode.is"
BUCKET = "largescaleplasma-2023"

MODES = {
    "smp": {
        "folder": "final_somascan_smp",
        "prefix": "Proteomics_SMP_PC0",
        "analysis_role": "primary_manufacturer_recommended_normalization",
    },
    "raw": {
        "folder": "final_somascan_raw",
        "prefix": "Proteomics_PC0",
        "analysis_role": "non_normalized_sensitivity",
    },
}


def s3_client():
    session = botocore.session.Session(profile=PROFILE)
    return session.create_client(
        "s3",
        endpoint_url=ENDPOINT,
        region_name="us-east-1",
        config=Config(
            signature_version="s3v4",
            retries={"max_attempts": 10, "mode": "standard"},
            connect_timeout=30,
            read_timeout=120,
            s3={"addressing_style": "path"},
        ),
    )


def main():
    with TARGETS.open("r", encoding="utf-8-sig", newline="") as handle:
        targets = list(csv.DictReader(handle, delimiter="\t"))
    if len(targets) != 9:
        raise RuntimeError(f"Expected 9 exact SomaScan assays, found {len(targets)}")

    client = s3_client()
    rows = []
    for mode_name, mode in MODES.items():
        for target in targets:
            seqid = target["SomaScan_SeqId"]
            prefix = f"{mode['folder']}/{mode['prefix']}_{seqid}_"
            response = client.list_objects_v2(Bucket=BUCKET, Prefix=prefix, MaxKeys=10)
            objects = response.get("Contents", [])
            if len(objects) != 1:
                raise RuntimeError(
                    f"Expected exactly one object for {mode_name}/{seqid}, found {len(objects)}"
                )
            item = objects[0]
            head = client.head_object(Bucket=BUCKET, Key=item["Key"])
            modified = head["LastModified"].astimezone(timezone.utc).isoformat()
            rows.append(
                {
                    "gene_symbol": target["gene_symbol"],
                    "SomaScan_SeqId": seqid,
                    "UniProt_ID": target["UniProt_ID"],
                    "normalization": mode_name,
                    "analysis_role": mode["analysis_role"],
                    "bucket": BUCKET,
                    "object_key": item["Key"],
                    "content_length_bytes": str(head["ContentLength"]),
                    "size_GiB": f"{head['ContentLength'] / (1024 ** 3):.6f}",
                    "etag": head["ETag"].strip('"'),
                    "last_modified_utc": modified,
                    "storage_class": item.get("StorageClass", "STANDARD"),
                    "inventory_status": "exact_remote_object_confirmed",
                    "source_mapping": TARGETS.name,
                }
            )

    fields = list(rows[0].keys())
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(OUTPUT)
    total = sum(int(row["content_length_bytes"]) for row in rows)
    print(f"Inventory written: {OUTPUT}")
    print(f"Objects: {len(rows)}; total: {total / (1024 ** 3):.3f} GiB")


if __name__ == "__main__":
    main()
