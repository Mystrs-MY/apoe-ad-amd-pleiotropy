"""
C6_extract_ukbppp.py
Extract UKB-PPP .tar files (which contain : in filenames, invalid on Windows),
merge all 22+chrX chromosomes, and write one clean merged .txt.gz per protein.

Output: ukbppp_merged/{gene}_merged.txt.gz
Each merged file has: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
"""

import os, sys, gzip, tarfile, glob, io, re
from pathlib import Path

# Config
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESOURCE_ROOT = Path(os.environ.get("A1_RESOURCE_ROOT", PROJECT_ROOT / "data" / "external"))
TAR_DIR = str(RESOURCE_ROOT / "ukbppp_proteins")
OUT_DIR = str(RESOURCE_ROOT / "ukbppp_proteins" / "ukbppp_merged")
os.makedirs(OUT_DIR, exist_ok=True)

# Find all .tar files
tar_files = sorted(glob.glob(os.path.join(TAR_DIR, '*.tar')))
print(f"Found {len(tar_files)} .tar files")

# Protein name map (extract gene symbol from filename)
# e.g., c3_p01024_oid30776_v1_inflammation_ii.tar -> C3

for i, tar_path in enumerate(tar_files):
    basename = os.path.basename(tar_path).replace('.tar', '')
    gene = basename.split('_')[0].upper()
    out_path = os.path.join(OUT_DIR, f'{gene}_merged.txt.gz')

    # Skip if already processed
    if os.path.exists(out_path):
        print(f"[{i+1}/{len(tar_files)}] {gene:12s} SKIP (already merged)")
        continue

    print(f"[{i+1}/{len(tar_files)}] {gene:12s} extracting...", end=' ', flush=True)

    try:
        with tarfile.open(tar_path, 'r') as tar:
            members = [m for m in tar.getmembers() if m.name.endswith('.gz')]

            # Read header from first file
            first_f = tar.extractfile(members[0])
            with gzip.GzipFile(fileobj=first_f) as gz:
                header = gz.readline()  # keep as bytes

            # Merge all chromosome files, skip header on subsequent files
            total_lines = 0
            with gzip.open(out_path, 'wb') as outf:
                outf.write(header)

                for member in sorted(members, key=lambda m: m.name):
                    f = tar.extractfile(member)
                    with gzip.GzipFile(fileobj=f) as gz:
                        # Skip header for all files
                        gz.readline()
                        while True:
                            chunk = gz.read(65536)  # 64KB chunks
                            if not chunk:
                                break
                            outf.write(chunk)
                            total_lines += chunk.count(b'\n')

        # Count SNPs
        snp_count = 0
        with gzip.open(out_path, 'rt') as f:
            next(f)  # skip header
            for _ in f:
                snp_count += 1

        file_size_mb = os.path.getsize(out_path) / (1024*1024)
        print(f"OK ({snp_count:,} SNPs, {file_size_mb:.0f} MB)")

    except Exception as e:
        print(f"FAIL: {e}")

# Summary
merged_files = sorted(glob.glob(os.path.join(OUT_DIR, '*_merged.txt.gz')))
total_gb = sum(os.path.getsize(f) for f in merged_files) / (1024**3)
print(f"\n{'='*60}")
print(f"DONE: {len(merged_files)} proteins merged ({total_gb:.1f} GB)")
print(f"Output: {OUT_DIR}/")
