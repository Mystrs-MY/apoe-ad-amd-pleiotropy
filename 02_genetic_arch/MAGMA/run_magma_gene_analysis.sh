#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RESOURCE="${A1_RESOURCE_ROOT_WSL:-/mnt/d/AD_AMD/Resource}"
MAGMA="$ROOT/A1_protein_upgrade/tools/magma_v1.10/extracted/magma"
GENE_LOC="$ROOT/A1_protein_upgrade/tools/magma_v1.10/gene_location/Rev.NCBI37.3.gene.loc"
LD_PREFIX="$RESOURCE/EUR/EUR"
OUT="$ROOT/02_genetic_arch/MAGMA"
INPUT_DIR="$OUT/input"
RESULT_DIR="$OUT/results"
LOG_DIR="$OUT/logs"
ANNOT_PREFIX="$OUT/NCBI37_3_EUR_window35_10"

mkdir -p "$INPUT_DIR" "$RESULT_DIR" "$LOG_DIR"

for file in "$MAGMA" "$GENE_LOC" "$LD_PREFIX.bed" "$LD_PREFIX.bim" "$LD_PREFIX.fam"; do
  if [[ ! -s "$file" ]]; then
    echo "Missing required input: $file" >&2
    exit 1
  fi
done

chmod +x "$MAGMA"
"$MAGMA" --version | tee "$LOG_DIR/magma_version.txt"

if [[ ! -s "$ANNOT_PREFIX.genes.annot" ]]; then
  "$MAGMA" \
    --annotate window=35,10 \
    --snp-loc "$LD_PREFIX.bim" \
    --gene-loc "$GENE_LOC" \
    --out "$ANNOT_PREFIX" \
    > "$LOG_DIR/annotation.stdout.log" 2> "$LOG_DIR/annotation.stderr.log"
fi

declare -A GWAS_FILES=(
  [AD]="AD_Wightman_cleaned_hg19.tsv.gz"
  [Dry_AMD]="AMD_Dry_R12_cleaned_hg19.tsv.gz"
  [Wet_AMD]="AMD_Wet_R12_cleaned_hg19.tsv.gz"
  [Any_AMD]="AMD_H7_R12_cleaned_hg19.tsv.gz"
)

for trait in AD Dry_AMD Wet_AMD Any_AMD; do
  gwas="$RESOURCE/GWAS/${GWAS_FILES[$trait]}"
  pval="$INPUT_DIR/${trait}.pval.tsv"
  prefix="$RESULT_DIR/$trait"

  if [[ ! -s "$gwas" ]]; then
    echo "Missing GWAS: $gwas" >&2
    exit 1
  fi

  if [[ ! -s "$pval" ]]; then
    tmp="$pval.tmp"
    gzip -dc "$gwas" | awk -F '\t' 'BEGIN {OFS="\t"}
      NR==1 {
        sub(/\r$/, "", $10);
        if ($1!="SNP" || $9!="P" || $10!="N") {
          print "Unexpected GWAS columns: expected SNP in column 1, P in column 9, N in column 10" > "/dev/stderr";
          exit 2
        }
        print "SNP", "P", "N";
        next
      }
      {
        sub(/\r$/, "", $10);
        if ($1!="" && $9!="" && $10!="") print $1, $9, $10;
      }
    ' > "$tmp"
    mv "$tmp" "$pval"
  fi

  if [[ ! -s "$prefix.genes.out" ]]; then
    "$MAGMA" \
      --bfile "$LD_PREFIX" \
      --gene-annot "$ANNOT_PREFIX.genes.annot" \
      --pval "$pval" use=SNP,P ncol=N duplicate=error \
      --gene-model snp-wise=mean \
      --genes-only \
      --out "$prefix" \
      > "$LOG_DIR/${trait}.stdout.log" 2> "$LOG_DIR/${trait}.stderr.log"
  fi
done

{
  printf 'sha256\tpath\n'
  sha256sum \
    "$MAGMA" "$GENE_LOC" \
    "$LD_PREFIX.bed" "$LD_PREFIX.bim" "$LD_PREFIX.fam" \
    "$RESOURCE/GWAS/AD_Wightman_cleaned_hg19.tsv.gz" \
    "$RESOURCE/GWAS/AMD_Dry_R12_cleaned_hg19.tsv.gz" \
    "$RESOURCE/GWAS/AMD_Wet_R12_cleaned_hg19.tsv.gz" \
    "$RESOURCE/GWAS/AMD_H7_R12_cleaned_hg19.tsv.gz" \
    | awk 'BEGIN {OFS="\t"} {hash=$1; $1=""; sub(/^  /, ""); print hash, $0}'
} > "$LOG_DIR/input_sha256.tsv"

cat > "$LOG_DIR/run_manifest.tsv" <<EOF
field	value
MAGMA_version	1.10_linux
execution_environment	Ubuntu_WSL2
LD_reference	1000_Genomes_EUR
LD_reference_samples	503
gene_location	Rev.NCBI37.3.gene.loc
gene_location_source	STAR_Protocols_MAGMA_analysis_protocol
gene_window	upstream_35kb_downstream_10kb
gene_model	snp-wise_mean
sample_size_mode	per_variant_N_column
duplicate_SNP_policy	error
EOF

echo "MAGMA gene analysis completed."
