#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RESOURCE="${A1_RESOURCE_ROOT_WSL:-$ROOT/data/external}"
MAGMA="${MAGMA_BIN:-magma}"
GENE_LOC="${A1_MAGMA_GENE_LOC:-$ROOT/data/external/Rev.NCBI37.3.gene.loc}"
LD_PREFIX="$RESOURCE/EUR/EUR"
OUT="$ROOT/02_genetic_arch/MAGMA"
RUN_TAG="${A1_MAGMA_RUN_TAG:-correctedN_20260903}"
INPUT_DIR="$OUT/input${RUN_TAG:+_$RUN_TAG}"
RESULT_DIR="$OUT/results${RUN_TAG:+_$RUN_TAG}"
LOG_DIR="$OUT/logs${RUN_TAG:+_$RUN_TAG}"
ANNOT_PREFIX="$OUT/NCBI37_3_EUR_window35_10"

mkdir -p "$INPUT_DIR" "$RESULT_DIR" "$LOG_DIR"

if [[ "$MAGMA" != */* ]]; then
  MAGMA="$(command -v "$MAGMA" || true)"
fi
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

  # Rebuild a cached MAGMA input whenever its corrected GWAS source is newer.
  if [[ ! -s "$pval" || "$gwas" -nt "$pval" ]]; then
    tmp="$pval.tmp"
    gzip -dc "$gwas" | awk -F '\t' 'BEGIN {OFS="\t"}
      NR==1 {
        for (i=1; i<=NF; i++) {
          h=$i; sub(/\r$/, "", h);
          if (h=="SNP") snp_col=i;
          if (h=="P") p_col=i;
          if (h=="N_TOTAL") n_col=i;
        }
        if (!snp_col || !p_col || !n_col) {
          print "Unexpected GWAS columns: SNP, P, and N_TOTAL are required" > "/dev/stderr";
          exit 2
        }
        print "SNP", "P", "N";
        next
      }
      {
        snp=$snp_col; p=$p_col; n=$n_col; sub(/\r$/, "", n);
        if (snp!="" && p!="" && n!="") print snp, p, n;
      }
    ' > "$tmp"
    mv "$tmp" "$pval"
  fi

  # A gene result is reusable only when it is newer than the exact p-value input.
  if [[ ! -s "$prefix.genes.out" || "$pval" -nt "$prefix.genes.out" ]]; then
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
  printf 'resource_id\tsha256\n'
  printf 'MAGMA_binary\t%s\n' "$(sha256sum "$MAGMA" | awk '{print $1}')"
  printf 'NCBI37_3_gene_location\t%s\n' "$(sha256sum "$GENE_LOC" | awk '{print $1}')"
  printf '1000G_EUR_bed\t%s\n' "$(sha256sum "$LD_PREFIX.bed" | awk '{print $1}')"
  printf '1000G_EUR_bim\t%s\n' "$(sha256sum "$LD_PREFIX.bim" | awk '{print $1}')"
  printf '1000G_EUR_fam\t%s\n' "$(sha256sum "$LD_PREFIX.fam" | awk '{print $1}')"
  for trait in AD_Wightman AMD_Dry_R12 AMD_Wet_R12 AMD_H7_R12; do
    case "$trait" in
      AD_Wightman) file="$RESOURCE/GWAS/AD_Wightman_cleaned_hg19.tsv.gz" ;;
      AMD_Dry_R12) file="$RESOURCE/GWAS/AMD_Dry_R12_cleaned_hg19.tsv.gz" ;;
      AMD_Wet_R12) file="$RESOURCE/GWAS/AMD_Wet_R12_cleaned_hg19.tsv.gz" ;;
      AMD_H7_R12) file="$RESOURCE/GWAS/AMD_H7_R12_cleaned_hg19.tsv.gz" ;;
    esac
    printf '%s\t%s\n' "$trait" "$(sha256sum "$file" | awk '{print $1}')"
  done
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
sample_size_mode	per_variant_total_analysed_sample_size
sample_size_source_column	N_TOTAL
duplicate_SNP_policy	error
EOF

echo "MAGMA gene analysis completed."
