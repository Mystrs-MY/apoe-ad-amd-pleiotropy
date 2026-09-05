#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_TAG="${A1_MIXER_RUN_TAG:-correctedN_20260903}"
INPUT_DIR="${SCRIPT_DIR}/inputs_${RUN_TAG}"
OUTPUT_DIR="${SCRIPT_DIR}/results_${RUN_TAG}"
PYTHON_BIN="${MIXER_PYTHON:-python}"
MIXER_PY="${MIXER_SCRIPT:-${SCRIPT_DIR}/gsa-mixer/precimed/mixer.py}"
LIB_BGMG="${MIXER_LIB:-${SCRIPT_DIR}/gsa-mixer/src/build/lib/libbgmg.so}"
BIM="${MIXER_BIM_TEMPLATE:-${SCRIPT_DIR}/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim}"
LD="${MIXER_LD_TEMPLATE:-${SCRIPT_DIR}/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld}"
EXTRACT="${MIXER_EXTRACT:-${SCRIPT_DIR}/reference/ldsc/w_hm3.snplist}"
SEED="${A1_MIXER_SEED:-123}"

mkdir -p "${OUTPUT_DIR}"
for required in "${PYTHON_BIN}" "${MIXER_PY}" "${LIB_BGMG}" "${EXTRACT}"; do
  [[ -e "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 1; }
done

declare -A FILES=(
  [AD_Wightman]="${INPUT_DIR}/AD_Wightman.mixer.gz"
  [AMD_Dry]="${INPUT_DIR}/AMD_Dry.mixer.gz"
  [AMD_Wet]="${INPUT_DIR}/AMD_Wet.mixer.gz"
  [AMD_Any]="${INPUT_DIR}/AMD_Any.mixer.gz"
)

run_fit1() {
  local trait="$1"
  local out="${OUTPUT_DIR}/${trait}"
  [[ -s "${out}.json" ]] && { echo "Reusing completed ${out}.json"; return; }
  "${PYTHON_BIN}" "${MIXER_PY}" fit1 \
    --bim-file "${BIM}" \
    --ld-file "${LD}" \
    --trait1-file "${FILES[${trait}]}" \
    --out "${out}" \
    --extract "${EXTRACT}" \
    --chr2use 1-22 \
    --seed "${SEED}" \
    --lib "${LIB_BGMG}"
}

run_fit2() {
  local amd="$1"
  local out="${OUTPUT_DIR}/AD_vs_${amd}"
  [[ -s "${out}.json" ]] && { echo "Reusing completed ${out}.json"; return; }
  "${PYTHON_BIN}" "${MIXER_PY}" fit2 \
    --bim-file "${BIM}" \
    --ld-file "${LD}" \
    --trait1-file "${FILES[AD_Wightman]}" \
    --trait2-file "${FILES[${amd}]}" \
    --trait1-params-file "${OUTPUT_DIR}/AD_Wightman.json" \
    --trait2-params-file "${OUTPUT_DIR}/${amd}.json" \
    --out "${out}" \
    --extract "${EXTRACT}" \
    --chr2use 1-22 \
    --seed "${SEED}" \
    --lib "${LIB_BGMG}"
}

for trait in AD_Wightman AMD_Dry AMD_Wet AMD_Any; do
  run_fit1 "${trait}"
done
for trait in AMD_Dry AMD_Wet AMD_Any; do
  run_fit2 "${trait}"
done

{
  echo "analysis_run_tag=${RUN_TAG}"
  echo "seed=${SEED}"
  echo "MiXeR_python=externally_configured"
  echo "MiXeR_script=externally_configured"
  echo "input_directory=02_genetic_arch/MiXeR/inputs_${RUN_TAG}"
} > "${OUTPUT_DIR}/run_metadata.txt"

echo "MiXeR completed: ${OUTPUT_DIR}"
