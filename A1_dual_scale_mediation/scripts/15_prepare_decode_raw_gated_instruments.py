#!/usr/bin/env python3
"""Reuse the audited mapping workflow for the gated two-assay raw sensitivity."""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_SCRIPT = Path(__file__).with_name("10_prepare_decode_smp_instruments.py")
spec = importlib.util.spec_from_file_location("decode_mapping_workflow", SOURCE_SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load {SOURCE_SCRIPT}")
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)

prefix = "decode_raw_gated_10074_128_8687_26"
workflow.CANDIDATES = ROOT / "data_processed" / f"{prefix}_pqtl_candidates_p5e8.tsv"
workflow.MAPPED_OUTPUT = ROOT / "data_processed" / f"{prefix}_pqtl_candidates_mapped_qc.tsv.gz"
workflow.ELIGIBLE_OUTPUT = ROOT / "data_processed" / f"{prefix}_pqtl_candidates_for_clumping.tsv"
workflow.SUMMARY_OUTPUT = ROOT / "tables" / f"{prefix}_instrument_mapping_qc.tsv"
workflow.APOE_AUDIT_OUTPUT = ROOT / "tables" / f"{prefix}_coordinate_and_APOE_exclusion_audit.tsv"
workflow.LOG_OUTPUT = ROOT / "logs" / f"{prefix}_instrument_mapping.log"

workflow.main()
