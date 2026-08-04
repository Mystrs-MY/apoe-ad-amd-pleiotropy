#!/usr/bin/env python3
"""Extract direct APOE-variant alpha estimates for the literature-prioritized panel."""

from __future__ import annotations

import gzip
import io
import math
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "config" / "resources.yml").read_text(encoding="utf-8"))
RESOURCE_DIRS = [Path(path) for path in CONFIG["paths"]["ukb_ppp_tar_dirs"]]
PROVENANCE = ROOT / "tables" / "Table_Literature_Prioritized_Protein_Provenance.tsv"
OUTPUT = ROOT / "tables" / "APOE_variant_to_literature_proteins_alpha.tsv"
MANIFEST = ROOT / "data_processed" / "ukbppp_local_assay_manifest.tsv"
TARGET_MAPPING = ROOT / "tables" / "APOE_linkable_target_assay_mapping.tsv"

VARIANTS = {
    "rs429358": {"chr": "19", "position_hg19": 45411941, "effect_allele": "C", "other_allele": "T"},
    "rs7412": {"chr": "19", "position_hg19": 45412079, "effect_allele": "T", "other_allele": "C"},
}


def build_manifest() -> pd.DataFrame:
    rows = []
    pattern = re.compile(
        r"^(?P<gene>[^_]+)_(?P<uniprot>[^_]+)_(?P<oid>OID\d+)_v(?P<version>\d+)_(?P<panel>.+)\.tar$",
        re.IGNORECASE,
    )
    paths = sorted(path for resource_dir in RESOURCE_DIRS for path in resource_dir.glob("*.tar"))
    for path in paths:
        match = pattern.match(path.name)
        if not match:
            rows.append({
                "gene_symbol": "mapping_unresolved", "UniProt_ID": "mapping_unresolved",
                "Olink_target_ID": "mapping_unresolved", "assay_version": "mapping_unresolved",
                "Olink_panel": "mapping_unresolved", "tar_path": str(path), "mapping_status": "filename_parse_failed",
            })
            continue
        values = match.groupdict()
        rows.append({
            "gene_symbol": values["gene"].upper(), "UniProt_ID": values["uniprot"].upper(),
            "Olink_target_ID": values["oid"].upper(), "assay_version": f"v{values['version']}",
            "Olink_panel": values["panel"], "tar_path": str(path), "mapping_status": "parsed_from_source_filename",
            "resource_scope": "synapse_primary_panel_download" if "syn51365303" in str(path) else "legacy_local_resource",
        })
    manifest = pd.DataFrame(rows)
    duplicated = manifest.loc[
        manifest["mapping_status"].eq("parsed_from_source_filename")
        & manifest["gene_symbol"].duplicated(keep=False),
        "gene_symbol",
    ].unique()
    if len(duplicated):
        raise ValueError(
            "Multiple Olink assays found for gene symbols; assay-level mapping is required: "
            + ", ".join(sorted(duplicated))
        )
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, sep="\t", index=False)
    return manifest


def read_variant_rows(tar_path: Path) -> dict[str, dict[str, str]]:
    found: dict[str, dict[str, str]] = {}
    with tarfile.open(tar_path, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile() and "chr19_" in member.name]
        if len(members) != 1:
            raise ValueError(f"Expected one chr19 member in {tar_path.name}; found {len(members)}")
        member_stream = archive.extractfile(members[0])
        if member_stream is None:
            raise OSError(f"Could not read {members[0].name}")
        with gzip.GzipFile(fileobj=member_stream) as compressed:
            text = io.TextIOWrapper(compressed, encoding="utf-8")
            header = text.readline().strip().split()
            required = {"CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "LOG10P"}
            if not required.issubset(header):
                raise ValueError(f"Unexpected UKB-PPP columns in {members[0].name}")
            for line in text:
                values = dict(zip(header, line.strip().split()))
                identifier_parts = values["ID"].split(":")
                if len(identifier_parts) < 2:
                    continue
                position_hg19_from_id = int(identifier_parts[1])
                for rsid, specification in VARIANTS.items():
                    if position_hg19_from_id == specification["position_hg19"]:
                        found[rsid] = values
                if len(found) == len(VARIANTS):
                    break
    return found


def p_from_log10(value: str) -> float:
    log10p = float(value)
    return 0.0 if log10p > 323 else 10 ** (-log10p)


def verified_values(values: pd.Series) -> set[str]:
    invalid = {"", "NA", "NOT_REPORTED", "MAPPING_UNRESOLVED"}
    return {str(value).strip().upper() for value in values if str(value).strip().upper() not in invalid}


def build_target_mapping(evidence: pd.DataFrame, manifest_by_gene: dict[str, pd.Series]) -> pd.DataFrame:
    primary = evidence[evidence["inclusion_status"] == "primary_high_confidence_panel"].copy()
    key_columns = ["gene_symbol", "protein_name", "protein_form_or_isoform"]
    targets = primary[key_columns].drop_duplicates().sort_values(key_columns).reset_index(drop=True)
    rows = []
    for index, target in targets.iterrows():
        mask = np.logical_and.reduce([primary[column].eq(target[column]) for column in key_columns])
        source = primary.loc[mask]
        gene = target["gene_symbol"]
        assay = manifest_by_gene.get(gene)
        source_oids = verified_values(source["assay_target_ID"])
        source_uniprots = verified_values(source["UniProt_ID"])
        source_platforms = ";".join(sorted(verified_values(source["proteomic_platform"])))
        isoform_specific = "isoform" in (
            f"{target['protein_name']} {target['protein_form_or_isoform']}".lower()
        )
        canonical_target = str(target["protein_name"]).strip().upper() == gene
        if assay is None:
            status = "assay_unavailable"
            confidence = "not_mapped"
            basis = "no_exact_gene_prefix_Olink_assay"
            eligible = False
        else:
            local_oid = str(assay["Olink_target_ID"]).upper()
            local_uniprot = str(assay["UniProt_ID"]).upper()
            exact_oid = local_oid in source_oids
            exact_uniprot = local_uniprot in source_uniprots
            same_olink_platform = any("OLINK" in value for value in verified_values(source["proteomic_platform"]))
            if isoform_specific:
                status = "mapping_unresolved_isoform_specific"
                confidence = "not_mapped"
                basis = "local_Olink_assay_is_not_isoform_specific"
                eligible = False
            elif exact_oid:
                status = "eligible_exact_Olink_target_ID"
                confidence = "high"
                basis = "exact_OID_match"
                eligible = True
            elif exact_uniprot and canonical_target:
                status = "eligible_exact_UniProt_cross_platform"
                confidence = "moderate"
                basis = "exact_UniProt_and_canonical_protein;platform_heterogeneity_retained"
                eligible = True
            elif same_olink_platform and canonical_target:
                status = "eligible_same_Olink_platform_canonical_target"
                confidence = "moderate"
                basis = "same_Olink_platform_and_canonical_target;source_assay_ID_not_OID"
                eligible = True
            else:
                status = "mapping_unresolved_cross_platform"
                confidence = "not_mapped"
                basis = "gene_symbol_only_or_source_assay_identity_unresolved"
                eligible = False
        rows.append({
            "target_entry_id": f"L2T{index + 1:03d}",
            "gene_symbol": gene,
            "protein_name": target["protein_name"],
            "protein_form_or_isoform": target["protein_form_or_isoform"],
            "source_record_ids": ";".join(sorted(set(source["record_id"]))),
            "source_PMIDs": ";".join(sorted(verified_values(source["PMID"]))),
            "source_platforms": source_platforms,
            "source_assay_target_IDs": ";".join(sorted(source_oids)) or "not_reported",
            "source_UniProt_IDs": ";".join(sorted(source_uniprots)) or "not_reported",
            "local_Olink_target_ID": assay["Olink_target_ID"] if assay is not None else "mapping_unresolved",
            "local_Olink_UniProt_ID": assay["UniProt_ID"] if assay is not None else "mapping_unresolved",
            "local_Olink_panel": assay["Olink_panel"] if assay is not None else "mapping_unresolved",
            "mapping_status": status,
            "mapping_confidence": confidence,
            "mapping_basis": basis,
            "eligible_for_alpha_beta_triangulation": str(eligible).lower(),
            "manual_verification_status": "rule_based_mapping_verified_from_provenance_and_synapse_filename" if eligible else "requires_manual_assay_review",
        })
    mapping = pd.DataFrame(rows)
    mapping.to_csv(TARGET_MAPPING, sep="\t", index=False)
    return mapping


def main() -> None:
    if not MANIFEST.exists() or not TARGET_MAPPING.exists():
        raise FileNotFoundError("Run 21_build_name_matched_mapping.py before alpha extraction.")
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("NA")
    manifest = manifest[manifest["mapping_status"] == "parsed_from_source_filename"].copy()
    manifest_by_oid = {row["Olink_target_ID"]: row for _, row in manifest.iterrows()}
    evidence = pd.read_csv(PROVENANCE, sep="\t", dtype=str).fillna("NA")
    primary = evidence[evidence["inclusion_status"] == "primary_high_confidence_panel"].copy()
    panel_genes = sorted(primary["gene_symbol"].unique())
    target_mapping = pd.read_csv(TARGET_MAPPING, sep="\t", dtype=str).fillna("NA")
    eligible_mapping = target_mapping[target_mapping["eligible_for_primary"] == "true"].copy()
    eligible_genes = set(eligible_mapping["standardized_gene_symbol"])
    strict_genes = set(target_mapping.loc[
        target_mapping["eligible_for_strict_sensitivity"] == "true", "standardized_gene_symbol"
    ])
    selected_by_gene: dict[str, pd.Series] = {}
    for gene, group in eligible_mapping.groupby("standardized_gene_symbol"):
        selected_oids = sorted(set(group.loc[group["selected_assay"] != "none", "selected_assay"]))
        if len(selected_oids) != 1:
            raise ValueError(f"Expected one selected assay for eligible gene {gene}; found {selected_oids}")
        selected_by_gene[gene] = manifest_by_oid[selected_oids[0]]

    extracted_by_gene: dict[str, dict[str, dict[str, str]]] = {}
    rows = []
    for gene in panel_genes:
        assay = selected_by_gene.get(gene)
        if assay is not None and gene not in extracted_by_gene:
            extracted_by_gene[gene] = read_variant_rows(Path(assay["tar_path"]))
        gene_mapping = eligible_mapping[eligible_mapping["standardized_gene_symbol"] == gene]
        mapping_eligible = gene in eligible_genes
        target_entry_ids = ";".join(gene_mapping["literature_target_entry"]) if mapping_eligible else "NA"
        mapping_statuses = ";".join(sorted(set(gene_mapping["mapping_status"]))) if mapping_eligible else "mapping_unresolved_or_assay_unavailable"
        mapping_confidence = (
            "high" if any(gene_mapping["mapping_confidence"] == "high")
            else "moderate" if mapping_eligible else "not_mapped"
        )
        for rsid, specification in VARIANTS.items():
            source = extracted_by_gene.get(gene, {}).get(rsid)
            if assay is None:
                status = "alpha_unavailable_local_UKB_PPP_subset"
                reason = "No exact gene-prefix Olink assay was available in the verified local and Synapse-selected UKB-PPP resources."
            elif source is None:
                status = "direct_variant_absent_in_assay_GWAS"
                reason = f"{rsid} was not found at the verified hg19 position in the assay chr19 file."
            else:
                status = "direct_variant_available"
                reason = "NA"

            if source is not None:
                original_effect = source["ALLELE1"].upper()
                original_other = source["ALLELE0"].upper()
                requested_effect = specification["effect_allele"]
                if original_effect == requested_effect and original_other == specification["other_allele"]:
                    beta = float(source["BETA"])
                    flip = "false"
                elif original_other == requested_effect and original_effect == specification["other_allele"]:
                    beta = -float(source["BETA"])
                    flip = "true"
                else:
                    beta = math.nan
                    status = "alleles_not_harmonizable"
                    reason = f"Observed {original_effect}/{original_other}; expected {requested_effect}/{specification['other_allele']}."
                    flip = "not_applicable"
                se = float(source["SE"])
                p_value = p_from_log10(source["LOG10P"])
                eaf = float(source["A1FREQ"])
                sample_size = source["N"]
                source_id = source["ID"]
                info = source["INFO"]
            else:
                original_effect = original_other = flip = "NA"
                beta = se = p_value = eaf = math.nan
                sample_size = source_id = info = "NA"

            rows.append({
                "gene_symbol": gene,
                "protein_name": gene,
                "protein_form_or_isoform": "Olink assay-level target",
                "variant": rsid,
                "chromosome": specification["chr"],
                "position_hg19": specification["position_hg19"],
                "position_hg38_from_GENPOS": source["GENPOS"] if source is not None else "NA",
                "requested_effect_allele": specification["effect_allele"],
                "requested_other_allele": specification["other_allele"],
                "original_effect_allele": original_effect,
                "original_other_allele": original_other,
                "allele_flipped": flip,
                "beta": beta,
                "SE": se,
                "P_value": p_value,
                "effect_allele_frequency_original": eaf,
                "sample_size": sample_size,
                "imputation_INFO": info,
                "source_variant_ID": source_id,
                "proxy_used": "false",
                "proxy_r2": "NA",
                "proteomic_platform": "UKB-PPP Olink Explore 3072",
                "UniProt_ID": assay["UniProt_ID"] if assay is not None else "mapping_unresolved",
                "assay_target_ID": assay["Olink_target_ID"] if assay is not None else "mapping_unresolved",
                "assay_version": assay["assay_version"] if assay is not None else "mapping_unresolved",
                "Olink_panel": assay["Olink_panel"] if assay is not None else "mapping_unresolved",
                "effect_unit": "per genetically predicted Olink NPX unit as reported by UKB-PPP",
                "alpha_source": str(assay["tar_path"]) if assay is not None else "NA",
                "availability_status": status,
                "exclusion_reason": reason,
                "eligible_for_two_step_mapping": str(mapping_eligible).lower(),
                "eligible_for_strict_sensitivity": str(gene in strict_genes).lower(),
                "linked_Layer2_target_entry_IDs": target_entry_ids,
                "target_mapping_status": mapping_statuses,
                "target_mapping_confidence": mapping_confidence,
                "manual_verification_status": "direct_hg19_ID_and_hg38_GENPOS_coordinates_plus_alleles_verified_against_local_chr19_rsid_map" if source is not None else "not_available",
            })

    output = pd.DataFrame(rows).replace({np.nan: "NA"})
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)
    print(f"Primary protein target entries: {len(target_mapping)}")
    print(f"Primary unique genes: {len(panel_genes)}")
    print(f"Direct alpha rows: {(output['availability_status'] == 'direct_variant_available').sum()}")
    print(f"Genes with any local assay: {output.loc[output['assay_target_ID'] != 'mapping_unresolved', 'gene_symbol'].nunique()}")
    print(f"Genes with eligible target-to-assay mapping: {output.loc[output['eligible_for_two_step_mapping'] == 'true', 'gene_symbol'].nunique()}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
