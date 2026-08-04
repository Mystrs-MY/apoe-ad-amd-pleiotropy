param(
    [switch]$SkipDownloads
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$scripts = Join-Path $root 'scripts'
$rawPrefix = 'decode_raw_gated_10074_128_8687_26'
$rawCandidates = Join-Path $root "data_processed\${rawPrefix}_pqtl_candidates_for_clumping.tsv"

function Invoke-Step {
    param([string]$Label, [scriptblock]$Action)
    Write-Host "[$Label]"
    & $Action
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE" }
}

Invoke-Step 'Build public-supplement assay gate' {
    python (Join-Path $scripts '06_build_decode_gate.py')
}
Invoke-Step 'Inventory exact deCODE objects' {
    python (Join-Path $scripts '07_inventory_decode_target_objects.py')
}
if (-not $SkipDownloads) {
    Invoke-Step 'Download/checksum nine SMP objects' {
        python (Join-Path $scripts '08_download_decode_target_objects.py') --normalization smp
    }
}
Invoke-Step 'Extract SMP direct alpha and pQTL candidates' {
    pwsh -File (Join-Path $scripts '09_extract_decode_alpha_pqtl.ps1') -Normalization smp
}
Invoke-Step 'Map SMP candidates and prepare instruments' {
    python (Join-Path $scripts '10_prepare_decode_smp_instruments.py')
}
Invoke-Step 'Re-estimate SMP genome-wide beta' {
    Rscript (Join-Path $scripts '11_reestimate_decode_smp_beta.R') --mode=genome_wide
}
Invoke-Step 'Re-estimate SMP cis-only beta' {
    Rscript (Join-Path $scripts '11_reestimate_decode_smp_beta.R') --mode=cis_only
}
Invoke-Step 'Build SMP alpha and initial mediation' {
    python (Join-Path $scripts '12_build_decode_alpha_and_mediation.py')
}
Invoke-Step 'Audit SMP PAV and shared instruments' {
    python (Join-Path $scripts '13_audit_decode_smp_instruments.py')
}
Invoke-Step 'Run SMP target-gene PAV-filtered cis sensitivity' {
    Rscript (Join-Path $scripts '14_pav_filtered_decode_cis_sensitivity.R')
}
Invoke-Step 'Rebuild final SMP mediation' {
    python (Join-Path $scripts '12_build_decode_alpha_and_mediation.py')
}

if (-not $SkipDownloads) {
    Invoke-Step 'Download/checksum two gated raw objects' {
        python (Join-Path $scripts '08_download_decode_target_objects.py') --normalization raw --assays 10074_128 8687_26
    }
}
Invoke-Step 'Extract two gated raw assays' {
    & (Join-Path $scripts '09_extract_decode_alpha_pqtl.ps1') -Normalization raw -Assays @('10074_128', '8687_26')
}
Invoke-Step 'Map gated raw candidates' {
    python (Join-Path $scripts '15_prepare_decode_raw_gated_instruments.py')
}
Invoke-Step 'Re-estimate gated raw genome-wide beta' {
    Rscript (Join-Path $scripts '11_reestimate_decode_smp_beta.R') --mode=genome_wide `
        "--candidate-file=$($rawCandidates.Replace('\', '/'))" "--output-prefix=$rawPrefix" `
        --normalization-label=raw_non_normalized --planned-assays=2 `
        --inference-role=outcome_triggered_normalization_robustness_not_independent_confirmation
}
Invoke-Step 'Re-estimate gated raw cis-only beta' {
    Rscript (Join-Path $scripts '11_reestimate_decode_smp_beta.R') --mode=cis_only `
        "--candidate-file=$($rawCandidates.Replace('\', '/'))" "--output-prefix=$rawPrefix" `
        --normalization-label=raw_non_normalized --planned-assays=2 `
        --inference-role=outcome_triggered_normalization_robustness_not_independent_confirmation
}
Invoke-Step 'Compare raw and SMP normalizations' {
    python (Join-Path $scripts '16_compare_decode_raw_smp.py')
}
Invoke-Step 'Refresh attrition provenance' {
    python (Join-Path $scripts '01_build_phase0_assets.py')
}
Invoke-Step 'Validate design and completed deCODE gate' {
    python (Join-Path $scripts '03_validate_design.py')
}

Write-Host 'deCODE same-platform extension completed.'
