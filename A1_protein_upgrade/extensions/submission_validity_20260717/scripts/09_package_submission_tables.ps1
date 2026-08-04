[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path (Split-Path -Parent $MyInvocation.MyCommand.Path) '..')).Path
$source = Join-Path $root 'tables'
$projectRoot = (Resolve-Path (Join-Path $root '..\..\..')).Path
$destination = Join-Path $projectRoot 'tables_submission\supplementary_tables'
New-Item -ItemType Directory -Path $destination -Force | Out-Null

$mapping = [ordered]@{
    'Kunkle_AD_instrument_coverage.tsv'                    = 'TableS33a_Kunkle_AD_Instrument_Coverage.tsv'
    'Kunkle_AD_protein_beta_sensitivity.tsv'               = 'TableS33b_Kunkle_AD_Protein_Beta.tsv'
    'Kunkle_AD_two_step_mediation_sensitivity.tsv'         = 'TableS33c_Kunkle_AD_Two_Step_Mediation.tsv'
    'Kunkle_AD_aggregate_comparison_with_Wightman.tsv'     = 'TableS33d_Kunkle_vs_Wightman_Aggregate_Comparison.tsv'
    'PAV_epitope_instrument_audit.tsv'                     = 'TableS34a_PAV_Epitope_Instrument_Audit.tsv'
    'PAV_filtered_protein_beta_sensitivity.tsv'            = 'TableS34b_PAV_Filtered_Protein_Beta.tsv'
    'PAV_filtered_two_step_mediation_sensitivity.tsv'      = 'TableS34c_PAV_Filtered_Two_Step_Mediation.tsv'
    'PAV_filtered_aggregate_comparison.tsv'                = 'TableS34d_PAV_Filtered_Aggregate_Comparison.tsv'
    'Lu2026_cross_dataset_candidate_feasibility.tsv'       = 'TableS35a_Lu2026_Candidate_Feasibility.tsv'
    'Lu2026_exact_assay_extension_gate.tsv'                = 'TableS35b_Lu2026_Exact_Assay_Gate.tsv'
    'Lu2026_APOE_variant_to_protein_alpha.tsv'             = 'TableS35c_Lu2026_Exact_Assay_APOE_Alpha.tsv'
    'Lu2026_exact_assay_beta_sensitivity.tsv'              = 'TableS36a_Lu2026_Exact_Assay_Protein_Beta.tsv'
    'Lu2026_exact_assay_two_step_mediation_sensitivity.tsv'= 'TableS36b_Lu2026_Exact_Assay_Two_Step_Mediation.tsv'
    'Lu2026_exact_assay_run_summary.tsv'                   = 'TableS36c_Lu2026_Exact_Assay_Run_Summary.tsv'
}

$manifest = foreach ($entry in $mapping.GetEnumerator()) {
    $sourcePath = Join-Path $source $entry.Key
    if (-not (Test-Path -LiteralPath $sourcePath)) {
        throw "Required extension table is missing: $sourcePath"
    }
    $destinationPath = Join-Path $destination $entry.Value
    Copy-Item -LiteralPath $sourcePath -Destination $destinationPath -Force
    $hash = Get-FileHash -Algorithm SHA256 -LiteralPath $destinationPath
    [pscustomobject]@{
        source_file = $entry.Key
        submission_file = $entry.Value
        bytes = (Get-Item -LiteralPath $destinationPath).Length
        sha256 = $hash.Hash
    }
}

$manifest | Export-Csv -LiteralPath (Join-Path $root 'logs\submission_table_manifest.tsv') -Delimiter "`t" -NoTypeInformation
Write-Host "Packaged $($manifest.Count) extension tables into $destination"
