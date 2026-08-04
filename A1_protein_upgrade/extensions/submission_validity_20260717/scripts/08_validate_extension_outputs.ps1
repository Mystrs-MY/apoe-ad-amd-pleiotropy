[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path (Split-Path -Parent $MyInvocation.MyCommand.Path) '..')).Path
$tables = Join-Path $root 'tables'
$logPath = Join-Path $root 'logs\extension_validation_summary.tsv'
$checks = [System.Collections.Generic.List[object]]::new()

function Add-Check {
    param(
        [string]$Name,
        [bool]$Passed,
        [string]$Observed,
        [string]$Expected
    )
    $checks.Add([pscustomobject]@{
        check = $Name
        passed = $Passed
        observed = $Observed
        expected = $Expected
    })
    if (-not $Passed) {
        throw "Validation failed: $Name; observed=$Observed; expected=$Expected"
    }
}

function Read-Tsv {
    param([string]$Name)
    $path = Join-Path $tables $Name
    Add-Check "file:$Name" (Test-Path -LiteralPath $path) $path 'file exists'
    $rows = @(Import-Csv -LiteralPath $path -Delimiter "`t")
    Add-Check "rows:$Name" ($rows.Count -gt 0) ([string]$rows.Count) '>0'
    return $rows
}

function Is-Corrected {
    param([string]$Fdr, [string]$Bonferroni)
    if ($Fdr -in @('', 'NA', $null) -or $Bonferroni -in @('', 'NA', $null)) {
        return $false
    }
    return ([double]$Fdr -le 0.05 -or [double]$Bonferroni -le 0.05)
}

$kCoverage = Read-Tsv 'Kunkle_AD_instrument_coverage.tsv'
$kBeta = Read-Tsv 'Kunkle_AD_protein_beta_sensitivity.tsv'
$kMediation = Read-Tsv 'Kunkle_AD_two_step_mediation_sensitivity.tsv'
$kAggregate = Read-Tsv 'Kunkle_AD_aggregate_comparison_with_Wightman.tsv'
$pavAudit = Read-Tsv 'PAV_epitope_instrument_audit.tsv'
$pavBeta = Read-Tsv 'PAV_filtered_protein_beta_sensitivity.tsv'
$pavMediation = Read-Tsv 'PAV_filtered_two_step_mediation_sensitivity.tsv'
$pavAggregate = Read-Tsv 'PAV_filtered_aggregate_comparison.tsv'
$luFeasibility = Read-Tsv 'Lu2026_cross_dataset_candidate_feasibility.tsv'
$luGate = Read-Tsv 'Lu2026_exact_assay_extension_gate.tsv'
$luAlpha = Read-Tsv 'Lu2026_APOE_variant_to_protein_alpha.tsv'
$luBeta = Read-Tsv 'Lu2026_exact_assay_beta_sensitivity.tsv'
$luMediation = Read-Tsv 'Lu2026_exact_assay_two_step_mediation_sensitivity.tsv'
$luSummary = Read-Tsv 'Lu2026_exact_assay_run_summary.tsv'

$kPrimaryBeta = @($kBeta | Where-Object { $_.method_role -eq 'primary' })
$kProteinMediation = @($kMediation | Where-Object { $_.row_type -eq 'protein' })
Add-Check 'Kunkle protein coverage rows' ($kCoverage.Count -eq 25) $kCoverage.Count '25'
Add-Check 'Kunkle prespecified instrument-pair total' ((($kCoverage | Measure-Object prespecified_instrument_pairs -Sum).Sum) -eq 409) (($kCoverage | Measure-Object prespecified_instrument_pairs -Sum).Sum) '409'
Add-Check 'Kunkle instrument pairs present in outcome' ((($kCoverage | Measure-Object present_in_Kunkle -Sum).Sum) -eq 385) (($kCoverage | Measure-Object present_in_Kunkle -Sum).Sum) '385'
Add-Check 'Kunkle instruments retained after harmonization' ((($kCoverage | Measure-Object retained_after_harmonization -Sum).Sum) -eq 335) (($kCoverage | Measure-Object retained_after_harmonization -Sum).Sum) '335'
Add-Check 'Kunkle primary protein beta rows' ($kPrimaryBeta.Count -eq 25) $kPrimaryBeta.Count '25'
Add-Check 'Kunkle protein mediation rows' ($kProteinMediation.Count -eq 50) $kProteinMediation.Count '50'
Add-Check 'Kunkle corrected protein beta count' (@($kPrimaryBeta | Where-Object { Is-Corrected $_.P_FDR_observed $_.P_Bonferroni_planned }).Count -eq 0) (@($kPrimaryBeta | Where-Object { Is-Corrected $_.P_FDR_observed $_.P_Bonferroni_planned }).Count) '0'
Add-Check 'Kunkle corrected mediation count' (@($kProteinMediation | Where-Object { Is-Corrected $_.indirect_P_FDR $_.indirect_P_Bonferroni }).Count -eq 0) (@($kProteinMediation | Where-Object { Is-Corrected $_.indirect_P_FDR $_.indirect_P_Bonferroni }).Count) '0'

$unresolved = @($pavAudit | Where-Object { $_.annotation_status -like 'unresolved*' })
$highRisk = @($pavAudit | Where-Object { $_.high_confidence_epitope_risk -eq 'TRUE' })
Add-Check 'PAV audit instrument-pair rows' ($pavAudit.Count -eq 409) $pavAudit.Count '409'
Add-Check 'PAV unresolved rows' ($unresolved.Count -eq 332) $unresolved.Count '332'
Add-Check 'PAV high-confidence risk rows' ($highRisk.Count -eq 2) $highRisk.Count '2'
$highRiskKeys = @($highRisk | ForEach-Object { "$($_.gene_symbol):$($_.SNP)" } | Sort-Object)
Add-Check 'PAV high-confidence risk identities' (($highRiskKeys -join ',') -eq 'CSF2:rs25882,TREM2:rs188904277') ($highRiskKeys -join ',') 'CSF2:rs25882,TREM2:rs188904277'

$csf2 = @($pavBeta | Where-Object { $_.method_role -eq 'primary' -and $_.gene_symbol -eq 'CSF2' -and $_.outcome -eq 'wet_AMD' })
Add-Check 'PAV-filtered CSF2 wet-AMD row' ($csf2.Count -eq 1) $csf2.Count '1'
Add-Check 'PAV-filtered CSF2 wet-AMD nominal P' ([math]::Abs([double]$csf2[0].P_value - 0.0919267565207214) -lt 1e-12) $csf2[0].P_value '0.0919267565207214'
$trem2Path = @($pavMediation | Where-Object { $_.row_type -eq 'protein' -and $_.variant -eq 'rs429358' -and $_.gene_symbol -eq 'TREM2' -and $_.outcome -eq 'AD' })
Add-Check 'PAV-filtered TREM2 mediation row' ($trem2Path.Count -eq 1) $trem2Path.Count '1'
Add-Check 'PAV-filtered TREM2 planned Bonferroni P' ([math]::Abs([double]$trem2Path[0].indirect_P_Bonferroni_planned - 0.0292058775597498) -lt 1e-12) $trem2Path[0].indirect_P_Bonferroni_planned '0.0292058775597498'

$luSourceCandidates = @($luFeasibility | Where-Object { $_.Lu_candidate_status -eq 'candidate_supported_in_at_least_two_dataset_outcomes_and_two_independent_datasets' })
$luEligible = @($luFeasibility | Where-Object { $_.eligible_for_Lu2026_exact_assay_extension -eq 'TRUE' })
Add-Check 'Lu source candidates' ($luSourceCandidates.Count -eq 8) $luSourceCandidates.Count '8'
Add-Check 'Lu source candidates supported in at least two independent datasets' (@($luSourceCandidates | Where-Object { [int]$_.n_supporting_datasets -ge 2 }).Count -eq 8) (@($luSourceCandidates | Where-Object { [int]$_.n_supporting_datasets -ge 2 }).Count) '8'
Add-Check 'Lu exact-assay candidates' ($luEligible.Count -eq 5) $luEligible.Count '5'
Add-Check 'Lu direct APOE alpha rows' ($luAlpha.Count -eq 10) $luAlpha.Count '10'
Add-Check 'Lu run-summary layers' ($luSummary.Count -eq 2) $luSummary.Count '2'
foreach ($row in $luSummary) {
    Add-Check "Lu beta estimability:$($row.analysis_set)" ([int]$row.estimable_beta_tests -eq 20) $row.estimable_beta_tests '20'
    Add-Check "Lu mediation estimability:$($row.analysis_set)" ([int]$row.mediation_tests -eq 40) $row.mediation_tests '40'
    Add-Check "Lu corrected beta signals:$($row.analysis_set)" ([int]$row.beta_FDR_signals -eq 0 -and [int]$row.beta_Bonferroni_signals -eq 0) "$($row.beta_FDR_signals)/$($row.beta_Bonferroni_signals)" '0/0'
    Add-Check "Lu corrected mediation signals:$($row.analysis_set)" ([int]$row.mediation_FDR_signals -eq 0 -and [int]$row.mediation_Bonferroni_signals -eq 0) "$($row.mediation_FDR_signals)/$($row.mediation_Bonferroni_signals)" '0/0'
    Add-Check "Lu primary panel unchanged:$($row.analysis_set)" ($row.primary_panel_modified -eq 'FALSE') $row.primary_panel_modified 'FALSE'
}

$checks | Export-Csv -LiteralPath $logPath -Delimiter "`t" -NoTypeInformation
Write-Host "Validated $($checks.Count) extension output checks."
