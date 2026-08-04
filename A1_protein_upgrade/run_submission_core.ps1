param(
    [switch]$SkipDownload,
    [switch]$RefreshLiterature,
    [switch]$SkipPrimary,
    [switch]$SkipPWAS,
    [switch]$SkipSubmissionValidityExtension
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$UpgradeRoot = $PSScriptRoot
$ProjectRoot = (Resolve-Path (Join-Path $UpgradeRoot '..')).Path
$Rscript = (Get-Command Rscript -ErrorAction Stop).Source

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Label,
        [Parameter(Mandatory = $true)][scriptblock]$Action
    )
    Write-Host "`n== $Label ==" -ForegroundColor Cyan
    & $Action
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

Set-Location $ProjectRoot

if (-not $SkipPrimary) {
    $primaryArgs = @('-ExecutionPolicy', 'Bypass', '-File', (Join-Path $UpgradeRoot 'run_all.ps1'))
    if ($SkipDownload) { $primaryArgs += '-SkipDownload' }
    if ($RefreshLiterature) { $primaryArgs += '-RefreshLiterature' }
    Invoke-Checked 'Primary literature-prioritized protein workflow' {
        & powershell @primaryArgs
    }
}

if (-not $SkipPWAS) {
    Invoke-Checked 'Prespecified five-protein PWAS crosswalk' {
        & powershell -ExecutionPolicy Bypass -File (Join-Path $UpgradeRoot 'run_PWAS2026_extension.ps1')
    }
}

Invoke-Checked 'Submission-facing supplementary tables' {
    & $Rscript (Join-Path $UpgradeRoot 'scripts\39_prepare_submission_facing_supplement_tables.R')
}

Invoke-Checked 'Targeted pre-submission sensitivities' {
    & $Rscript (Join-Path $UpgradeRoot 'scripts\40_targeted_submission_sensitivities.R')
}

if (-not $SkipSubmissionValidityExtension) {
    $extensionEntry = Join-Path $UpgradeRoot 'extensions\submission_validity_20260717\run_extension.ps1'
    $extensionArgs = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $extensionEntry)
    if ($SkipDownload) { $extensionArgs += '-SkipDownloads' }
    Invoke-Checked 'Kunkle, PAV and exact-assay external candidate sensitivities' {
        & powershell @extensionArgs
    }
}

$deCODEExportScript = Join-Path $ProjectRoot 'A1_dual_scale_mediation\scripts\18_export_decode_submission_assets.py'
if (Test-Path -LiteralPath $deCODEExportScript) {
    Invoke-Checked 'Sanitized deCODE submission tables and source data' {
        & py -3.12 $deCODEExportScript
    }
}

Invoke-Checked 'Submission consistency QA' {
    & powershell -ExecutionPolicy Bypass -File (Join-Path $UpgradeRoot 'run_submission_QA.ps1')
}

Invoke-Checked 'Public release safety QA' {
    & py -3.12 (Join-Path $ProjectRoot 'workflow\validate_release.py')
}

Write-Host "`nA1 submission-facing core workflow completed." -ForegroundColor Green
