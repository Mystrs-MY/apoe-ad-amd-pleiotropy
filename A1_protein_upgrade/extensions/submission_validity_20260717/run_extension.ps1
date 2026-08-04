[CmdletBinding()]
param(
    [switch]$SkipDownloads,
    [switch]$SkipKunkle,
    [switch]$SkipLu
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$upgradeRoot = (Resolve-Path (Join-Path $root '..\..')).Path
$scripts = Join-Path $root 'scripts'
$luRaw = Join-Path $root 'data_raw\lu2026_ukbppp_exact_assays'
$luInventory = Join-Path $root 'data_raw\synapse_inventory.tsv'
$luSelection = Join-Path $root 'config\lu2026_exact_assay_download_selection.tsv'
$synapseDownloader = Join-Path $upgradeRoot 'scripts\14_inventory_and_download_synapse_ukbppp.py'

function Invoke-Checked {
    param([string]$Program, [string[]]$Arguments)
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code $LASTEXITCODE`: $Program $($Arguments -join ' ')"
    }
}

if (-not $SkipDownloads) {
    Invoke-Checked 'powershell' @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $scripts '01_download_kunkle_stage1.ps1'))
    Invoke-Checked 'powershell' @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $scripts '02_download_external_supplements.ps1'))
}

if (-not $SkipKunkle) {
    Invoke-Checked 'Rscript' @((Join-Path $scripts '03_kunkle_ad_sensitivity.R'))
}

Invoke-Checked 'Rscript' @((Join-Path $scripts '04_pav_epitope_audit.R'))
Invoke-Checked 'Rscript' @((Join-Path $scripts '05_lu2026_feasibility_gate.R'))

if (-not $SkipLu) {
    $selectionRows = Import-Csv -LiteralPath $luSelection -Delimiter "`t" |
        Where-Object { $_.selection_status -eq 'selected_for_download' }
    if ($selectionRows.Count -ne 5) {
        throw "Expected exactly five gate-passing Lu 2026 assays; observed $($selectionRows.Count)."
    }

    $verifiedTars = @(Get-ChildItem -LiteralPath $luRaw -Filter '*.tar' -File -ErrorAction SilentlyContinue)
    if ($verifiedTars.Count -lt 5) {
        if (-not $env:SYNAPSE_AUTH_TOKEN) {
            throw 'Five verified Lu assay archives are not present and SYNAPSE_AUTH_TOKEN is not set.'
        }
        Invoke-Checked 'python' @(
            $synapseDownloader,
            '--entity', 'syn51365303',
            '--output', $luRaw,
            '--inventory', $luInventory,
            '--select', $luSelection,
            '--reuse-inventory'
        )
    }

    Invoke-Checked 'python' @((Join-Path $scripts '06_prepare_lu2026_ukbppp.py'))
    Invoke-Checked 'Rscript' @((Join-Path $scripts '07_lu2026_exact_assay_sensitivity.R'))
}

if (-not $SkipLu) {
    Invoke-Checked 'powershell' @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $scripts '08_validate_extension_outputs.ps1'))
    Invoke-Checked 'powershell' @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $scripts '09_package_submission_tables.ps1'))
}

Write-Host 'A1 bounded submission-validity extension completed.'
