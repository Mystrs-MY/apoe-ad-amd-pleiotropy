param(
    [switch]$SkipMain,
    [switch]$SkipSupplementary
)

$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
$env:A1_PROJECT_ROOT = $Root

function Invoke-RFigure {
    param([string]$RelativePath)
    $script = Join-Path $Root $RelativePath
    if (-not (Test-Path -LiteralPath $script)) { throw "Missing figure script: $script" }
    Write-Host "Running $RelativePath" -ForegroundColor Cyan
    Push-Location (Split-Path -Parent $script)
    try {
        & Rscript $script
        if ($LASTEXITCODE -ne 0) { throw "Figure script failed: $RelativePath" }
    }
    finally {
        Pop-Location
    }
}

if (-not $SkipMain) {
    Invoke-RFigure '01_meta\redraw_figure2_vector.R'
}

if (-not $SkipSupplementary) {
    Write-Host 'Supplementary multi-panel layout scripts are intentionally excluded; use the committed source-data tables for numerical verification.' -ForegroundColor DarkGray
}

Write-Host 'Figure workflow completed.' -ForegroundColor Green
