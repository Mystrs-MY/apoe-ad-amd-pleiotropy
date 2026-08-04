$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$ExtensionRoot = Join-Path $ProjectRoot 'A1_protein_upgrade\extensions\PWAS2026_crosswalk_extension'
$ScriptRoot = Join-Path $ExtensionRoot 'scripts'
$LogRoot = Join-Path $ExtensionRoot 'logs\pipeline_runs'
$Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
if (-not $Rscript) {
    $Rscript = Get-ChildItem 'C:\Program Files\R' -Filter Rscript.exe -Recurse -ErrorAction Stop |
        Sort-Object FullName -Descending | Select-Object -First 1 -ExpandProperty FullName
}
$Python = (Get-Command python -ErrorAction Stop).Source
$RunStamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$RunLog = Join-Path $LogRoot "PWAS2026_crosswalk_$RunStamp.log"

New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null

function Invoke-LoggedStep {
    param(
        [Parameter(Mandatory = $true)][string]$Label,
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )

    "[$(Get-Date -Format o)] START $Label" | Tee-Object -FilePath $RunLog -Append
    & $Executable @Arguments 2>&1 | Tee-Object -FilePath $RunLog -Append
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
    "[$(Get-Date -Format o)] END $Label" | Tee-Object -FilePath $RunLog -Append
}

"A1 frozen five-protein PWAS 2026 crosswalk extension" | Set-Content -LiteralPath $RunLog -Encoding UTF8
"ProjectRoot=$ProjectRoot" | Add-Content -LiteralPath $RunLog -Encoding UTF8
"Seed=20260714" | Add-Content -LiteralPath $RunLog -Encoding UTF8

Invoke-LoggedStep '01 audit frozen archives' $Python @((Join-Path $ScriptRoot '01_audit_PWAS5_archives.py'))
Invoke-LoggedStep '02 extract APOE alpha and pQTL candidates' $Python @((Join-Path $ScriptRoot '02_extract_PWAS5_alpha_and_pqtl.py'))
Invoke-LoggedStep '03 protein-to-outcome main instruments' $Rscript @((Join-Path $ScriptRoot '03_run_PWAS5_beta_main.R'))
Invoke-LoggedStep '04 protein-to-outcome cis instruments' $Rscript @((Join-Path $ScriptRoot '04_run_PWAS5_beta_cis.R'))
Invoke-LoggedStep '05 two-step mediation' $Rscript @((Join-Path $ScriptRoot '05_run_PWAS5_mediation.R'))
Invoke-LoggedStep '06 incremental aggregate sensitivity' $Rscript @((Join-Path $ScriptRoot '06_run_PWAS5_incremental_aggregate.R'))
Invoke-LoggedStep '07 plot source panels for composite Fig. S11' $Rscript @((Join-Path $ScriptRoot '07_plot_PWAS5_extension.R'), $ProjectRoot)

"[$(Get-Date -Format o)] PIPELINE COMPLETE" | Add-Content -LiteralPath $RunLog -Encoding UTF8
Write-Host "PWAS5 pipeline complete. Log: $RunLog"
