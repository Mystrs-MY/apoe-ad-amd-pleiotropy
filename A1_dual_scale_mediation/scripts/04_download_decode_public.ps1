param(
    [string]$OutputDir = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$workRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $workRoot "data_raw\decode_public"
}
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

$url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10567571/supplementaryFiles"
$zipPath = Join-Path $OutputDir "PMC10567571_supplementaryFiles.zip"
$partPath = "$zipPath.part"
$xlsxName = "41586_2023_6563_MOESM3_ESM.xlsx"
$xlsxPath = Join-Path $OutputDir $xlsxName
$expectedZipSha256 = "44a3123715e58b2306742ba56b041d6778c2b15c29a287b370275d4c17b47482"
$expectedXlsxSha256 = "280d6acb1c0c0dd2dafb975b143744608358e88da8fb4527e19de94546f96042"

function Assert-Hash {
    param([string]$Path, [string]$Expected)
    $actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $Expected) {
        throw "SHA-256 mismatch for $Path. Expected $Expected; observed $actual"
    }
}

if (Test-Path -LiteralPath $zipPath -PathType Leaf) {
    Assert-Hash -Path $zipPath -Expected $expectedZipSha256
} else {
    & curl.exe -L --fail --retry 8 --retry-all-errors --retry-delay 5 `
        --continue-at - --output $partPath $url
    if ($LASTEXITCODE -ne 0) {
        throw "Europe PMC download failed with exit code $LASTEXITCODE; resumable file retained at $partPath"
    }
    Assert-Hash -Path $partPath -Expected $expectedZipSha256
    Move-Item -LiteralPath $partPath -Destination $zipPath
}

if (Test-Path -LiteralPath $xlsxPath -PathType Leaf) {
    Assert-Hash -Path $xlsxPath -Expected $expectedXlsxSha256
} else {
    & tar -xf $zipPath -C $OutputDir $xlsxName
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to extract $xlsxName"
    }
    Assert-Hash -Path $xlsxPath -Expected $expectedXlsxSha256
}

Write-Output "deCODE public supplement ready: $xlsxPath"
Write-Output "Source: $url"
Write-Output "ZIP SHA-256: $expectedZipSha256"
Write-Output "XLSX SHA-256: $expectedXlsxSha256"
