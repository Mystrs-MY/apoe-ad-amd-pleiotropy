param(
    [string]$DestinationRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'data_raw\kunkle_2019'),
    [int]$RetryCount = 100
)

$ErrorActionPreference = "Stop"
$url = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt"
$expectedBytes = 569380789L
$final = Join-Path $DestinationRoot "Kunkle_etal_Stage1_results.txt"
$partial = "$final.part"
$log = Join-Path $DestinationRoot "download.log"

New-Item -ItemType Directory -Force -Path $DestinationRoot | Out-Null

if ((Test-Path -LiteralPath $final) -and ((Get-Item -LiteralPath $final).Length -eq $expectedBytes)) {
    "$(Get-Date -Format s) already complete: $final" | Add-Content -LiteralPath $log
    exit 0
}

if ((Test-Path -LiteralPath $partial) -and ((Get-Item -LiteralPath $partial).Length -eq $expectedBytes)) {
    Move-Item -LiteralPath $partial -Destination $final -Force
    $hash = (Get-FileHash -LiteralPath $final -Algorithm SHA256).Hash
    "$(Get-Date -Format s) completed existing partial bytes=$expectedBytes sha256=$hash" | Add-Content -LiteralPath $log
    exit 0
}

"$(Get-Date -Format s) starting/resuming $url" | Add-Content -LiteralPath $log
$curlLog = Join-Path $DestinationRoot "curl_stderr.log"
$curlArgs = @(
    "--fail", "--location", "--continue-at", "-", "--retry", "$RetryCount", "--retry-all-errors",
    "--retry-delay", "5", "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "60",
    "--output", $partial, $url
)
$curlProcess = Start-Process -FilePath "curl.exe" -ArgumentList $curlArgs -Wait -PassThru -NoNewWindow `
    -RedirectStandardError $curlLog
if (Test-Path -LiteralPath $curlLog) {
    Get-Content -LiteralPath $curlLog | Add-Content -LiteralPath $log
}
if ($curlProcess.ExitCode -ne 0) {
    "$(Get-Date -Format s) curl exit code $($curlProcess.ExitCode)" | Add-Content -LiteralPath $log
    exit $curlProcess.ExitCode
}

$actualBytes = (Get-Item -LiteralPath $partial).Length
if ($actualBytes -ne $expectedBytes) {
    "$(Get-Date -Format s) size mismatch: expected=$expectedBytes actual=$actualBytes" | Add-Content -LiteralPath $log
    exit 2
}

Move-Item -LiteralPath $partial -Destination $final -Force
$hash = (Get-FileHash -LiteralPath $final -Algorithm SHA256).Hash
"$(Get-Date -Format s) complete bytes=$actualBytes sha256=$hash" | Add-Content -LiteralPath $log
