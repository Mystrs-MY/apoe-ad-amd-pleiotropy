param(
    [string]$DestinationRoot = (Join-Path (Split-Path -Parent $PSScriptRoot) 'data_raw'),
    [int]$RetryCount = 100
)

$ErrorActionPreference = "Stop"

function Test-ZipArchive {
    param([Parameter(Mandatory = $true)][string]$Path)

    try {
        Add-Type -AssemblyName System.IO.Compression.FileSystem
        $archive = [System.IO.Compression.ZipFile]::OpenRead($Path)
        $entryCount = $archive.Entries.Count
        $archive.Dispose()
        return ($entryCount -gt 0)
    }
    catch {
        return $false
    }
}

function Get-ResumableZip {
    param(
        [Parameter(Mandatory = $true)][string]$Url,
        [Parameter(Mandatory = $true)][string]$OutputDirectory,
        [Parameter(Mandatory = $true)][string]$FileName
    )

    New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
    $final = Join-Path $OutputDirectory $FileName
    $partial = "$final.part"
    $log = Join-Path $OutputDirectory "download.log"

    if ((Test-Path -LiteralPath $final) -and (Test-ZipArchive -Path $final)) {
        "$(Get-Date -Format s) already complete: $final" | Add-Content -LiteralPath $log
        return
    }

    if ((Test-Path -LiteralPath $partial) -and (Test-ZipArchive -Path $partial)) {
        Move-Item -LiteralPath $partial -Destination $final -Force
        $item = Get-Item -LiteralPath $final
        $hash = (Get-FileHash -LiteralPath $final -Algorithm SHA256).Hash
        "$(Get-Date -Format s) completed existing partial bytes=$($item.Length) sha256=$hash" | Add-Content -LiteralPath $log
        return
    }

    "$(Get-Date -Format s) starting/resuming $Url" | Add-Content -LiteralPath $log
    $curlLog = Join-Path $OutputDirectory "curl_stderr.log"
    $curlArgs = @(
        "--fail", "--location", "--continue-at", "-", "--retry", "$RetryCount", "--retry-all-errors",
        "--retry-delay", "5", "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "60",
        "--output", $partial, $Url
    )
    $curlProcess = Start-Process -FilePath "curl.exe" -ArgumentList $curlArgs -Wait -PassThru -NoNewWindow `
        -RedirectStandardError $curlLog
    if (Test-Path -LiteralPath $curlLog) {
        Get-Content -LiteralPath $curlLog | Add-Content -LiteralPath $log
    }
    if ($curlProcess.ExitCode -ne 0) {
        "$(Get-Date -Format s) curl exit code $($curlProcess.ExitCode)" | Add-Content -LiteralPath $log
        throw "Download failed with curl exit code $($curlProcess.ExitCode): $Url"
    }

    if (-not (Test-ZipArchive -Path $partial)) {
        "$(Get-Date -Format s) ZIP integrity check failed: $partial" | Add-Content -LiteralPath $log
        throw "Downloaded file is not a readable ZIP archive: $partial"
    }

    Move-Item -LiteralPath $partial -Destination $final -Force
    $item = Get-Item -LiteralPath $final
    $hash = (Get-FileHash -LiteralPath $final -Algorithm SHA256).Hash
    "$(Get-Date -Format s) complete bytes=$($item.Length) sha256=$hash" | Add-Content -LiteralPath $log
}

Get-ResumableZip `
    -Url "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10567571/supplementaryFiles" `
    -OutputDirectory (Join-Path $DestinationRoot "eldjarn_2023") `
    -FileName "PMC10567571_SupplementaryFiles.zip"

Get-ResumableZip `
    -Url "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13190297/supplementaryFiles" `
    -OutputDirectory (Join-Path $DestinationRoot "lu_2026") `
    -FileName "PMC13190297_SupplementaryFiles.zip"
