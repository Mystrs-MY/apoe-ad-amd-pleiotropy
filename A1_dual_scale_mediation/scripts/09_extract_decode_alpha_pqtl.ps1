param(
    [ValidateSet('smp', 'raw')]
    [string]$Normalization = 'smp',
    [string[]]$Assays
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$inventoryPath = Join-Path $root 'tables\decode_target_object_inventory.tsv'
$inputDir = Join-Path $root "data_raw\decode_somascan\$Normalization"
$outputDir = Join-Path $root 'data_processed'
$logDir = Join-Path $root 'logs'
$awkScript = Join-Path $PSScriptRoot '09_extract_decode_rows.awk'
$gzip = 'C:\rtools43\usr\bin\gzip.exe'
$gawk = 'C:\rtools43\usr\bin\gawk.exe'

foreach ($required in @($inventoryPath, $awkScript, $gzip, $gawk)) {
    if (-not (Test-Path -LiteralPath $required)) { throw "Missing required input: $required" }
}

$rows = Import-Csv -Delimiter "`t" -LiteralPath $inventoryPath |
    Where-Object { $_.normalization -eq $Normalization }
if ($Assays) {
    $rows = @($rows | Where-Object { $Assays -contains $_.SomaScan_SeqId })
    $missing = @($Assays | Where-Object { $_ -notin $rows.SomaScan_SeqId })
    if ($missing.Count) { throw "Requested assays absent from inventory: $($missing -join ',')" }
    $expectedCount = $Assays.Count
    $subsetSuffix = '_gated_' + (($Assays | Sort-Object) -join '_')
} else {
    $expectedCount = 9
    $subsetSuffix = ''
}
if ($rows.Count -ne $expectedCount) { throw "Expected $expectedCount $Normalization inventory rows, found $($rows.Count)" }

New-Item -ItemType Directory -Force -Path $outputDir, $logDir | Out-Null
$alphaPath = Join-Path $outputDir "decode_${Normalization}${subsetSuffix}_direct_APOE_alpha.tsv"
$pqtlPath = Join-Path $outputDir "decode_${Normalization}${subsetSuffix}_pqtl_candidates_p5e8.tsv"
$logPath = Join-Path $logDir "decode_${Normalization}${subsetSuffix}_row_extraction.tsv"
$alphaTmp = "$alphaPath.tmp"
$pqtlTmp = "$pqtlPath.tmp"
$logTmp = "$logPath.tmp"

$sourceHeader = "Chrom`tPos`tName`trsids`teffectAllele`totherAllele`tBeta`tPval`tminus_log10_pval`tSE`tN`tImpMAF"
[System.IO.File]::WriteAllText(
    $alphaTmp,
    "normalization`tanalysis_role`tgene_symbol`tSomaScan_SeqId`tUniProt_ID`tsource_file`trequested_variant`t$sourceHeader`n",
    [System.Text.UTF8Encoding]::new($false)
)
[System.IO.File]::WriteAllText(
    $pqtlTmp,
    "normalization`tanalysis_role`tgene_symbol`tSomaScan_SeqId`tUniProt_ID`tsource_file`tpasses_decode_studywide_p1_8e9`t$sourceHeader`n",
    [System.Text.UTF8Encoding]::new($false)
)
[System.IO.File]::WriteAllText(
    $logTmp,
    "normalization`tgene_symbol`tSomaScan_SeqId`tdirect_alpha_rows`tpqtl_candidate_rows`n",
    [System.Text.UTF8Encoding]::new($false)
)

foreach ($row in $rows) {
    $fileName = Split-Path -Leaf $row.object_key
    $source = Join-Path $inputDir $fileName
    if (-not (Test-Path -LiteralPath $source)) { throw "Missing completed source file: $source" }
    if ((Get-Item -LiteralPath $source).Length -ne [int64]$row.content_length_bytes) {
        throw "Size mismatch for $source"
    }
    & $gzip -t -- $source
    if ($LASTEXITCODE -ne 0) { throw "gzip integrity test failed: $source" }

    $role = $row.analysis_role
    $sourceArg = $source.Replace('\', '/')
    $alphaArg = $alphaTmp.Replace('\', '/')
    $pqtlArg = $pqtlTmp.Replace('\', '/')
    $result = & $gzip -dc -- $source | & $gawk `
        -v "normalization=$Normalization" `
        -v "analysis_role=$role" `
        -v "gene=$($row.gene_symbol)" `
        -v "seqid=$($row.SomaScan_SeqId)" `
        -v "uniprot=$($row.UniProt_ID)" `
        -v "source_file=$sourceArg" `
        -v "alpha_out=$alphaArg" `
        -v "pqtl_out=$pqtlArg" `
        -v p_primary=5e-8 `
        -v p_strict=1.8e-9 `
        -f $awkScript
    if ($LASTEXITCODE -ne 0) { throw "Row extraction failed: $source" }
    [System.IO.File]::AppendAllText($logTmp, "$result`n", [System.Text.UTF8Encoding]::new($false))
    Write-Output $result
}

Move-Item -Force -LiteralPath $alphaTmp -Destination $alphaPath
Move-Item -Force -LiteralPath $pqtlTmp -Destination $pqtlPath
Move-Item -Force -LiteralPath $logTmp -Destination $logPath
Write-Output "Alpha output: $alphaPath"
Write-Output "pQTL candidates: $pqtlPath"
Write-Output "Extraction log: $logPath"
