param(
    [switch]$SkipDownload,
    [switch]$RefreshLiterature
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Project = Split-Path -Parent $Root
$Python = 'python'
$Rscript = 'Rscript'
$Inventory = Join-Path $Root 'data_processed\syn51365303_inventory.tsv'
$Selection = Join-Path $Root 'config\synapse_primary_panel_selection.tsv'
$ResourceRoot = if ($env:A1_RESOURCE_ROOT) { $env:A1_RESOURCE_ROOT } else { Join-Path $Project 'data\external' }
$DownloadDir = if ($env:UKB_PPP_DOWNLOAD_DIR) {
    $env:UKB_PPP_DOWNLOAD_DIR
} else {
    Join-Path $ResourceRoot 'UKB-PPP\syn51365303_European_discovery'
}

function Invoke-Step {
    param([string]$Label, [scriptblock]$Action)
    Write-Host "`n== $Label ==" -ForegroundColor Cyan
    & $Action
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE" }
}

Set-Location $Project

if ($RefreshLiterature) {
    Invoke-Step 'Literature search' { & $Python "$Root\scripts\01_search_literature.py" }
    Invoke-Step 'Title/abstract screening' { & $Python "$Root\scripts\02_title_abstract_screen.py" }
    Invoke-Step 'Open full-text retrieval' { & $Python "$Root\scripts\03_fetch_open_fulltext.py" }
    Invoke-Step 'Supplement inventory' { & $Python "$Root\scripts\04_inventory_supplement_tables.py" }
    Invoke-Step 'Provenance construction' { & $Python "$Root\scripts\05_build_literature_provenance.py" }
    Invoke-Step 'Screening/mapping exclusions' { & $Python "$Root\scripts\06_build_screening_mapping_exclusions.py" }
}

if (-not $SkipDownload) {
    if (-not $env:SYNAPSE_AUTH_TOKEN) {
        throw 'Set SYNAPSE_AUTH_TOKEN in the current environment. The token must not be written into this script.'
    }
    if (-not (Test-Path -LiteralPath $Inventory)) {
        Invoke-Step 'Synapse inventory' {
            & $Python "$Root\scripts\14_inventory_and_download_synapse_ukbppp.py" `
                --output $DownloadDir --inventory $Inventory --inventory-only
        }
    }
    Invoke-Step 'Select primary-panel Synapse files' {
        & $Python "$Root\scripts\15_select_primary_panel_synapse_files.py" `
            --provenance "$Root\tables\Table_Literature_Prioritized_Protein_Provenance.tsv" `
            --inventory $Inventory `
            --local-dir (Join-Path $ResourceRoot 'ukbppp_proteins') `
            --local-dir $DownloadDir `
            --output $Selection
    }
    Invoke-Step 'Resumable Synapse download' {
        & $Python "$Root\scripts\14_inventory_and_download_synapse_ukbppp.py" `
            --output $DownloadDir --inventory $Inventory --select $Selection `
            --reuse-inventory --retries 12 --retry-wait 15
    }
}

Invoke-Step 'GRCh37 protein coordinates' { & $Python "$Root\scripts\16_fetch_primary_gene_coordinates_grch37.py" }
Invoke-Step 'Protein-name-matched assay mapping' { & $Python "$Root\scripts\21_build_name_matched_mapping.py" }
Invoke-Step 'APOE alpha extraction' { & $Python "$Root\scripts\07_extract_apoe_alpha.py" }
Invoke-Step 'Genome-wide significant pQTL extraction' { & $Python "$Root\scripts\08_extract_reestimable_pqtl.py" }
Invoke-Step 'Main protein beta re-estimation' { & $Rscript "$Root\scripts\09_reestimate_literature_panel_beta.R" }
Invoke-Step 'APOE total effects' { & $Python "$Root\scripts\10_extract_apoe_total_effects.py" }
Invoke-Step 'Cis-only beta sensitivity' { & $Rscript "$Root\scripts\13_cis_only_beta_sensitivity.R" }
Invoke-Step 'Main two-step mediation' { & $Rscript "$Root\scripts\11_two_step_mediation.R" }
Invoke-Step 'Cis-only two-step mediation' { & $Rscript "$Root\scripts\11_two_step_mediation.R" --analysis-set=cis_only_sensitivity }
Invoke-Step 'Strict-annotation two-step mediation' { & $Rscript "$Root\scripts\11_two_step_mediation.R" --analysis-set=strict_annotation_sensitivity }
Invoke-Step 'Finalize provenance and comparisons' { & $Python "$Root\scripts\12_finalize_panel_tables.py" }
Invoke-Step 'Name-matched panel QA' { & $Python "$Root\scripts\22_validate_name_matched_panel.py" }
Invoke-Step 'Materialize Supplementary Tables S18-S29' { & $Python "$Root\scripts\23_package_supplementary_tables.py" }
Invoke-Step 'Integrate four verified priority studies' { & $Python "$Root\scripts\28_integrate_priority_studies_and_status.py" }
Invoke-Step 'Prepare quantitative figure source data' { & $Rscript "$Root\scripts\17_prepare_independent_protein_figure_data.R" }

Write-Host "`nA1 protein upgrade completed." -ForegroundColor Green
