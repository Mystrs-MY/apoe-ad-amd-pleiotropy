$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$UpgradeRoot = $PSScriptRoot
$ProjectRoot = (Resolve-Path (Join-Path $UpgradeRoot '..')).Path
$LogRoot = Join-Path $UpgradeRoot 'logs'
$RunStamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$RunLog = Join-Path $LogRoot "submission_core_QA_$RunStamp.log"
New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null

$failures = [System.Collections.Generic.List[string]]::new()

function Require-File {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        $failures.Add("Missing required file: $Path")
    }
}

function Reject-NamePattern {
    param([string]$Root, [string]$Pattern)
    if (-not (Test-Path -LiteralPath $Root)) { return }
    Get-ChildItem -LiteralPath $Root -File -Recurse | Where-Object { $_.Name -match $Pattern } | ForEach-Object {
        $failures.Add("Forbidden submission asset: $($_.FullName)")
    }
}

$mainDir = Join-Path $ProjectRoot 'figures_submission\main'
$suppFigDir = Join-Path $ProjectRoot 'figures_submission\supplementary_figures'
$suppTableDir = Join-Path $ProjectRoot 'tables_submission\supplementary_tables'
$mainTableDir = Join-Path $ProjectRoot 'tables_submission\main_tables'

# Submission artwork is maintained outside the minimal public code package.

foreach ($name in @('Table1_APOE_Isoform_Defining_Variant_Effects.csv', 'Table2_APOE_Exclusion_MR.csv')) {
    Require-File (Join-Path $mainTableDir $name)
}

Reject-NamePattern $mainTableDir '^Table1_APOE_Wald_Ratio\.csv$'

foreach ($name in @(
    'TableS1a_Wald_Ratio_Full.csv',
    'TableS1b_Wald_Ratio_Display.csv',
    'TableS2a_GW_MR_Full.csv',
    'TableS2b_GW_MR_With_vs_Without_APOE.csv',
    'TableS2c_APOE_Exclusion_Window_Sensitivity.tsv',
    'TableS3a_MiXeR_Univariate.csv',
    'TableS3b_MiXeR_Bivariate.csv',
    'TableS3c_MiXeR_Venn.csv',
    'TableS3d_HDL_Results.csv',
    'TableS4_Protein_Selection.csv',
    'TableS5a_LAVA_AD_vs_DryAMD.csv',
    'TableS5b_LAVA_AD_vs_WetAMD.csv',
    'TableS5c_LAVA_AD_vs_AnyAMD.csv',
    'TableS5d_HyPrColoc_APOE_Result.csv',
    'TableS5e_MAGMA_Locus_Summary.csv',
    'TableS5f_HyPrColoc_Prior_Sensitivity.tsv',
    'TableS5g_MAGMA_Full_Gene_Results.tsv.gz',
    'TableS6_rs429358_Protein_Effects.csv',
    'TableS7_Protein_MR_Mediation.csv',
    'TableS8_Pathway_Mediation_Summary.csv',
    'TableS9_Pathway_MVMR.csv',
    'TableS10_Ridge_MVMR_Detailed.csv',
    'TableS11_Covariance_Mapping_Bootstrap_Sensitivity.tsv',
    'TableS11b_Shared_Protein_Instrument_Audit.tsv',
    'TableS11c_Leave_One_Protein_Out_Aggregate.tsv',
    'TableS12_Study_Characteristics.csv',
    'TableS13_Subgroup_Analyses.csv',
    'TableS14_Sensitivity_Analyses.csv',
    'TableS15_GRADE_SoF.csv',
    'TableS16_Excluded_Studies.csv',
    'TableS17_GRADE_Evidence_Profile.csv',
    'TableS18_Literature_Study_Screening.tsv',
    'TableS19_Protein_Provenance_Master.tsv',
    'TableS20_Cross_Platform_Mapping.tsv',
    'TableS21_Protein_to_UKB_PPP_Assay_Mapping.tsv',
    'TableS22_APOE_Variant_to_Protein_Alpha.tsv',
    'TableS23_Same_Assay_Protein_Beta.tsv',
    'TableS24_APOE_Linkable_Eligibility_Flow.tsv',
    'TableS25_Expanded_Primary_Two_Step_Mediation.tsv',
    'TableS26_Cis_Only_Two_Step_Mediation.tsv',
    'TableS27_Expanded_Strict_Legacy_Comparison.tsv',
    'TableS28_Excluded_Name_Matched_Panel_Genes.tsv',
    'TableS29_Figure5_Source_Data_and_QA_Manifest.tsv',
    'TableS30_Five_Protein_PWAS_Data_Sources_and_Crosswalk.tsv',
    'TableS31_Five_Protein_PWAS_APOE_Alpha_and_Beta.tsv',
    'TableS32_Five_Protein_PWAS_Mediation_and_Aggregate_Sensitivity.tsv',
    'TableS33a_Kunkle_AD_Instrument_Coverage.tsv',
    'TableS33b_Kunkle_AD_Protein_Beta.tsv',
    'TableS33c_Kunkle_AD_Two_Step_Mediation.tsv',
    'TableS33d_Kunkle_vs_Wightman_Aggregate_Comparison.tsv',
    'TableS34a_PAV_Epitope_Instrument_Audit.tsv',
    'TableS34b_PAV_Filtered_Protein_Beta.tsv',
    'TableS34c_PAV_Filtered_Two_Step_Mediation.tsv',
    'TableS34d_PAV_Filtered_Aggregate_Comparison.tsv',
    'TableS35a_Lu2026_Candidate_Feasibility.tsv',
    'TableS35b_Lu2026_Exact_Assay_Gate.tsv',
    'TableS35c_Lu2026_Exact_Assay_APOE_Alpha.tsv',
    'TableS36a_Lu2026_Exact_Assay_Protein_Beta.tsv',
    'TableS36b_Lu2026_Exact_Assay_Two_Step_Mediation.tsv',
    'TableS36c_Lu2026_Exact_Assay_Run_Summary.tsv',
    'TableS37a_deCODE_Exact_Assay_Gate.tsv',
    'TableS37b_deCODE_APOE_Alpha.tsv',
    'TableS38a_deCODE_Genome_Wide_Two_Step_Mediation.tsv',
    'TableS38b_deCODE_Cis_Only_Two_Step_Mediation.tsv',
    'TableS38c_deCODE_Cis_PAV_Filtered_Two_Step_Mediation.tsv',
    'TableS39a_deCODE_Shared_Instrument_Audit.tsv',
    'TableS39b_deCODE_PAV_Epitope_Audit.tsv',
    'TableS39c_deCODE_Raw_Normalization_Sensitivity.tsv',
    'TableS39d_deCODE_File_Integrity_and_Validation.tsv'
)) {
    Require-File (Join-Path $suppTableDir $name)
}

Reject-NamePattern $mainDir '^Figure_6\.'
Reject-NamePattern $suppFigDir '^FigS3[a-d]_'
Reject-NamePattern $suppFigDir '^FigS11[a-d]_'
Reject-NamePattern $suppFigDir 'C3_full|C3_SuSiE|C3_reciprocal|C3_evidence'
Reject-NamePattern $suppTableDir 'Drug_Target|Drug-Target|C3_full|C3_coloc|C3_conditional|^TableS([4-9][0-9]|[1-9][0-9]{2,})'

$manuscriptDir = Join-Path $ProjectRoot 'manuscripts'
$chineseManuscript = Get-ChildItem -LiteralPath $manuscriptDir -File -Filter '*V1.4.md' |
    Where-Object { $_.Name -notlike '*CellGenomics*' } |
    Sort-Object Length -Descending |
    Select-Object -First 1 -ExpandProperty FullName
$manuscripts = @(
    $chineseManuscript,
    (Join-Path $ProjectRoot 'manuscripts\Article1_CellGenomics_English_Manuscript_V1.4.md'),
    (Join-Path $ProjectRoot 'manuscripts\Article1_CellGenomics_Supplementary_Information_V1.4.md'),
    (Join-Path $ProjectRoot 'manuscripts\submission_materials\Article1_CellGenomics_Cover_Letter_V1.4.md'),
    (Join-Path $ProjectRoot 'manuscripts\submission_materials\Article1_CellPress_Editorial_Elements_V1.4.md'),
    (Join-Path $ProjectRoot 'manuscripts\submission_materials\Article1_CellGenomics_Highlights_and_In_Brief_V1.4.md'),
    (Join-Path $ProjectRoot 'manuscripts\submission_materials\Article1_CellGenomics_Key_Resources_Table_V1.4.tsv')
)
$forbiddenText = 'Exploratory target MR|exploratory target prioritization|exploratory target-MR|drug-target MR|drug target MR|药物靶点|探索性靶点|转化优先级|转化预警|C3 full-region|C3 完整位点|OID30776|rs11569470|rs2230199|LINCS|CMap|FAERS|VigiBase|Figure 6(?![a-zA-Z])|Fig\. 6(?![a-zA-Z])|__NEW_'
foreach ($manuscript in $manuscripts) {
    Require-File $manuscript
    if (Test-Path -LiteralPath $manuscript) {
        $matches = Select-String -LiteralPath $manuscript -Pattern $forbiddenText -AllMatches
        foreach ($match in $matches) {
            $failures.Add("Forbidden manuscript text: ${manuscript}:$($match.LineNumber): $($match.Line.Trim())")
        }
    }
}

$summary = @(
    "A1 submission-facing core QA",
    "Timestamp: $(Get-Date -Format o)",
    "Project: $ProjectRoot",
    "Failures: $($failures.Count)"
)
$summary | Set-Content -LiteralPath $RunLog -Encoding UTF8
$failures | Add-Content -LiteralPath $RunLog -Encoding UTF8

if ($failures.Count -gt 0) {
    $failures | ForEach-Object { Write-Host $_ -ForegroundColor Red }
    throw "Submission QA failed with $($failures.Count) issue(s). See $RunLog"
}

Write-Host "Submission QA passed. Log: $RunLog" -ForegroundColor Green
