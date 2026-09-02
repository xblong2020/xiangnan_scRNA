$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$RscriptExe = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
$Figure7RLib = Join-Path $ProjectRoot 'data\processed\driver\figure7_external_validation\rlib'

if (-not (Test-Path -LiteralPath $RscriptExe)) {
    throw "Rscript 4.5.0 not found: $RscriptExe"
}

New-Item -ItemType Directory -Force -Path $Figure7RLib | Out-Null
Set-Location -LiteralPath $ProjectRoot

function Invoke-Figure7Step {
    param([string]$Label, [string]$Script)
    $started = Get-Date
    Write-Host "[Figure7] START $Label"
    & $RscriptExe $Script
    if ($LASTEXITCODE -ne 0) {
        throw "Figure 7 step failed ($LASTEXITCODE): $Label"
    }
    $elapsed = (Get-Date) - $started
    Write-Host ("[Figure7] DONE  {0} ({1:n1} min)" -f $Label, $elapsed.TotalMinutes)
}

Invoke-Figure7Step '01 preflight and frozen-signature audit' 'scripts\figure7_00_preflight_audit.R'
Invoke-Figure7Step '02 cohort-specific expression preparation and coverage audit' 'scripts\figure7_01_prepare_bulk_expression.R'
Invoke-Figure7Step '03 rank-based axis scores, controls and 1000x matched random signatures' 'scripts\figure7_02_calculate_bulk_axis_scores.R'
Invoke-Figure7Step '04 Figure 7A cohort flow' 'scripts\plot_figure7a_cohort_flow.R'
Invoke-Figure7Step '05 Figure 7B frozen signature mapping' 'scripts\plot_figure7b_bulk_signature_mapping.R'
Invoke-Figure7Step '06 tumour-normal effects and meta-analysis' 'scripts\figure7_03_analyze_tumour_normal.R'
Invoke-Figure7Step '07 Figure 7C forest plot' 'scripts\plot_figure7c_tumour_normal_forest.R'
Invoke-Figure7Step '08 clinicopathological association models' 'scripts\figure7_04_analyze_clinical_associations.R'
Invoke-Figure7Step '09 Figure 7D clinical heatmap' 'scripts\plot_figure7d_clinical_heatmap.R'
Invoke-Figure7Step '10 multivariable Cox, EPV, PH and collinearity diagnostics' 'scripts\figure7_05_fit_multivariable_cox.R'
Invoke-Figure7Step '11 Figure 7E Cox forest plot' 'scripts\plot_figure7e_multivariable_cox_forest.R'
Invoke-Figure7Step '12 internal and locked external prediction validation' 'scripts\figure7_06_evaluate_incremental_prediction.R'
Invoke-Figure7Step '13 Figure 7F incremental prediction' 'scripts\plot_figure7f_incremental_prediction.R'
Invoke-Figure7Step '14 prespecified locked survival grouping' 'scripts\figure7_07_prepare_survival_groups.R'
Invoke-Figure7Step '15 Figure 7G survival curves' 'scripts\plot_figure7g_survival_curves.R'
Invoke-Figure7Step '16 sensitivity and specificity analyses' 'scripts\figure7_08_run_sensitivity_analyses.R'
Invoke-Figure7Step '17 Figure 7H sensitivity summary' 'scripts\plot_figure7h_sensitivity_summary.R'
Invoke-Figure7Step '18 automatic validation' 'scripts\validate_figure7_external_validation.R'
Invoke-Figure7Step '19 A-to-H review preview' 'scripts\plot_figure7_external_validation_preview.R'
Invoke-Figure7Step '20 final report and terminal summary' 'scripts\generate_figure7_external_validation_report.R'
