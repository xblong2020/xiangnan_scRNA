$ErrorActionPreference = 'Stop'

$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
if (-not (Test-Path -LiteralPath $Rscript)) { throw "Rscript not found: $Rscript" }

$steps = @(
  'figure5_00_preflight_audit.R',
  'figure5_00b_frozen_signature_audit.R',
  'figure5_01_calculate_three_axis_scores.R',
  'figure5_02_orient_pseudotime.R',
  'figure5_03_build_patient_pseudobulk.R',
  'plot_figure5a_temporal_framework.R',
  'figure5_04_fit_three_axis_gam.R',
  'plot_figure5b_three_axis_gam.R',
  'figure5_05_prepare_temporal_heatmap.R',
  'plot_figure5c_temporal_heatmap.R',
  'figure5_06_calculate_temporal_landmarks.R',
  'figure5_07_bootstrap_temporal_landmarks.R',
  'plot_figure5d_temporal_landmarks.R',
  'figure5_08_analyze_precedence_probability.R',
  'plot_figure5e_precedence_matrix.R',
  'figure5_09_analyze_method_concordance.R',
  'plot_figure5f_method_concordance.R',
  'figure5_10_analyze_patient_temporal_order.R',
  'plot_figure5g_patient_forest.R',
  'plot_figure5h_overlapping_phase_model.R',
  'validate_figure5_temporal_positioning.R',
  'plot_figure5_temporal_positioning_preview.R',
  'generate_figure5_temporal_positioning_report.R'
)

Push-Location $ProjectRoot
try {
  foreach ($step in $steps) {
    Write-Host "[Figure 5] Running $step"
    $arguments = @((Join-Path $PSScriptRoot $step))
    if ($step -eq 'figure5_07_bootstrap_temporal_landmarks.R') { $arguments += @('--n-bootstrap', '1000', '--seed', '20260805') }
    & $Rscript @arguments
    if ($LASTEXITCODE -ne 0) { throw "$step failed with exit code $LASTEXITCODE" }
  }
} finally {
  Pop-Location
}

Write-Host '[Figure 5] Complete. See reports/figure5_temporal_positioning_report.md and metadata/driver/figure5_temporal_positioning/figure5_total_report.json.'
