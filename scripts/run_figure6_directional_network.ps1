$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$Rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
$MetadataDir = Join-Path $ProjectRoot 'metadata\driver\figure6_directional_network'
$DeltaExport = Join-Path $MetadataDir 'figure6_celloracle_programme_deltas_by_cell.tsv.gz'

if (-not (Test-Path $Rscript)) { throw "Rscript not found: $Rscript" }
& $Rscript (Join-Path $PSScriptRoot 'figure6_00_preflight_audit.R')
if (-not (Test-Path $DeltaExport)) {
  $WslProject = '/mnt/c/Users/Administrator/OneDrive/文档/湘南学院单细胞'
  wsl.exe -e bash -lc "source ~/miniconda3/etc/profile.d/conda.sh && conda activate celloracle_pip && cd '$WslProject' && export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8 NUMBA_NUM_THREADS=8 && python scripts/figure6_00_export_celloracle_programme_deltas.py"
  if ($LASTEXITCODE -ne 0) { throw "CellOracle Figure 6 export failed" }
}

$Scripts = @(
  'figure6_01_build_perturbation_matrix.R', 'plot_figure6a_perturbation_response_matrix.R',
  'figure6_02_prepare_representative_vector_fields.R', 'plot_figure6b_representative_vector_fields.R',
  'figure6_03_calculate_directional_asymmetry.R', 'plot_figure6c_directional_asymmetry.R',
  'figure6_04_compare_virtual_perturbation_methods.R', 'plot_figure6d_cross_method_concordance.R',
  'figure6_05_analyze_target_pathway_overlap.R', 'plot_figure6e_target_pathway_overlap.R',
  'figure6_06_fit_competing_network_models.R', 'plot_figure6f_competing_models.R',
  'figure6_09_negative_controls_and_confounders.R',
  'figure6_07_build_evidence_graded_network.R', 'plot_figure6g_evidence_graded_network.R',
  'figure6_08_prepare_literature_comparison.R', 'plot_figure6h_foxm1_cebpb_comparison.R',
  'figure6_10_finalize_report_and_preview.R',
  'validate_figure6_directional_network.R'
)
foreach ($Script in $Scripts) {
  Write-Host "[Figure 6] $Script"
  & $Rscript (Join-Path $PSScriptRoot $Script)
  if ($LASTEXITCODE -ne 0) { throw "Figure 6 script failed: $Script" }
}

Write-Host "Figure 6 complete."
Write-Host "Scripts: $PSScriptRoot\figure6_*.R; $PSScriptRoot\plot_figure6*.R"
Write-Host "Figures: $ProjectRoot\figures\driver\figure6*"
Write-Host "Source data: $MetadataDir"
Write-Host "Preview: $ProjectRoot\figures\driver\figure6_directional_network_preview"
Write-Host "Report: $ProjectRoot\reports\figure6_directional_network_report.md"
