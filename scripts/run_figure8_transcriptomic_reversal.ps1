$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
$Python = Join-Path $ProjectRoot '.venv-drugreflector\Scripts\python.exe'
$MetadataDir = Join-Path $ProjectRoot 'metadata\driver\figure8_transcriptomic_reversal'
$VariantInput = Join-Path $ProjectRoot 'data\processed\driver\figure8_transcriptomic_reversal\figure8_signature_variants_wide.tsv.gz'

if (-not (Test-Path -LiteralPath $Rscript)) { throw "Rscript not found: $Rscript" }
if (-not (Test-Path -LiteralPath $Python)) { throw "DrugReflector Python not found: $Python" }

function Invoke-Figure8R {
    param([Parameter(Mandatory = $true)][string]$ScriptName)
    $Path = Join-Path $PSScriptRoot $ScriptName
    Write-Host "[Figure8] R: $ScriptName"
    & $Rscript $Path
    if ($LASTEXITCODE -ne 0) { throw "Figure 8 R step failed: $ScriptName" }
}

Set-Location -LiteralPath $ProjectRoot
Invoke-Figure8R 'figure8_00_preflight_audit.R'

Write-Host '[Figure8] DrugReflector: export frozen model gene order'
& $Python (Join-Path $PSScriptRoot 'figure8_drugreflector_inference.py') --mode export --metadata-dir $MetadataDir --seed 20260805
if ($LASTEXITCODE -ne 0) { throw 'DrugReflector model-gene export failed' }

Invoke-Figure8R 'figure8_01_prepare_signature_composition.R'
Invoke-Figure8R 'plot_figure8a_target_state_definition.R'
Invoke-Figure8R 'plot_figure8b_signature_composition.R'
Invoke-Figure8R 'plot_figure8c_reversal_workflow.R'

Write-Host '[Figure8] DrugReflector: preregistered signature-variant inference'
& $Python (Join-Path $PSScriptRoot 'figure8_drugreflector_inference.py') --mode variants --input $VariantInput --metadata-dir $MetadataDir --seed 20260805
if ($LASTEXITCODE -ne 0) { throw 'DrugReflector variant inference failed' }

Invoke-Figure8R 'figure8_02_analyze_drugreflector_stability.R'
Invoke-Figure8R 'plot_figure8d_drugreflector_stability.R'
Invoke-Figure8R 'figure8_03_analyze_cross_method_concordance.R'
Invoke-Figure8R 'plot_figure8e_cross_method_concordance.R'
Invoke-Figure8R 'figure8_04_analyze_mechanism_classes.R'
Invoke-Figure8R 'plot_figure8f_mechanism_classes.R'
Invoke-Figure8R 'figure8_05_validate_external_perturbation_signatures.R'
Invoke-Figure8R 'plot_figure8g_external_signature_validation.R'
Invoke-Figure8R 'figure8_07_run_random_signature_benchmark.R'
Invoke-Figure8R 'figure8_08_analyze_toxicity_stress_penalties.R'
Invoke-Figure8R 'figure8_06_build_integrated_reversal_score.R'
Invoke-Figure8R 'plot_figure8h_integrated_prioritization.R'
Invoke-Figure8R 'validate_figure8_transcriptomic_reversal.R'
Invoke-Figure8R 'figure8_09_build_preview_report.R'
Invoke-Figure8R 'validate_figure8_transcriptomic_reversal.R'

Write-Host '[Figure8] Completed. See reports\figure8_transcriptomic_reversal_report.md'
