$ErrorActionPreference = 'Stop'
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
$Python = Join-Path $ProjectRoot '.venv-drugreflector\Scripts\python.exe'
$Meta = Join-Path $ProjectRoot 'metadata\driver\figure8_transcriptomic_reversal_v2_mainfigure'
$Data = Join-Path $ProjectRoot 'data\processed\driver\figure8_transcriptomic_reversal_v2_mainfigure'
$env:LC_ALL = 'Chinese'
Set-Location -LiteralPath $ProjectRoot

function Invoke-Figure8V2R {
    param([Parameter(Mandatory = $true)][string]$Name, [string[]]$Arguments = @())
    & $Rscript (Join-Path $PSScriptRoot $Name) @Arguments
    if ($LASTEXITCODE -ne 0) { throw "Figure 8 v2 R step failed: $Name" }
}

$Baseline = Join-Path $Meta 'figure8_v2_protected_figure1_7_hash_before.tsv'
if (-not (Test-Path -LiteralPath $Baseline)) {
    Invoke-Figure8V2R 'figure8_v2_00_freeze_audit.R' @('--workers=4')
}
Invoke-Figure8V2R 'figure8_v2_01_build_continuous_signature.R'

& $Python (Join-Path $PSScriptRoot 'figure8_v2_drugreflector_inference.py') --mode variants --input (Join-Path $Data 'figure8_v2_signature_variants_wide.tsv.gz') --metadata-dir $Meta --seed 20260805
if ($LASTEXITCODE -ne 0) { throw 'Figure 8 v2 DrugReflector variant inference failed' }
Invoke-Figure8V2R 'figure8_v2_02_analyze_drugreflector.R'
Invoke-Figure8V2R 'figure8_v2_03_matched_random.R' @('--stage=prepare')
& $Python (Join-Path $PSScriptRoot 'figure8_v2_drugreflector_inference.py') --mode random --input (Join-Path $Data 'figure8_v2_matched_random_signatures_wide.tsv.gz') --watchlist (Join-Path $Meta 'figure8_v2_candidate_watchlist.tsv') --metadata-dir $Meta --top-n 200 --batch-size 25 --seed 20260805
if ($LASTEXITCODE -ne 0) { throw 'Figure 8 v2 matched-null inference failed' }
Invoke-Figure8V2R 'figure8_v2_03_matched_random.R' @('--stage=summarize')
Invoke-Figure8V2R 'figure8_v2_04_fetch_external_resources.R'
Invoke-Figure8V2R 'figure8_v2_05_cross_framework.R'
Invoke-Figure8V2R 'figure8_v2_06_moa_network.R'
Invoke-Figure8V2R 'figure8_v2_07_prism.R'
Invoke-Figure8V2R 'figure8_v2_08_nuisance_literature.R'
Invoke-Figure8V2R 'figure8_v2_09_integrate_evidence.R'
Invoke-Figure8V2R 'figure8_v2_10_plot_mainfigure.R'
Invoke-Figure8V2R 'figure8_v2_11_plot_extended_data.R'
Invoke-Figure8V2R 'figure8_v2_12_validate_report.R'

