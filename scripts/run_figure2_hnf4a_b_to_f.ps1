[CmdletBinding()]
param(
    [switch]$ForceScTenifold
)

$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
$OutputEncoding = [Console]::OutputEncoding
$ProjectRoot = Split-Path -Parent $PSScriptRoot
$Python = 'python'
$Rscript = 'C:\Program Files\R\R-4.5.0\bin\Rscript.exe'
$TargetTf = 'HNF4A'
Set-Location -LiteralPath $ProjectRoot

if (-not (Test-Path -LiteralPath $Rscript)) {
    throw "Required Rscript not found: $Rscript"
}

function Invoke-Step {
    param([string]$Name, [scriptblock]$Command)
    Write-Host "[$Name]"
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Step failed ($LASTEXITCODE): $Name"
    }
}

Invoke-Step '1 preflight' { & $Python scripts/preflight_figure2_hnf4a.py --target-tf $TargetTf }
Invoke-Step '2 Figure 2B data' { & $Python scripts/prepare_figure2b_hnf4a_data.py }
Invoke-Step '3 Figure 2B plot' { & $Rscript scripts/plot_figure2b_hnf4a_baseline.R --target-tf $TargetTf }
Invoke-Step '4 HNF4A CellOracle data' { & $Python scripts/prepare_figure2c_hnf4a_data.py --target-tf $TargetTf }
Invoke-Step '5 Figure 2C perturbation' { & $Rscript scripts/plot_figure2c_hnf4a_perturbation.R --target-tf $TargetTf }
Invoke-Step '6 Figure 2C inner product' { & $Rscript scripts/plot_figure2c_hnf4a_inner_product.R --target-tf $TargetTf }
Invoke-Step '7 Figure 2D pseudotime' { & $Rscript scripts/plot_figure2d_hnf4a_pseudotime_inner_product.R --target-tf $TargetTf }
Invoke-Step '8a identity-high input' { & $Python scripts/prepare_figure2e_hnf4a_sctenifoldknk.py --target-tf $TargetTf }

$ScTenifoldResult = 'metadata\driver\figure2e_hnf4a_sctenifoldknk\figure2e_hnf4a_normal_reference_perturbation_genes.tsv'
if ($ForceScTenifold -or -not (Test-Path -LiteralPath $ScTenifoldResult)) {
    Invoke-Step '8b identity-high scTenifoldKnk' {
        & $Rscript scripts/run_figure2e_hnf4a_sctenifoldknk.R --target-tf $TargetTf
    }
} else {
    Write-Host '[8b identity-high scTenifoldKnk] Existing HNF4A-specific result retained; use -ForceScTenifold to rerun.'
}

Invoke-Step '9 Figure 2E' { & $Rscript scripts/plot_figure2e_hnf4a_sctenifoldknk.R --target-tf $TargetTf }
Invoke-Step '10 Figure 2F' { & $Rscript scripts/plot_figure2f_hnf4a_pathway_enrichment.R --target-tf $TargetTf }
Invoke-Step '11 validation' { & $Python scripts/validate_figure2_hnf4a_b_to_f.py }
Invoke-Step '12 review montage' { & $Python scripts/make_figure2_hnf4a_montage.py }

$Validation = Get-Content -Raw -Encoding UTF8 -LiteralPath 'metadata\driver\figure2_hnf4a_b_to_f_validation_report.json' | ConvertFrom-Json
$Preflight = Get-Content -Raw -Encoding UTF8 -LiteralPath 'metadata\driver\figure2_hnf4a_preflight_report.json' | ConvertFrom-Json
$D = Get-Content -Raw -Encoding UTF8 -LiteralPath 'metadata\driver\figure2d_hnf4a\figure2d_hnf4a_report.json' | ConvertFrom-Json
$E = Get-Content -Raw -Encoding UTF8 -LiteralPath 'metadata\driver\figure2e_hnf4a\figure2e_hnf4a_report.json' | ConvertFrom-Json
$F = Get-Content -Raw -Encoding UTF8 -LiteralPath 'metadata\driver\figure2f_hnf4a\figure2f_hnf4a_report.json' | ConvertFrom-Json

Write-Host 'Scripts:'
Get-ChildItem -LiteralPath scripts -File | Where-Object { $_.Name -match 'figure2.*hnf4a|hnf4a.*figure2' } |
    Select-Object -ExpandProperty FullName
Write-Host 'Main figures:'
Get-ChildItem -Path figures\driver -Recurse -File |
    Where-Object { $_.Name -match '^figure2[b-f]_hnf4a_.*\.(pdf|png|svg|tiff)$' } |
    Select-Object -ExpandProperty FullName
Write-Host "Main statistics: UMAP rho=$($D.umap.spearman_rho); Figure 2E significant genes=$($E.n_significant_excluding_target); Figure 2F significant pathways=$($F.n_significant_pathways)"
Write-Host "Warnings/review-risk flags: $($Validation.review_risk_flags -join ' | ')"
Write-Host "SOX4 original files unmodified: $($Validation.all_required_checks_passed); hashes checked=$($Validation.n_sox4_hashes_checked)"
