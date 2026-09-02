[CmdletBinding()]
param(
    [string]$PythonExe = "",
    [string]$ValidationPythonExe = "python",
    [string]$RscriptExe = "C:\Program Files\R\R-4.5.0\bin\Rscript.exe",
    [switch]$ForceNetworks
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $ProjectRoot

if (-not $PythonExe) {
    $ProjectPython = Join-Path $ProjectRoot ".venv-scvi\Scripts\python.exe"
    $PythonExe = if (Test-Path -LiteralPath $ProjectPython) { $ProjectPython } else { "python" }
}
if (-not (Test-Path -LiteralPath $RscriptExe)) {
    throw "Required project Rscript was not found: $RscriptExe"
}

function Invoke-Step {
    param([string]$Name, [scriptblock]$Command)
    Write-Host "`n[$Name]" -ForegroundColor Cyan
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Name failed with exit code $LASTEXITCODE"
    }
}

Invoke-Step "1/16 preflight" { & $PythonExe scripts\preflight_figure3_egr1.py }
Invoke-Step "2/16 Figure 3A evidence" { & $PythonExe scripts\prepare_figure3a_stress_transition_selection.py }
Invoke-Step "3/16 Figure 3A plot" { & $RscriptExe scripts\plot_figure3a_stress_transition_selection.R }
Invoke-Step "4/16 Figure 3B common baseline" { & $PythonExe scripts\prepare_figure3b_egr1_data.py }
Invoke-Step "5/16 Figure 3B plot" { & $RscriptExe scripts\plot_figure3b_egr1_baseline.R }
Invoke-Step "6/16 EGR1 CellOracle data" { & $PythonExe scripts\prepare_figure3c_egr1_data.py }
Invoke-Step "7/16 Figure 3C perturbation" { & $RscriptExe scripts\plot_figure3c_egr1_perturbation.R }
Invoke-Step "8/16 Figure 3C inner product" { & $RscriptExe scripts\plot_figure3c_egr1_inner_product.R }
Invoke-Step "9/16 Figure 3D statistics" { & $PythonExe scripts\prepare_figure3d_egr1_statistics.py }
Invoke-Step "9b/16 Figure 3D plot" { & $RscriptExe scripts\plot_figure3d_egr1_pseudotime_inner_product.R }
Invoke-Step "10/16 stress-transition inputs" { & $PythonExe scripts\prepare_figure3e_egr1_sctenifoldknk.py }

$MainConsensus = "metadata\driver\figure3e_egr1\figure3e_egr1_stressed_regenerative_consensus_perturbation_genes.tsv"
if ($ForceNetworks -or -not (Test-Path -LiteralPath $MainConsensus)) {
    Invoke-Step "11/16 main EGR1 scTenifoldKnk" {
        & $RscriptExe scripts\run_figure3e_egr1_sctenifoldknk.R `
            --input_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative `
            --output_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/results `
            --metadata_dir=metadata/driver/figure3e_egr1 `
            --subset=stressed_regenerative `
            --seeds=15071990,15071991,15071992 `
            --nc_nnet=10 --nc_ncells=500 --ncores=8
    }
} else {
    Write-Host "`n[11/16 main EGR1 scTenifoldKnk] Reusing existing consensus: $MainConsensus" -ForegroundColor DarkCyan
}

$DeterminismReport = "metadata\driver\figure3e_egr1_determinism\figure3e_egr1_same_seed_determinism_report.json"
if ($ForceNetworks -or -not (Test-Path -LiteralPath $DeterminismReport)) {
    Invoke-Step "11a/16 same-seed determinism repeat" {
        & $RscriptExe scripts\run_figure3e_egr1_sctenifoldknk.R `
            --input_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative `
            --output_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/stressed_regenerative/determinism_repeat `
            --metadata_dir=metadata/driver/figure3e_egr1_determinism `
            --subset=stressed_regenerative --seeds=15071990 `
            --nc_nnet=10 --nc_ncells=500 --ncores=8
    }
    Invoke-Step "11a2/16 same-seed determinism audit" {
        & $PythonExe scripts\audit_figure3e_egr1_determinism.py
    }
} else {
    Write-Host "[11a/16 determinism] Reusing independent same-seed audit." -ForegroundColor DarkCyan
}

$SensitivitySpecs = @(
    @{ Subset = "stressed_injured"; Cells = 247 },
    @{ Subset = "intermediate_pseudotime"; Cells = 500 },
    @{ Subset = "malignant_like"; Cells = 500 }
)
foreach ($Spec in $SensitivitySpecs) {
    $Subset = $Spec.Subset
    $Consensus = "metadata\driver\figure3e_egr1_sensitivity\figure3e_egr1_${Subset}_consensus_perturbation_genes.tsv"
    if ($ForceNetworks -or -not (Test-Path -LiteralPath $Consensus)) {
        Invoke-Step "11b/16 sensitivity network: $Subset" {
            & $RscriptExe scripts\run_figure3e_egr1_sctenifoldknk.R `
                "--input_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/$Subset" `
                "--output_dir=data/processed/driver/figure3e_egr1_sctenifoldknk/$Subset/results" `
                --metadata_dir=metadata/driver/figure3e_egr1_sensitivity `
                "--subset=$Subset" --seeds=15071990 --nc_nnet=3 `
                "--nc_ncells=$($Spec.Cells)" --ncores=8
        }
    } else {
        Write-Host "[11b/16 sensitivity network] Reusing $Subset consensus." -ForegroundColor DarkCyan
    }
}

Invoke-Step "12/16 Figure 3E plot" { & $RscriptExe scripts\plot_figure3e_egr1_sctenifoldknk.R }
Invoke-Step "12b/16 Figure 3E sensitivity summary" { & $PythonExe scripts\summarize_figure3e_egr1_sensitivity.py }
Invoke-Step "12c/16 Figure 3E sensitivity plot" { & $RscriptExe scripts\plot_figure3e_egr1_sensitivity.R }
Invoke-Step "13/16 Figure 3F strict ORA" { & $RscriptExe scripts\plot_figure3f_egr1_pathway_enrichment.R }
Invoke-Step "14/16 three-axis audit" { & $PythonExe scripts\audit_three_axis_figure_consistency.py }
Invoke-Step "15/16 review montage" { & $PythonExe scripts\make_figure3_egr1_montage.py }
Invoke-Step "15b/16 validation" { & $ValidationPythonExe scripts\validate_figure3_egr1.py }
Invoke-Step "16/16 final report" { & $PythonExe scripts\generate_figure3_egr1_report.py }

$DReport = Get-Content -LiteralPath "metadata\driver\figure3d_egr1\figure3d_egr1_report.json" -Raw | ConvertFrom-Json
$EReport = Get-Content -LiteralPath "metadata\driver\figure3e_egr1\figure3e_egr1_report.json" -Raw | ConvertFrom-Json
$FReport = Get-Content -LiteralPath "metadata\driver\figure3f_egr1\figure3f_egr1_report.json" -Raw | ConvertFrom-Json
$Validation = Get-Content -LiteralPath "metadata\driver\figure3_egr1_validation\figure3_egr1_a_to_f_validation_report.json" -Raw | ConvertFrom-Json
$Shared = Import-Csv -Delimiter "`t" "metadata\driver\three_axis_figure_consistency\figure2_figure3_figure4_shared_colour_limits.tsv"

Write-Host "`nFigure 3 EGR1 A-F complete" -ForegroundColor Green
Write-Host "Scripts: scripts\*figure3*egr1*"
Write-Host "Main figures: figures\driver\figure3*_egr1* and figures\driver\figure3a_stress_transition"
Write-Host "Source data: metadata\driver\figure3*_egr1* and data\processed\driver\figure3e_egr1_sctenifoldknk"
Write-Host "CellOracle source: metadata\driver\celloracle_module6_8\celloracle_module6_8_cell_shift_summary.tsv.gz (EGR1=0)"
Write-Host "Stress-transition network: stressed_injured + regenerative_progenitor, 3,000 genes"
Write-Host "Figure 3E significant genes: $($EReport.n_significant_excluding_target)"
Write-Host "Figure 3F significant pathways: $($FReport.n_significant_pathways)"
Write-Host "Figure 3D peak absolute stage: $($DReport.observed_umap_pattern.absolute_effect_peak_stage)"
Write-Host "Three-axis shared symmetric limit: $($Shared[0].shared_symmetric_limit)"
Write-Host "Review-risk flags: $($Validation.review_risk_flags.Count)"
Write-Host "SCI assessment: $($Validation.sci_main_figure_assessment)"
Write-Host "Protected SOX4/HNF4A unchanged: $($Validation.protected_assets_unchanged)"
