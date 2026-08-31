# =============================================================================
# SMR and HEIDI analysis
# 2 inguinal hernia GWAS datasets x 4 GTEx v8 tissues = 8 analyses
#
# Requirements:
#   - SMR v1.3.1 executable available in PATH, or set the SMR_EXE
#     environment variable to the full executable path.
#   - 1000 Genomes European LD reference panel.
#   - UKB and FinnGen GWAS summary files prepared by
#     03_GWAS_Preparation.R.
#   - GTEx v8 BESD-lite files for the four tissues listed below.
# =============================================================================

param(
    [switch]$PreflightOnly
)
$ErrorActionPreference = "Stop"

# ===== Repository paths =====
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

# By default, inputs are expected under the repository data/ directory.
# Local paths can be overridden with environment variables so that large
# third-party datasets do not need to be copied into the repository.

if ($env:SMR_GWAS_DIR) {
    $gwasDir = $env:SMR_GWAS_DIR
} else {
    $gwasDir = Join-Path $repoRoot "data/GWAS"
}

if ($env:SMR_EQTL_DIR) {
    $eqtlDir = $env:SMR_EQTL_DIR
} else {
    $eqtlDir = Join-Path $repoRoot "data/GTEx_v8/eQTL_besd_lite"
}

if ($env:SMR_LD_DIR) {
    $ldDir = $env:SMR_LD_DIR
} else {
    $ldDir = Join-Path $repoRoot "data/1000Genomes"
}

$outDir = Join-Path $repoRoot "processed_results/06_SMR_HEIDI"

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

# ===== SMR executable =====
# If SMR_EXE is defined, use it; otherwise expect this executable in PATH.
if ($env:SMR_EXE) {
    $smr = $env:SMR_EXE
} else {
    $smr = "smr-1.3.1-win.exe"
}

# ===== Check SMR executable =====
if (Test-Path $smr) {
    $smr = (Resolve-Path $smr).Path
} elseif (-not (Get-Command $smr -ErrorAction SilentlyContinue)) {
    throw "SMR executable not found: $smr"
}

# ===== LD reference =====
$ldReference = Join-Path $ldDir "EUR_1000G"

# ===== GWAS datasets =====
# Individual GWAS files can also be overridden directly.
if ($env:SMR_UKB_GWAS) {
    $ukbGwas = $env:SMR_UKB_GWAS
} else {
    $ukbGwas = Join-Path $gwasDir "UKB/UKB_gwas_for_smr.txt"
}

if ($env:SMR_FINNGEN_GWAS) {
    $finnGwas = $env:SMR_FINNGEN_GWAS
} else {
    $finnGwas = Join-Path $gwasDir "FinnGen/Finn_gwas_for_smr.txt"
}

$studies = @(
    @{
        Name = "UKB"
        File = $ukbGwas
    },
    @{
        Name = "FinnGen"
        File = $finnGwas
    }
)

# ===== GTEx v8 tissues =====
$tissues = @(
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Muscle_Skeletal",
    "Cells_Cultured_fibroblasts"
)

# ===== Check LD reference files =====
foreach ($ext in @(".bed", ".bim", ".fam")) {
    $file = "${ldReference}${ext}"
    if (-not (Test-Path $file)) {
        throw "Missing LD reference file: $file"
    }
}

# ===== Run all 8 SMR/HEIDI analyses =====
foreach ($study in $studies) {

    $studyName = $study.Name
    $gwasFile  = $study.File

    if (-not (Test-Path $gwasFile)) {
        throw "Missing GWAS summary file: $gwasFile"
    }

    foreach ($tissue in $tissues) {

        $eqtlPrefix = Join-Path $eqtlDir "${tissue}.lite"

        foreach ($ext in @(".besd", ".epi", ".esi")) {
            $file = "${eqtlPrefix}${ext}"
            if (-not (Test-Path $file)) {
                throw "Missing GTEx BESD-lite file: $file"
            }
        }

        $outPrefix = Join-Path $outDir "${studyName}_${tissue}"

        Write-Host ""
        Write-Host "============================================================"
        Write-Host "Running SMR/HEIDI: $studyName x $tissue"
        Write-Host "============================================================"

if ($PreflightOnly) {
    Write-Host "Preflight OK: $studyName x $tissue"
    continue
}

        & $smr `
            --bfile $ldReference `
            --gwas-summary $gwasFile `
            --beqtl-summary $eqtlPrefix `
            --maf 0.01 `
            --diff-freq 0.2 `
            --peqtl-smr 5e-8 `
            --peqtl-heidi 1.57e-3 `
            --ld-upper-limit 0.9 `
            --ld-lower-limit 0.05 `
            --heidi-min-m 3 `
            --heidi-max-m 20 `
            --cis-wind 1000000 `
            --thread-num 5 `
            --diff-freq-prop 0.05 `
            --out $outPrefix

        if ($LASTEXITCODE -ne 0) {
            throw "SMR failed for $studyName x $tissue with exit code $LASTEXITCODE"
        }
    }
}

Write-Host ""

if ($PreflightOnly) {
    Write-Host "Preflight passed: all required inputs for the 8 SMR/HEIDI analyses were found."
} else {
    Write-Host "All 8 SMR/HEIDI analyses completed successfully."
}