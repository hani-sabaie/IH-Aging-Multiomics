# Reproducible GCTA-COJO analysis for UK Biobank and FinnGen SMAD3 loci.

$ErrorActionPreference = "Stop"

# ===== Repository root =====
$repoRoot = Split-Path -Parent $PSScriptRoot

# ===== External executable =====
$gcta = if ($env:GCTA_EXE) {
    $env:GCTA_EXE
} else {
    "gcta64"
}

# ===== Input/reference directories =====
$refdir = if ($env:REF_1000G_DIR) {
    $env:REF_1000G_DIR
} else {
    Join-Path $repoRoot "data/1000Genomes"
}

$ukbDir = Join-Path $repoRoot "data/GWAS/UKB"
$finDir = Join-Path $repoRoot "data/GWAS/FinnGen"

$outDir = Join-Path $repoRoot "processed_results/07_GCTA_COJO"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

# ===== Input files =====
$ukbCojo = Join-Path $ukbDir "UKB_SMAD3_for_GCTA.tsv"
$finCojo = Join-Path $finDir "Finn_SMAD3_for_GCTA.tsv"

$ukbRef = Join-Path $refdir "UKB_SMAD3_1000G_EUR"
$finRef = Join-Path $refdir "Finn_SMAD3_1000G_EUR"

# ===== Required-input checks =====
$required = @(
    $ukbCojo,
    $finCojo,
    "$ukbRef.bed",
    "$ukbRef.bim",
    "$ukbRef.fam",
    "$finRef.bed",
    "$finRef.bim",
    "$finRef.fam"
)

foreach ($f in $required) {
    if (-not (Test-Path $f)) {
        throw "Required input not found: $f"
    }
}

# ===== UK Biobank =====
$ukbOut = Join-Path $outDir "UKB_SMAD3_GCTA"

& $gcta `
    --bfile $ukbRef `
    --cojo-file $ukbCojo `
    --cojo-slct `
    --cojo-p 5e-5 `
    --out $ukbOut

if ($LASTEXITCODE -ne 0) {
    throw "GCTA-COJO failed for UK Biobank."
}

# ===== FinnGen =====
$finOut = Join-Path $outDir "Finn_SMAD3_GCTA"

& $gcta `
    --bfile $finRef `
    --cojo-file $finCojo `
    --cojo-slct `
    --cojo-p 5e-5 `
    --out $finOut

if ($LASTEXITCODE -ne 0) {
    throw "GCTA-COJO failed for FinnGen."
}

Write-Host "Completed GCTA-COJO for UK Biobank and FinnGen."
