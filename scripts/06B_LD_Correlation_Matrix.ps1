# Reproducible construction of locus-specific 1000 Genomes EUR LD references
# for UK Biobank and FinnGen SMAD3 analyses.

$ErrorActionPreference = "Stop"

# ===== Repository root =====
$repoRoot = Split-Path -Parent $PSScriptRoot

# ===== External executables =====
$plink2 = if ($env:PLINK2_EXE) {
    $env:PLINK2_EXE
} else {
    "plink2"
}

$plink = if ($env:PLINK_EXE) {
    $env:PLINK_EXE
} else {
    "plink"
}

# ===== Input/reference directories =====
$refdir = if ($env:REF_1000G_DIR) {
    $env:REF_1000G_DIR
} else {
    Join-Path $repoRoot "data/1000Genomes"
}

$ukbDir = Join-Path $repoRoot "data/GWAS/UKB"
$finDir = Join-Path $repoRoot "data/GWAS/FinnGen"

$pfilep  = Join-Path $refdir "chr15_phase3"
$eurKeep = Join-Path $refdir "EUR_samples.keep"

$ukbRs  = Join-Path $ukbDir "SMAD3_UKB_rs.tsv"
$finnRs = Join-Path $finDir "SMAD3_Finn_rs.tsv"

# ===== Required-input checks =====
$required = @(
    "$pfilep.pgen",
    "$pfilep.pvar.zst",
    "$pfilep.psam",
    $eurKeep,
    $ukbRs,
    $finnRs
)

foreach ($f in $required) {
    if (-not (Test-Path $f)) {
        throw "Required input not found: $f"
    }
}

# ===== UK Biobank =====
$ukbPrefix = Join-Path $refdir "UKB_SMAD3_1000G_EUR"

& $plink2 `
    --pfile $pfilep vzs `
    --extract $ukbRs `
    --max-alleles 2 `
    --keep $eurKeep `
    --keep-allele-order `
    --make-bed `
    --out $ukbPrefix

if ($LASTEXITCODE -ne 0) {
    throw "PLINK2 failed while constructing the UKB LD reference."
}

& $plink `
    --bfile $ukbPrefix `
    --keep-allele-order `
    --r square `
    --out "${ukbPrefix}_LD"

if ($LASTEXITCODE -ne 0) {
    throw "PLINK failed while calculating the UKB LD matrix."
}

# ===== FinnGen =====
$finnPrefix = Join-Path $refdir "Finn_SMAD3_1000G_EUR"

& $plink2 `
    --pfile $pfilep vzs `
    --extract $finnRs `
    --max-alleles 2 `
    --keep $eurKeep `
    --keep-allele-order `
    --make-bed `
    --out $finnPrefix

if ($LASTEXITCODE -ne 0) {
    throw "PLINK2 failed while constructing the FinnGen LD reference."
}

& $plink `
    --bfile $finnPrefix `
    --keep-allele-order `
    --r square `
    --out "${finnPrefix}_LD"

if ($LASTEXITCODE -ne 0) {
    throw "PLINK failed while calculating the FinnGen LD matrix."
}

Write-Host "Completed UKB and FinnGen SMAD3 LD-reference construction."
