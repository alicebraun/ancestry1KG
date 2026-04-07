#!/bin/bash
#SBATCH --job-name=pca_1kg
#SBATCH --output=tmp/pca_%j.out
#SBATCH --error=tmp/pca_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G

set -euo pipefail

echo "=============================================="
echo "  PCA using PLINK --pca                       "
echo "=============================================="

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <plink_input_prefix> <outname>"
  echo "  <plink_input_prefix> : prefix for the pruned PLINK files"
  echo "                         (expects {prefix}.geno.05.pruned.bed/bim/fam)"
  echo "  <outname>            : base name for PCA output files"
  exit 1
fi

plink_input="$1"
outname="$2"

mkdir -p tmp

# ---- Validate input files ---- #
for ext in bed bim fam; do
  if [ ! -f "${plink_input}.geno.05.pruned.${ext}" ]; then
    echo "Error: ${plink_input}.geno.05.pruned.${ext} not found!"
    exit 1
  fi
done

# ---- Run PLINK PCA (10 components) ---- #
echo "Running PLINK PCA on ${plink_input}.geno.05.pruned ..."
plink --bfile "${plink_input}.geno.05.pruned" \
      --pca 10 \
      --allow-no-sex \
      --out "${outname}"

# ---- Normalise eigenvec header ---- #
# PLINK 1.9 outputs no header; PLINK 2 outputs "#FID IID PC1 ...".
# Downstream R scripts (3_plots.R) expect column names C1, C2, ..., C10.
first_field=$(head -1 "${outname}.eigenvec" | awk '{print $1}')

if echo "${first_field}" | grep -qiE '^#?FID$'; then
  # PLINK2: has a header — strip leading '#' and rename PC -> C
  sed '1s/^#//; 1s/\bPC\([0-9]*\)\b/C\1/g' "${outname}.eigenvec" \
    > tmp/${outname}_eigenvec_tmp
else
  # PLINK1: no header — prepend one
  { echo "FID IID C1 C2 C3 C4 C5 C6 C7 C8 C9 C10"
    cat "${outname}.eigenvec"
  } > tmp/${outname}_eigenvec_tmp
fi

mv tmp/${outname}_eigenvec_tmp "${outname}.eigenvec"

echo ""
echo "PCA complete."
echo "  Eigenvectors : ${outname}.eigenvec  (columns FID IID C1-C10)"
echo "  Eigenvalues  : ${outname}.eigenval"
echo ""
echo "Next step — generate plots:"
echo "  Rscript 3_plots.R \$(pwd) <Q_annotated_prefix> ${outname}.eigenvec"