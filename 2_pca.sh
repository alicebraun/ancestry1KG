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
echo "  PCA using RICOPILI pcaer                    "
echo "=============================================="

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <plink_input_prefix> <outname>"
  echo "  <plink_input_prefix> : prefix for the pruned PLINK files"
  echo "                         (expects {prefix}.geno.05.pruned.bed/bim/fam)"
  echo "  <outname>            : identifier for pcaer run"
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

# ---- Initialise QC directory and run pcaer ---- #
# Use absolute path for bim since pcaer runs from qc/
abs_bim="$(pwd)/${plink_input}.geno.05.pruned.bim"
id_tager_2 --create --nn scz_${outname}_mix ${plink_input}.geno.05.pruned.fam
cd qc
pcaer "${abs_bim}" --out ${outname}
cd ..

# ---- Copy PCA output to working directory ---- #
mds_file=$(ls qc/${outname}*.menv.mds_cov 2>/dev/null | head -1 || true)

if [ -z "$mds_file" ]; then
  echo "Warning: no .menv.mds_cov file found in qc/ — check pcaer output manually."
else
  cp "$mds_file" .
  mds_base="$(basename "$mds_file")"
  echo ""
  echo "PCA complete."
  echo "  Copied: ${mds_base}"
  echo ""
  echo "Next step — generate plots:"
  echo "  Rscript 3_plots.R \$(pwd) <Q_annotated_prefix> ${mds_base}"
fi
