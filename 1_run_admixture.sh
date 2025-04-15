#!/bin/bash

set -euo pipefail

echo "=============================================="
echo "  1KG ADMIXTURE Ancestry Inference Pipeline   "
echo "=============================================="

# Prompt for user inputs
read -rp "Enter your PLINK input prefix (e.g. genotypes.bed, omitting .bed): " plink_input
read -rp "Enter an output prefix : " prefix

THREADS=16
K=26

# Check if input files exist
for ext in bed bim fam; do
  if [ ! -f "${plink_input}.${ext}" ]; then
    echo "Error: ${plink_input}.${ext} not found!"
    exit 1
  fi
done

# ---- Step 1–3: Merge, QC, Prune, Create .pop ---- #

echo "Merging with 1KG reference..."
plink --bfile 1KG_high_coverage_20130606_g1k_3202.merged \
      --bmerge "${plink_input}" \
      --make-bed \
      --out "1kg_${prefix}.merged" || {
        echo "PLINK merge failed (likely due to triallelic SNPs)."
        echo "See: https://www.cog-genomics.org/plink/1.9/data#merge"
        exit 1
      }

echo "Checking SNP missingness..."
plink --bfile "1kg_${prefix}.merged" --missing --out "1kg_${prefix}.merged"

missing_count=$(awk '$5 > 0.05' "1kg_${prefix}.merged.lmiss" | wc -l)
echo "Found $missing_count SNPs with missing rate > 5% — filtering..."

plink --bfile "1kg_${prefix}.merged" \
      --geno 0.05 \
      --make-bed \
      --out "1kg_${prefix}.geno.05.merged"

echo "Pruning SNPs..."
plink --bfile "1kg_${prefix}.geno.05.merged" \
      --indep-pairwise 50 10 0.1

plink --bfile "1kg_${prefix}.geno.05.merged" \
      --extract plink.prune.in \
      --make-bed \
      --out "1kg_${prefix}.geno.05.merged.pruned"

echo "Creating .pop file based on RICOPILI formatted FID..."
awk '{
  if ($1 ~ /con/ || $1 ~ /cas/) {
    print "-";
  } else {
    print $1;
  }
}' "1kg_${prefix}.geno.05.merged.pruned.fam" > "1kg_${prefix}.geno.05.merged.pruned.pop"

# ---- Step 4: Run ADMIXTURE ---- #

echo "Running ADMIXTURE with K=$K..."
admixture --supervised --seed 666 -C 10 -j${THREADS} \
          "1kg_${prefix}.geno.05.merged.pruned.bed" $K

echo "Finished."
echo "Outputs:"
echo "- ADMIXTURE Q file: 1kg_${prefix}.geno.05.merged.pruned.${K}.Q"
echo "- Pop file:         1kg_${prefix}.geno.05.merged.pruned.pop"
echo "- BED prefix:       1kg_${prefix}.geno.05.merged.pruned"
echo "Next step: Run the R script to assign labels and visualize ancestry."