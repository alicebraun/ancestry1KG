#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --output=tmp/ancestry1KG_%j.out
#SBATCH --error=tmp/ancestry1KG_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH --mem=32G

submit_dir="$(pwd)"

set -euo pipefail

echo "=============================================="
echo "  1KG ADMIXTURE Ancestry Inference Pipeline   "
echo "=============================================="

# ---- Parse inputs (non-interactive) ---- #
if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <plink_input_prefix> <output_prefix>"
  exit 1
fi

mkdir -p tmp
plink_input="$1"
prefix="$2"

THREADS=16
K=26

# ---- Validate input files ---- #
for ext in bed bim fam; do
  if [ ! -f "${plink_input}.${ext}" ]; then
    echo "Error: ${plink_input}.${ext} not found!"
    exit 1
  fi
done

# ---- Step 1: Merge ---- #
merged="1kg_${prefix}.merged"
if [ ! -f "${merged}.bed" ]; then
  echo "Merging with 1KG reference..."
  plink --bfile 1KG_high_coverage_20130606_g1k_3202.merged \
        --bmerge "${plink_input}" \
        --make-bed \
        --allow-no-sex \
        --out "$merged" || {
          echo "PLINK merge failed (likely due to triallelic SNPs)."
          exit 1
        }
else
  echo "Skipping merge — output ${merged}.bed already exists."
fi

# ---- Step 2: Missingness filter ---- #
filtered="${merged}.geno.05"
if [ ! -f "${filtered}.bed" ]; then
  echo "Checking SNP missingness..."
  plink --bfile "$merged" --missing --out "$merged"

  missing_count=$(awk '$5 > 0.05' "${merged}.lmiss" | wc -l)
  echo "Found $missing_count SNPs with missing rate > 5% — filtering..."

  plink --bfile "$merged" \
        --geno 0.05 \
        --make-bed \
        --allow-no-sex \
        --out "$filtered"
else
  echo "Skipping missingness filtering — output ${filtered}.bed already exists."
fi

# ---- Step 3: Pruning ---- #
pruned="${filtered}.pruned"
if [ ! -f "${pruned}.bed" ]; then
  echo "Pruning SNPs..."
  plink --bfile "$filtered" \
        --indep-pairwise 50 10 0.1

  plink --bfile "$filtered" \
        --extract plink.prune.in \
        --make-bed \
        --allow-no-sex \
        --out "$pruned"
else
  echo "Skipping pruning — output ${pruned}.bed already exists."
fi

# ---- Step 4: Create .pop file ---- #
popfile="${pruned}.pop"
if [ ! -f "$popfile" ]; then
  echo "Creating .pop file..."
  awk '{
    if ($1 ~ /con/ || $1 ~ /cas/) {
      print "-";
    } else {
      print $1;
    }
  }' "${pruned}.fam" > "$popfile"
else
  echo "Skipping .pop file creation — $popfile already exists."
fi

# ---- Step 5: Run ADMIXTURE ---- #
qfile="${pruned}.${K}.Q"
if [ ! -f "$qfile" ]; then
  echo "Running ADMIXTURE with K=$K..."
  admixture --supervised --seed 666 -C 1 -j${THREADS} "${pruned}.bed" $K
else
  echo "Skipping ADMIXTURE — $qfile already exists."
fi

echo "Finished running ADMIXTURE."
echo "Outputs:"
echo "- ADMIXTURE Q file: $qfile"
echo "- Pop file:         $popfile"
echo "- BED prefix:       $pruned"

# Run the R script
echo "Loading R..."
module load 2024
module load Anaconda3/2024.06-1

# Correct way to enable conda commands
source /sw/arch/RHEL9/EB_production/2024/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh
source activate rp_env

Rscript  2_ancestry_inference.R  "$submit_dir" "${prefix}.merged.geno.05.pruned" 

# Check output
if [ -f "1kg_${prefix}.merged.geno.05.pruned.Q.annotated" ]; then
  echo "R script completed successfully. Output: 1kg_${prefix}.merged.geno.05.pruned.Q.annotated"
else
  echo "R script failed or output file missing."
  exit 1
fi