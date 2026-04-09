#!/bin/bash
#SBATCH --job-name=align_1kg
#SBATCH --output=tmp/align_%j.out
#SBATCH --error=tmp/align_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=16G

set -euo pipefail

echo "=============================================="
echo "  1KG hg38 Alignment (impute_dirsub)          "
echo "=============================================="

if [ "$#" -ne 2 ]; then
  echo "Usage: sbatch $0 <plink_input_prefix> <outname>"
  echo "  <plink_input_prefix>  : full path prefix to input .bed/.bim/.fam"
  echo "  <outname>             : identifier for impute_dirsub run"
  exit 1
fi

plink_input="$1"
outname="$2"
ancestry_dir="$(pwd)"

mkdir -p tmp

# ---- Symlink 1KG reference files if not already present ---- #
REF_SOURCE="/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pgcdrc/scz/working/1KG_pca/2_output"
for ext in bed bim fam; do
  ref="1KG_high_coverage_20130606_g1k_3202.merged.${ext}"
  if [[ ! -f "$ref" && ! -L "$ref" ]]; then
    echo "Symlinking ${ref}..."
    ln -s "${REF_SOURCE}/${ref}" .
  fi
done

# ---- Validate input files ---- #
for ext in bed bim fam; do
  if [ ! -f "${plink_input}.${ext}" ]; then
    echo "Error: ${plink_input}.${ext} not found!"
    exit 1
  fi
done

# ---- Run impute_dirsub alignment ---- #
# Must run from the directory containing the PLINK files
plink_dir="$(dirname "$plink_input")"

plink_base="$(basename "$plink_input")"
cd "$plink_dir"

if [ -f "pi_sub/${plink_base}.hg38.ch.fl.bed" ]; then
  echo "Skipping impute_dirsub — pi_sub/${plink_base}.hg38.ch.fl.bed already exists."
else
  echo "Running impute_dirsub --onlyalign from: $plink_dir"
  impute_dirsub \
    --refdir /gpfs/work5/0/pgcdac/imputation_references/1KG_high_coverage_mm4_jan2025/ \
    --outname "$outname" \
    --onlyalign
fi

# ---- Wait for impute_dirsub sub-jobs to complete ---- #
# impute_dirsub spawns its own SLURM sub-jobs and returns immediately.
# Poll for flipping_done (last stage) before proceeding.
echo "Waiting for impute_dirsub sub-jobs (polling for flipping_done)..."
MAX_WAIT=25200  # 7 hours
WAITED=0
while [ ! -f "flipping_done" ]; do
  sleep 120
  WAITED=$((WAITED + 120))
  echo "  ...still waiting (${WAITED}s elapsed)"
  if [ "$WAITED" -ge "$MAX_WAIT" ]; then
    echo "Error: timed out waiting for flipping_done after ${MAX_WAIT}s"
    exit 1
  fi
done
echo "impute_dirsub complete."

# ---- Symlink aligned output files back into working directory ---- #
# impute_dirsub names output files after the input filename, not --outname.
# Only symlink the input dataset (not the reference which also gets aligned).
echo "Looking for aligned files in pi_sub/..."

aligned_files=$(ls "pi_sub/${plink_base}.hg38.ch.fl.bed" 2>/dev/null | sed 's/\.bed$//' || true)

if [ -z "$aligned_files" ]; then
  echo "Error: pi_sub/${plink_base}.hg38.ch.fl.bed not found."
  echo "Contents of pi_sub/:"
  ls pi_sub/ || true
  exit 1
fi

cd "$ancestry_dir"

echo "Symlinking aligned files into $(pwd)..."
for f in $aligned_files; do
  base="$(basename "$f")"
  for ext in bed bim fam; do
    ln -sf "pi_sub/${base}.${ext}" "${base}.${ext}"
    echo "  Linked: ${base}.${ext}"
  done
done

echo ""
echo "Done. Aligned file prefix(es):"
for f in $aligned_files; do
  echo "  $(basename "$f")"
done
echo ""
echo "Next step:"
echo "  sbatch 1_run_admixture.sh <aligned_prefix> <output_prefix>"
