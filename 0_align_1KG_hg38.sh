#!/bin/bash
#SBATCH --job-name=align_1kg
#SBATCH --output=tmp/align_%j.out
#SBATCH --error=tmp/align_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G

set -euo pipefail

echo "=============================================="
echo "  1KG hg38 Alignment (impute_dirsub)          "
echo "=============================================="

# ---- Parse inputs ---- #
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

echo "Running impute_dirsub --onlyalign from: $plink_dir"
cd "$plink_dir"

impute_dirsub \
  --refdir /gpfs/work5/0/pgcdac/imputation_references/1KG_high_coverage_mm4_jan2025/ \
  --outname "$outname" \
  --onlyalign \
  --serial

# ---- Symlink aligned output files into ancestry1KG directory ---- #
echo "Looking for aligned files in pi_sub/..."

aligned_files=$(ls pi_sub/${outname}*.hg38.ch.fl.bed 2>/dev/null | sed 's/\.bed$//' || true)

if [ -z "$aligned_files" ]; then
  echo "Error: No files matching pi_sub/${outname}*.hg38.ch.fl.bed found."
  echo "Contents of pi_sub/:"
  ls pi_sub/ || true
  exit 1
fi

cd "$ancestry_dir"

echo "Symlinking aligned files into $(pwd)..."
for f in $aligned_files; do
  base="$(basename "$f")"
  for ext in bed bim fam; do
    ln -sf "${plink_dir}/pi_sub/${base}.${ext}" "${base}.${ext}"
    echo "  Linked: ${base}.${ext}"
  done
done

echo ""
echo "Done. Aligned file prefix(es):"
for f in $aligned_files; do
  echo "  $(basename "$f")"
done
echo ""
echo "Next step — run the ADMIXTURE pipeline:"
echo "  sbatch 1_run_admixture.sh <aligned_prefix> <output_prefix>"
