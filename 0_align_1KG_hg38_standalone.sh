#!/bin/bash
#SBATCH --job-name=align_1kg
#SBATCH --output=tmp/align_%j.out
#SBATCH --error=tmp/align_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G

set -euo pipefail

# -----------------------------------------------------------------------
# Configuration — set these paths for your HPC environment
# -----------------------------------------------------------------------
# SLURM copies the script to its spool, so BASH_SOURCE[0] points to the wrong place.
# Use SLURM_SUBMIT_DIR (the directory sbatch/srun was called from) instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
REF_BIM="1KG_high_coverage_20130606_g1k_3202.merged.bim"  # symlink or copy in working directory
CHAIN_FILE="/gpfs/work5/0/pgcdac/ricopili_download/dependencies/liftover/hg19ToHg38.over.chain.gz"
LIFTOVER_BIN="/gpfs/work5/0/pgcdac/ricopili_download/dependencies/liftover/liftOver"

# Activate conda so Rscript is available
module load 2024
module load Anaconda3/2024.06-1
source /sw/arch/RHEL9/EB_production/2024/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh
conda activate admix_r
# -----------------------------------------------------------------------

echo "=============================================="
echo "  1KG hg38 Alignment (standalone)             "
echo "=============================================="

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: sbatch $0 <plink_input_prefix> <outname> [hg19|hg38]"
  echo "  <plink_input_prefix> : prefix of input .bed/.bim/.fam files"
  echo "  <outname>            : base name for output files"
  echo "  [hg19|hg38]          : genome build of input (default: hg19)"
  exit 1
fi

plink_input="$1"
outname="$2"
input_build="${3:-hg19}"
ancestry_dir="$(pwd)"

mkdir -p tmp

# ---- Validate inputs ---- #
for ext in bed bim fam; do
  [ -f "${plink_input}.${ext}" ] || { echo "Error: ${plink_input}.${ext} not found!"; exit 1; }
done
[ -f "${REF_BIM}" ] || { echo "Error: Reference BIM not found: ${REF_BIM}"; exit 1; }


# ====================================================================== #
# Step 1: LiftOver hg19 → hg38 (skip if input is already hg38)
# ====================================================================== #
if [ "${input_build}" == "hg19" ]; then
  echo ""
  echo "Step 1: LiftOver hg19 → hg38..."
  [ -f "${CHAIN_FILE}" ] || { echo "Error: Chain file not found: ${CHAIN_FILE}"; exit 1; }

  # Convert BIM to 0-based BED format required by liftOver
  awk '{
    chr = ($1 ~ /^chr/) ? $1 : "chr"$1
    print chr"\t"($4 - 1)"\t"$4"\t"$2
  }' "${plink_input}.bim" > tmp/${outname}_lo_input.bed

  ${LIFTOVER_BIN} \
    tmp/${outname}_lo_input.bed \
    "${CHAIN_FILE}" \
    tmp/${outname}_lo_output.bed \
    tmp/${outname}_lo_unmapped.bed

  # Variants that failed liftover
  grep -v '^#' tmp/${outname}_lo_unmapped.bed \
    | awk '{print $4}' > tmp/${outname}_unmapped.txt || true

  # Also exclude variants that lifted to unplaced/random contigs (e.g. chr4_GL000008v2_random)
  # PLINK only accepts autosomes 1-22, X, Y, XY, MT
  awk '{
    chr = $1; sub(/^chr/, "", chr)
    if (chr ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|XY|MT)$/) print $4"\t"chr
    else print $4 > "/dev/stderr"
  }' tmp/${outname}_lo_output.bed \
    > tmp/${outname}_update_chr.txt \
    2> tmp/${outname}_contig_exclude.txt

  # Merge unmapped + unplaced-contig variants into a single exclude list
  cat tmp/${outname}_unmapped.txt tmp/${outname}_contig_exclude.txt \
    | sort -u > tmp/${outname}_exclude_prelift.txt

  n_unmapped=$(wc -l < tmp/${outname}_exclude_prelift.txt)
  echo "  Variants excluded (unmapped or non-standard contig): ${n_unmapped}"

  # update-map: only for variants that passed the chr filter
  # update_chr.txt has 2 cols: SNP_ID chr  →  use $1 as key
  # lo_output.bed has 4 cols:  chr start end SNP_ID  →  $4 is SNP_ID, $3 is pos
  awk 'NR==FNR{keep[$1]=1; next} keep[$4]{print $4"\t"$3}' \
    tmp/${outname}_update_chr.txt \
    tmp/${outname}_lo_output.bed > tmp/${outname}_update_map.txt

  plink --bfile "${plink_input}" \
        --exclude tmp/${outname}_exclude_prelift.txt \
        --update-chr tmp/${outname}_update_chr.txt \
        --update-map tmp/${outname}_update_map.txt \
        --allow-no-sex \
        --make-bed \
        --out tmp/${outname}_lifted

  lifted_prefix="tmp/${outname}_lifted"
else
  echo "Step 1: Skipping liftOver (input build: ${input_build})."
  lifted_prefix="${plink_input}"
fi


# ====================================================================== #
# Step 2: Update variant IDs to chr:pos format
# ====================================================================== #
echo ""
echo "Step 2: Updating variant IDs to chr:pos..."

# Build an --update-name file: old_id <tab> new_id
paste \
  <(awk '{print $2}' "${lifted_prefix}.bim") \
  <(awk '{print $1":"$4}' "${lifted_prefix}.bim") \
  > tmp/${outname}_update_ids.txt

plink --bfile "${lifted_prefix}" \
      --update-name tmp/${outname}_update_ids.txt \
      --allow-no-sex \
      --make-bed \
      --out tmp/${outname}_renamed

# Remove duplicate chr:pos IDs introduced by liftover
# (keep only the first occurrence; exclude all others)
awk 'seen[$2]++ {print $2}' tmp/${outname}_renamed.bim \
  > tmp/${outname}_dup_ids.txt

n_dup=$(wc -l < tmp/${outname}_dup_ids.txt)
echo "  Duplicate chr:pos variants to exclude: ${n_dup}"

if [ "${n_dup}" -gt 0 ]; then
  plink --bfile tmp/${outname}_renamed \
        --exclude tmp/${outname}_dup_ids.txt \
        --allow-no-sex \
        --make-bed \
        --out tmp/${outname}_renamed_dedup
  mv tmp/${outname}_renamed_dedup.bed tmp/${outname}_renamed.bed
  mv tmp/${outname}_renamed_dedup.bim tmp/${outname}_renamed.bim
  mv tmp/${outname}_renamed_dedup.fam tmp/${outname}_renamed.fam
fi


# ====================================================================== #
# Step 3: Harmonize alleles against 1KG reference
# ====================================================================== #
echo ""
echo "Step 3: Harmonizing alleles against 1KG reference..."

Rscript "${SCRIPT_DIR}/0b_harmonize_alleles.R" \
  --input-bim   tmp/${outname}_renamed.bim \
  --ref-bim     "${REF_BIM}" \
  --out-flip    tmp/${outname}_flip.txt \
  --out-exclude tmp/${outname}_exclude.txt

plink --bfile  tmp/${outname}_renamed \
      --flip   tmp/${outname}_flip.txt \
      --exclude tmp/${outname}_exclude.txt \
      --allow-no-sex \
      --make-bed \
      --out "${outname}.hg38.ch.fl"

echo ""
echo "======================================================"
echo "Done. Aligned output prefix: ${outname}.hg38.ch.fl"
echo "======================================================"
echo ""
echo "Next step:"
echo "  sbatch 1_run_admixture.sh ${outname}.hg38.ch.fl <output_prefix>"
