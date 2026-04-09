#!/bin/bash
#SBATCH --job-name=split_ancestry
#SBATCH --output=tmp/split_%j.out
#SBATCH --error=tmp/split_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# ============================================================ #
#  4_split_by_ancestry.sh
#
#  Split a PLINK dataset by inferred superpopulation ancestry.
#
#  Usage:
#    sbatch 4_split_by_ancestry.sh <plink_prefix> <outname>
#
#  Arguments:
#    plink_prefix  — PLINK file prefix for the cohort to split
#                    (e.g. CLOZUK2_merged.hg38.ch.fl)
#    outname       — the outname used in step 1
#                    (e.g. clz2a → reads 1kg_clz2a.geno.05.pruned.Q.annotated)
#
#  Output:
#    {plink_prefix}_{SUPERPOP}.bed/bim/fam  for each group
#    keep files in split_keep/
#
#  Superpop groups: AFR AMR EAS EUR SAS MIX
#  Only target samples (Ancestry_Pheno == "-") are split.
#  1KG reference samples are excluded.
# ============================================================ #

set -euo pipefail

plink_input="${1:?Usage: $0 <plink_prefix> <outname>}"
outname="${2:?Usage: $0 <plink_prefix> <outname>}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
annotated="${SCRIPT_DIR}/1kg_${outname}.geno.05.pruned.Q.annotated"
keep_dir="${SCRIPT_DIR}/split_keep"

cd "$SCRIPT_DIR"
mkdir -p "$keep_dir" tmp

echo "PLINK input : $plink_input"
echo "Q.annotated : $annotated"
echo "Keep dir    : $keep_dir"

# Check inputs
if [[ ! -f "${plink_input}.bed" ]]; then
    echo "ERROR: ${plink_input}.bed not found" >&2; exit 1
fi
if [[ ! -f "$annotated" ]]; then
    echo "ERROR: $annotated not found" >&2; exit 1
fi

# ------------------------------------------------------------ #
#  Extract keep files per superpopulation
#  Header: FID IID Ancestry_Pheno Ancestry_Inf ...
#  Col 1=FID, Col 2=IID, Col 3=Ancestry_Pheno, Col 4=Ancestry_Inf
#  Only keep rows where Ancestry_Pheno == "-" (target samples)
#  Superpop = first word before "_" in Ancestry_Inf, or full value if no "_"
# ------------------------------------------------------------ #

echo "Extracting keep lists..."

awk -v keep_dir="$keep_dir" '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "Ancestry_Pheno") pheno_col = i
      if ($i == "Ancestry_Inf")   inf_col   = i
    }
    next
  }
  $pheno_col != "-" { next }
  {
    split($inf_col, a, "_");
    superpop = a[1];
    print $1, $2 > (keep_dir "/" superpop ".keep")
  }' "$annotated"

# Print summary
echo "Keep file counts:"
shopt -s nullglob
keep_files=("$keep_dir"/*.keep)
if [[ ${#keep_files[@]} -eq 0 ]]; then
    echo "  WARNING: no keep files written — check Ancestry_Pheno/Ancestry_Inf columns in $annotated"
else
    for f in "${keep_files[@]}"; do
        superpop=$(basename "$f" .keep)
        n=$(wc -l < "$f")
        echo "  ${superpop}: ${n} individuals"
    done
fi
shopt -u nullglob

# ------------------------------------------------------------ #
#  Split with PLINK
# ------------------------------------------------------------ #

shopt -s nullglob
for keepfile in "$keep_dir"/*.keep; do
    superpop=$(basename "$keepfile" .keep)
    n=$(wc -l < "$keepfile")

    if [[ "$n" -eq 0 ]]; then
        echo "Skipping ${superpop}: no individuals"
        continue
    fi

    outfile="${plink_input}_${superpop}"
    echo "Writing ${outfile} (n=${n})..."

    plink \
        --bfile "$plink_input" \
        --keep  "$keepfile" \
        --make-bed \
        --out   "$outfile" \
        --allow-no-sex \
        --silent
done

echo ""
echo "Done. Output files:"
ls -lh "${plink_input}"_*.bed 2>/dev/null || echo "  (none found)"
