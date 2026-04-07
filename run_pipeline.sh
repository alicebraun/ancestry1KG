#!/bin/bash
# ============================================================ #
#  ancestry1KG — full pipeline wrapper
#
#  Usage:
#    bash run_pipeline.sh --input <plink_prefix> \
#                         --outname <cohort_name> \
#                         --mode [ricopili|standalone] \
#                         [--build hg19|hg38]
#                         [--check]   # only run dependency checks
#
#  Example:
#    bash run_pipeline.sh --input CLOZUK2_merged \
#                         --outname clozuk2 \
#                         --mode ricopili
# ============================================================ #

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Defaults ---- #
MODE="ricopili"
BUILD="hg19"
INPUT=""
OUTNAME=""
CHECK_ONLY=false

# ---- Parse arguments ---- #
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)      INPUT="$2";   shift 2 ;;
    --outname)    OUTNAME="$2"; shift 2 ;;
    --mode)       MODE="$2";    shift 2 ;;
    --build)      BUILD="$2";   shift 2 ;;
    --check)      CHECK_ONLY=true; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ "$MODE" != "ricopili" && "$MODE" != "standalone" ]]; then
  echo "Error: --mode must be 'ricopili' or 'standalone'"
  exit 1
fi

# ============================================================ #
#  Dependency check
# ============================================================ #
PASS=true

check_ok()  { echo "  [OK]      $1"; }
check_warn() { echo "  [WARN]    $1"; }
check_fail() { echo "  [MISSING] $1"; PASS=false; }

echo ""
echo "=============================================="
echo "  Dependency check (mode: ${MODE})"
echo "=============================================="

# ---- Pipeline scripts ---- #
echo ""
echo "Pipeline scripts:"
for script in \
    "1_run_admixture.sh" \
    "1b_ancestry_inference.R" \
    "2_pca.sh" \
    "2_pca_standalone.sh" \
    "3_plots.R"; do
  if [ -f "${SCRIPT_DIR}/${script}" ]; then
    check_ok "${SCRIPT_DIR}/${script}"
  else
    check_fail "${SCRIPT_DIR}/${script}"
  fi
done

if [ "$MODE" == "ricopili" ]; then
  for script in "0_align_1KG_hg38.sh"; do
    [ -f "${SCRIPT_DIR}/${script}" ] && check_ok "${SCRIPT_DIR}/${script}" || check_fail "${SCRIPT_DIR}/${script}"
  done
else
  for script in "0_align_1KG_hg38_standalone.sh" "0b_harmonize_alleles.R"; do
    [ -f "${SCRIPT_DIR}/${script}" ] && check_ok "${SCRIPT_DIR}/${script}" || check_fail "${SCRIPT_DIR}/${script}"
  done
fi

# ---- Reference files ---- #
echo ""
echo "Reference files (expected in working directory):"
REF_BASE="1KG_high_coverage_20130606_g1k_3202.merged"
for ext in bed bim fam; do
  f="${REF_BASE}.${ext}"
  if [ -f "$f" ]; then
    size=$(du -sh "$f" | cut -f1)
    check_ok "${f}  (${size})"
  else
    check_fail "${f}"
  fi
done

# ---- System tools ---- #
echo ""
echo "System tools:"
if command -v plink &>/dev/null; then
  ver=$(plink --version 2>&1 | head -1)
  check_ok "plink    — ${ver}"
else
  check_fail "plink    — not found in PATH"
fi

if command -v admixture &>/dev/null; then
  ver=$(admixture 2>&1 | grep -i 'version' | head -1 || echo "found")
  check_ok "admixture — ${ver}"
else
  check_fail "admixture — not found in PATH"
fi

if [ "$MODE" == "ricopili" ]; then
  if command -v impute_dirsub &>/dev/null; then
    check_ok "impute_dirsub — $(which impute_dirsub)"
  else
    check_fail "impute_dirsub — not found in PATH (is RICOPILI loaded?)"
  fi
  if command -v id_tager_2 &>/dev/null; then
    check_ok "id_tager_2    — $(which id_tager_2)"
  else
    check_fail "id_tager_2    — not found in PATH"
  fi
  if command -v pcaer &>/dev/null; then
    check_ok "pcaer         — $(which pcaer)"
  else
    check_fail "pcaer         — not found in PATH"
  fi
else
  LIFTOVER_BIN="/gpfs/work5/0/pgcdac/ricopili_download/dependencies/liftover/liftOver"
  CHAIN_FILE="/gpfs/work5/0/pgcdac/ricopili_download/dependencies/liftover/hg19ToHg38.over.chain.gz"
  if [ -x "${LIFTOVER_BIN}" ]; then
    check_ok "liftOver      — ${LIFTOVER_BIN}"
  else
    check_fail "liftOver      — not found at ${LIFTOVER_BIN}"
  fi
  if [ -f "${CHAIN_FILE}" ]; then
    check_ok "chain file    — ${CHAIN_FILE}"
  else
    check_fail "chain file    — not found at ${CHAIN_FILE}"
  fi
fi

# ---- Conda environment ---- #
echo ""
echo "Conda environment (admix_r):"
if conda env list 2>/dev/null | grep -q "admix_r"; then
  check_ok "admix_r environment found"
  # Check key R packages
  for pkg in dplyr ggplot2 tidyr patchwork data.table scales; do
    if conda run -n admix_r Rscript -e "library(${pkg})" &>/dev/null 2>&1; then
      check_ok "  R package: ${pkg}"
    else
      check_warn "  R package: ${pkg} — not installed in admix_r"
    fi
  done
else
  check_fail "admix_r conda environment not found"
  echo ""
  echo "  To create it, run:"
  echo "    conda env create -f ${SCRIPT_DIR}/environment.yml"
  echo ""
  echo "  environment.yml contents:"
  echo "  ---"
  cat "${SCRIPT_DIR}/environment.yml" | sed 's/^/  /'
fi

# ---- Input files (if --input provided) ---- #
if [ -n "$INPUT" ]; then
  echo ""
  echo "Input PLINK files:"
  for ext in bed bim fam; do
    f="${INPUT}.${ext}"
    if [ -f "$f" ]; then
      size=$(du -sh "$f" | cut -f1)
      check_ok "${f}  (${size})"
    else
      check_fail "${f}"
    fi
  done
fi

# ---- Summary ---- #
echo ""
echo "=============================================="
if $PASS; then
  echo "  All checks passed."
else
  echo "  One or more dependencies are missing — fix above before running."
fi
echo "=============================================="
echo ""

if ! $PASS || $CHECK_ONLY; then
  exit $( $PASS && echo 0 || echo 1 )
fi

# ============================================================ #
#  Validate required args before submitting
# ============================================================ #
if [ -z "$INPUT" ] || [ -z "$OUTNAME" ]; then
  echo "Usage: bash run_pipeline.sh --input <plink_prefix> --outname <cohort_name> [--mode ricopili|standalone] [--build hg19|hg38]"
  exit 1
fi

mkdir -p tmp

echo "=============================================="
echo "  Submitting pipeline"
echo "  Input   : ${INPUT}"
echo "  Outname : ${OUTNAME}"
echo "  Mode    : ${MODE}"
echo "  Build   : ${BUILD}"
echo "=============================================="

# ============================================================ #
#  Step 0 — align to hg38
# ============================================================ #
if [ "$MODE" == "ricopili" ]; then
  ALIGN_SCRIPT="${SCRIPT_DIR}/0_align_1KG_hg38.sh"
  ALIGN_ARGS="${INPUT} ${OUTNAME}"
else
  ALIGN_SCRIPT="${SCRIPT_DIR}/0_align_1KG_hg38_standalone.sh"
  ALIGN_ARGS="${INPUT} ${OUTNAME} ${BUILD}"
fi

JOB0=$(sbatch --parsable "${ALIGN_SCRIPT}" ${ALIGN_ARGS})
echo "Submitted step 0 (align):      job ${JOB0}"

# ============================================================ #
#  Step 1 — ADMIXTURE + ancestry inference
# ============================================================ #
JOB1=$(sbatch --parsable \
  --dependency=afterok:${JOB0} \
  "${SCRIPT_DIR}/1_run_admixture.sh" \
  "${OUTNAME}.hg38.ch.fl" "${OUTNAME}")
echo "Submitted step 1 (admixture):  job ${JOB1} (after ${JOB0})"

# ============================================================ #
#  Step 2 — PCA
# ============================================================ #
if [ "$MODE" == "ricopili" ]; then
  PCA_SCRIPT="${SCRIPT_DIR}/2_pca.sh"
else
  PCA_SCRIPT="${SCRIPT_DIR}/2_pca_standalone.sh"
fi

JOB2=$(sbatch --parsable \
  --dependency=afterok:${JOB1} \
  "${PCA_SCRIPT}" \
  "1kg_${OUTNAME}" "${OUTNAME}")
echo "Submitted step 2 (PCA):        job ${JOB2} (after ${JOB1})"

# ============================================================ #
#  Step 3 — plots (inline SLURM script to resolve PCA file at runtime)
# ============================================================ #
PLOTS_WRAPPER="${SCRIPT_DIR}/tmp/run_plots_${OUTNAME}.sh"
cat > "${PLOTS_WRAPPER}" << EOF
#!/bin/bash
#SBATCH --job-name=plots_${OUTNAME}
#SBATCH --output=tmp/plots_%j.out
#SBATCH --error=tmp/plots_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G

module load 2024
module load Anaconda3/2024.06-1
source /sw/arch/RHEL9/EB_production/2024/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh
conda activate admix_r

# Locate PCA output — ricopili: .menv.mds_cov, standalone: .eigenvec
pca_file=\$(ls ${OUTNAME}*.menv.mds_cov ${OUTNAME}.eigenvec 2>/dev/null | head -1 || true)
if [ -z "\$pca_file" ]; then
  echo "Error: no PCA output file found for ${OUTNAME}"
  exit 1
fi

Rscript ${SCRIPT_DIR}/3_plots.R \$(pwd) 1kg_${OUTNAME}.geno.05.pruned "\$pca_file"
EOF
chmod +x "${PLOTS_WRAPPER}"

JOB3=$(sbatch --parsable \
  --dependency=afterok:${JOB2} \
  "${PLOTS_WRAPPER}")
echo "Submitted step 3 (plots):      job ${JOB3} (after ${JOB2})"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f tmp/ancestry1KG_${JOB1}.out"
