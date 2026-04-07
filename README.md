# ancestry1KG

Ancestry inference pipeline using the 1000 Genomes high-coverage reference and supervised [ADMIXTURE](https://dalexander.github.io/admixture/).  
Assigns individuals to one of 26 sub-populations (5 super-populations: AFR, AMR, EAS, EUR, SAS) and produces PCA-based visualisations.

---

## Prerequisites

### Reference files

The pipeline requires the 1KG high-coverage reference PLINK files in the **working directory** (symlinks are fine):

```
1KG_high_coverage_20130606_g1k_3202.merged.bed
1KG_high_coverage_20130606_g1k_3202.merged.bim
1KG_high_coverage_20130606_g1k_3202.merged.fam
```

> **SURFsnellius users:** files are at `/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/pgcdrc/scz/working/1KG_pca`

Based on 30x Illumina NovaSeq sequencing of 3,202 samples from the 1000 Genomes Project phase 3 sample set, along with 698 family members completing trios. Generated at the New York Genome Center (NHGRI Grant 3UM1HG008901). [AnVIL dataset](https://explore.anvilproject.org/datasets/63d97e58-a825-4b93-a9a7-0d9d165dc826)

### System tools

| Tool | Version tested | Source |
|------|---------------|--------|
| PLINK 1.9 | 1.90b | [cog-genomics.org/plink](https://www.cog-genomics.org/plink/) or `conda install -c bioconda plink` |
| ADMIXTURE | 1.3.0 | [dalexander.github.io/admixture](https://dalexander.github.io/admixture/) |
| R | ≥ 4.1 | `conda install -c conda-forge r-base` |
| liftOver | — | [UCSC](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver) — **standalone workflow only** |
| RICOPILI | — | [RICOPILI wiki](https://sites.google.com/a/broadinstitute.org/ricopili/) — **RICOPILI workflow only** |

### liftOver chain file (standalone workflow only)

Required for hg19 → hg38 conversion. Set the paths in `0_align_1KG_hg38_standalone.sh`:

```bash
LIFTOVER_BIN="/path/to/liftOver"
CHAIN_FILE="/path/to/hg19ToHg38.over.chain.gz"
```

Chain file available from UCSC: `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz`

### R environment

```bash
conda env create -f environment.yml
conda activate admix_r
```

R packages: `dplyr`, `ggplot2`, `tidyr`, `patchwork`, `data.table`, `scales`

---

## Workflows

Two parallel workflows are provided. Steps 1, 1b, and 3 are **identical** in both. Only the alignment (step 0) and PCA (step 2) differ.

### Workflow A — with RICOPILI

```
0_align_1KG_hg38.sh          (impute_dirsub --onlyalign)
        │
        ▼
1_run_admixture.sh            (merge + filter + prune + ADMIXTURE K=26)
        │   └── 1b_ancestry_inference.R  (auto-called)
        ▼
2_pca.sh                      (pcaer)
        ▼
3_plots.R                     (barplot + PCA panels)
```

### Workflow B — standalone (no RICOPILI)

```
0_align_1KG_hg38_standalone.sh   (liftOver + 0b_harmonize_alleles.R)
        │
        ▼
1_run_admixture.sh               (merge + filter + prune + ADMIXTURE K=26)
        │   └── 1b_ancestry_inference.R  (auto-called)
        ▼
2_pca_standalone.sh              (plink --pca)
        ▼
3_plots.R                        (barplot + PCA panels)
```

You can also run the full chain with the SLURM wrapper (see [Pipeline wrapper](#pipeline-wrapper)).

---

## Usage

### Step 0 — align to hg38

If your data are already hg38, skip this step and use your existing PLINK files directly in step 1.

**RICOPILI:**
```bash
sbatch 0_align_1KG_hg38.sh <plink_prefix> <outname>
# output: pi_sub/{plink_prefix}.hg38.ch.fl.bed/bim/fam, symlinked to working dir
```

**Standalone:**
```bash
sbatch 0_align_1KG_hg38_standalone.sh <plink_prefix> <outname>
# calls 0b_harmonize_alleles.R automatically
# output: {outname}.hg38.ch.fl.bed/bim/fam
```

### Step 1 — merge, filter, ADMIXTURE

```bash
sbatch 1_run_admixture.sh <outname>.hg38.ch.fl <output_prefix>
# output: 1kg_{output_prefix}.geno.05.pruned.Q.annotated
```

Steps performed automatically:
1. Merge with 1KG reference (`plink --bmerge`)
2. Missingness filter (`--geno 0.05`)
3. LD pruning (`--indep-pairwise 50 10 0.1`, `--mind 0.1`)
4. Create `.pop` file (population codes for 1KG; `-` for target samples)
5. Supervised ADMIXTURE K=26
6. Ancestry inference and annotation via `1b_ancestry_inference.R`

### Step 2 — PCA

**RICOPILI:**
```bash
sbatch 2_pca.sh 1kg_<output_prefix> <pca_outname>
# output: {pca_outname}.menv.mds_cov (copied to working dir - hopefully)
```

**Standalone:**
```bash
sbatch 2_pca_standalone.sh 1kg_<output_prefix> <pca_outname>
# output: {pca_outname}.eigenvec  (columns: FID IID C1-C10)
```

### Step 3 — plots

```bash
Rscript 3_plots.R <outname> <working_dir> <pca_file>
# e.g.:
Rscript 3_plots.R clz2a . clz2a.menv.mds_cov
Rscript 3_plots.R clz2a . clz2a.eigenvec
```

Outputs (in working dir):
- `1kg_{outname}.geno.05.pruned_barplot.jpeg` — ancestry proportions for inferred samples
- `1kg_{outname}.geno.05.pruned_pca_plot_PC1_vs_C{2,3,4,5}_panels.jpeg` — side-by-side 1KG / Inferred panels

---

## Pipeline wrapper

`run_pipeline.sh` submits all steps as a SLURM dependency chain:

```bash
# Check dependencies before submitting:
bash run_pipeline.sh --check

# Run RICOPILI workflow:
bash run_pipeline.sh --mode ricopili <plink_prefix> <outname>

# Run standalone workflow:
bash run_pipeline.sh --mode standalone <plink_prefix> <outname>
```

---

## 1000 Genomes population reference

| Code | Population | Superpopulation |
|------|-----------|----------------|
| AFR_ACB | African Caribbeans in Barbados | AFR |
| AFR_ASW | African Ancestry in Southwest US | AFR |
| AFR_ESN | Esan in Nigeria | AFR |
| AFR_GWD | Gambian in Western Divisions | AFR |
| AFR_LWK | Luhya in Webuye | AFR |
| AFR_MSL | Mende in Sierra Leone | AFR |
| AFR_YRI | Yoruba in Ibadan | AFR |
| AMR_CLM | Colombians in Medellín | AMR |
| AMR_MXL | Mexican Ancestry in Los Angeles | AMR |
| AMR_PEL | Peruvians in Lima | AMR |
| AMR_PUR | Puerto Ricans in Puerto Rico | AMR |
| EAS_CDX | Chinese Dai in Xishuangbanna | EAS |
| EAS_CHB | Han Chinese in Beijing | EAS |
| EAS_CHS | Southern Han Chinese | EAS |
| EAS_JPT | Japanese in Tokyo | EAS |
| EAS_KHV | Kinh in Ho Chi Minh City | EAS |
| EUR_CEU | Utah Residents (CEPH) with NW European ancestry | EUR |
| EUR_FIN | Finnish in Finland | EUR |
| EUR_GBR | British in England and Scotland | EUR |
| EUR_IBS | Iberian Population in Spain | EUR |
| EUR_TSI | Toscani in Italy | EUR |
| SAS_BEB | Bengali in Bangladesh | SAS |
| SAS_GIH | Gujarati Indian in Houston | SAS |
| SAS_ITU | Indian Telugu in the UK | SAS |
| SAS_PJL | Punjabi in Lahore | SAS |
| SAS_STU | Sri Lankan Tamil in the UK | SAS |
