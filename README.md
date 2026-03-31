# ancestry1KG

Inferring ancestry using the 1000 Genomes high coverage reference and admixture (https://dalexander.github.io/admixture/)

## Dependencies

### ADMIXTURE

```
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xvf admixture_linux-1.3.0.tar.gz
export PATH="$HOME/admixture/dist/admixture_linux-1.3.0:$PATH"
```

### Reference file

Prerequsities include a reference file, e.g. the high coverage 1KG reference file `1KG_high_coverage_20130606_g1k_3202.merged` containing 26 subpopulations and 5 superpopulations. `<br>`
Your reference file should contain information about the ancestry of individuals e.g. in the FID, which is needed to create the `.pop` file for supervised ADMIXTURE ancestry inference. `<br>`
For more details check out the [ADMIXTURE manual](https://dalexander.github.io/admixture/admixture-manual.pdf).

#### 1000 Genomes population overview

| Code    | Population Name                                 | Superpopulation        |
| ------- | ----------------------------------------------- | ---------------------- |
| AFR_ACB | African Caribbeans in Barbados                  | AFR (African)          |
| AFR_ASW | African Ancestry in Southwest US                | AFR (African)          |
| AFR_ESN | Esan in Nigeria                                 | AFR (African)          |
| AFR_GWD | Gambian in Western Divisions                    | AFR (African)          |
| AFR_LWK | Luhya in Webuye                                 | AFR (African)          |
| AFR_MSL | Mende in Sierra Leone                           | AFR (African)          |
| AFR_YRI | Yoruba in Ibadan                                | AFR (African)          |
| AMR_CLM | Colombians in Medellín                         | AMR (Admixed American) |
| AMR_MXL | Mexican Ancestry in Los Angeles                 | AMR (Admixed American) |
| AMR_PEL | Peruvians in Lima                               | AMR (Admixed American) |
| AMR_PUR | Puerto Ricans in Puerto Rico                    | AMR (Admixed American) |
| EAS_CDX | Chinese Dai in Xishuangbanna                    | EAS (East Asian)       |
| EAS_CHB | Han Chinese in Beijing                          | EAS (East Asian)       |
| EAS_CHS | Southern Han Chinese                            | EAS (East Asian)       |
| EAS_JPT | Japanese in Tokyo                               | EAS (East Asian)       |
| EAS_KHV | Kinh in Ho Chi Minh City                        | EAS (East Asian)       |
| EUR_CEU | Utah Residents (CEPH) with NW European ancestry | EUR (European)         |
| EUR_FIN | Finnish in Finland                              | EUR (European)         |
| EUR_GBR | British in England and Scotland                 | EUR (European)         |
| EUR_IBS | Iberian Population in Spain                     | EUR (European)         |
| EUR_TSI | Toscani in Italy                                | EUR (European)         |
| SAS_BEB | Bengali in Bangladesh                           | SAS (South Asian)      |
| SAS_GIH | Gujarati Indian in Houston                      | SAS (South Asian)      |
| SAS_ITU | Indian Telugu in the UK                         | SAS (South Asian)      |
| SAS_PJL | Punjabi in Lahore                               | SAS (South Asian)      |
| SAS_STU | Sri Lankan Tamil in the UK                      | SAS (South Asian)      |

## Prerequisites
### PLINK 1.9
To run the script you need PLINK 1.9 installed as a module on your HPC or manually and added to your bashrc.
You can download PLINK 1.9 via https://www.cog-genomics.org/plink/ or conda/mamba: `conda install bioconda::plink`

### R
Base R along with the `dplyr` and `ggplot2` libraries.

## Step 0: Prep PLINK files 
```bash
# replace PLINK prefix and output filename
 sbatch 0_align_1KG_hg38.sh <your_bfile> <your_outname> 
```
If you have [RICOPILI]()(https://github.com/Ripkelab/ricopili) installed, run `0_align_1KG_hg38.sh`. This script starts a module that aligns to hg38 based on a reference (in this case `1KG_high_coverage_mm4_jan2025/`), aligns snp names based on positions and alleles, standardize indel coding, align alleles, and resolves strandflips.

## Step 1: Merge with 1KG reference
```bash
# replace PLINK prefix and output filename
 sbatch 1_run_admixture.sh <your_bfile>.hg38.ch.fl <your_outname> 
```
This step performs merging, filtering for missingness (`--geno 0.05`), and pruning (`--indep-pairwise 50 10 0.1`) using PLINK 1.9.

## Step 2: Run Admixture
After merging and filtering ADMIXTURE is run using the supervised module.


## Step 3: Mapping ancestry in R

## Step 4: run PCA

## Step 5: Plots
```

```
