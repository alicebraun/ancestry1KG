# ancestry1KG

Inferring ancestry using the 1000 Genomes high coverage reference and admixture (https://dalexander.github.io/admixture/)

```
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xvf admixture_linux-1.3.0.tar.gz
export PATH="$HOME/admixture/dist/admixture_linux-1.3.0:$PATH"
```

# 1000 Genomes population overview
| Code     | Population Name                          | Superpopulation        |
|----------|-------------------------------------------|------------------------|
| AFR_ACB  | African Caribbeans in Barbados            | AFR (African)          |
| AFR_ASW  | African Ancestry in Southwest US          | AFR (African)          |
| AFR_ESN  | Esan in Nigeria                           | AFR (African)          |
| AFR_GWD  | Gambian in Western Divisions              | AFR (African)          |
| AFR_LWK  | Luhya in Webuye                           | AFR (African)          |
| AFR_MSL  | Mende in Sierra Leone                     | AFR (African)          |
| AFR_YRI  | Yoruba in Ibadan                          | AFR (African)          |
| AMR_CLM  | Colombians in Medellín                    | AMR (Admixed American) |
| AMR_MXL  | Mexican Ancestry in Los Angeles           | AMR (Admixed American) |
| AMR_PEL  | Peruvians in Lima                         | AMR (Admixed American) |
| AMR_PUR  | Puerto Ricans in Puerto Rico              | AMR (Admixed American) |
| EAS_CDX  | Chinese Dai in Xishuangbanna              | EAS (East Asian)       |
| EAS_CHB  | Han Chinese in Beijing                    | EAS (East Asian)       |
| EAS_CHS  | Southern Han Chinese                      | EAS (East Asian)       |
| EAS_JPT  | Japanese in Tokyo                         | EAS (East Asian)       |
| EAS_KHV  | Kinh in Ho Chi Minh City                  | EAS (East Asian)       |
| EUR_CEU  | Utah Residents (CEPH) with NW European ancestry | EUR (European)   |
| EUR_FIN  | Finnish in Finland                        | EUR (European)         |
| EUR_GBR  | British in England and Scotland           | EUR (European)         |
| EUR_IBS  | Iberian Population in Spain               | EUR (European)         |
| EUR_TSI  | Toscani in Italy                          | EUR (European)         |
| SAS_BEB  | Bengali in Bangladesh                     | SAS (South Asian)      |
| SAS_GIH  | Gujarati Indian in Houston                | SAS (South Asian)      |
| SAS_ITU  | Indian Telugu in the UK                   | SAS (South Asian)      |
| SAS_PJL  | Punjabi in Lahore                         | SAS (South Asian)      |
| SAS_STU  | Sri Lankan Tamil in the UK                | SAS (South Asian)      |

## Step 1: Merge your plink files with the 1KG reference

Note that the file needs to be QCed and on hg38

```bash
plink --bfile 1KG_high_coverage_20130606_g1k_3202.merged --bmerge $prefix --make-bed --out 1kg_$prefix.merged
plink --bfile 1kg_$prefix.merged --missing --out 1kg_$prefix.merged
echo "SNPs with a missing rate > 5%: $(awk '$5 > 0.05' 1kg_$prefix.merged.lmiss | wc -l)"
echo "excluding SNPs with a missing rate > 5%..."
plink --bfile 1kg_$prefix.merged --geno 0.05 --make-bed --out 1kg_$prefix.geno.05.merged
```

## Step 2: Pruning

The admixture manual recommends to prune the input files. `<br>`
The first command targets for removal each SNP that has an R2 value of greater than 01 with any other SNP within a 50-SNP sliding window (advanced by 10 SNPs each time).

```bash
plink --bfile 1kg_$prefix.geno.05.merged --indep-pairwise 50 10 0.1 
plink --bfile 1kg_$prefix.geno.05.merged --extract plink.prune.in --make-bed --out 1kg_$prefix.geno.05.merged.pruned
```

## Step 3: Generate .pop file

In order to run supervised ancestry inference a .pop file indicating the known ancestry is needed.

```bash
# count ancestries in the reference
echo "Ancestral populations in the reference file: $(awk '{print $1}' 1kg_$prefix.geno.05.merged.pruned.fam | grep -v 'con' | grep -v 'cas' | sort | uniq -c)"
awk '{
    # Exclude FIDs containing "con" or "cas" and output "-"
    if ($1 ~ /con/ || $1 ~ /cas/) {
        print "-";
    } else {
        # Otherwise, output the ancestry name based on the FID
        # In this case, output the FID or whatever ancestry you expect from FID
        # Extract ancestry from the FID (Assuming FID format: "AFR_ACB", "EUR_GBR", etc.)
        print $1;
    }
}' 1kg_$prefix.geno.05.merged.pruned.fam > 1kg_$prefix.geno.05.merged.pruned.pop
```

## Step 4: Run ADMIXTURE

```bash
# run admixture [options] inputFile K
admixture --supervised --seed 666 -C 10  -j16 1kg_$prefix.geno.05.merged.pruned.bed 26
```

## Step 5: Mapping

```R
# Load necessary library
library(dplyr)
library(tidyr)

# Define your prefix (from an environment variable or hardcoded)
prefix <- "your_prefix]"

# Build base path
base <- paste0("1kg_scz_", prefix, ".geno.05.merged.pruned")

# Build file names
fam_file <- paste0(base, ".fam")
q_file   <- paste0(base, ".26.Q")
pop_file <- paste0(base, ".pop")

# Read each file
fam <- read.table(fam_file, header = FALSE)
q   <- read.table(q_file, header = FALSE)
pop <- read.table(pop_file, header = FALSE)

# Combine them as a space-separated string (like for a shell command or input list)
df <- cbind(fam, q, pop)

# Define and assign all column names
base_cols <- c("FID", "IID", "PID", "MID", "Sex", "Pheno", )
q_cols <- paste0("Q", 1:num_q)
final_col <- "Ancestry_Pheno"
all_cols <- c(base_cols, q_cols, final_col)
colnames(df) <- all_cols

# Compute average Q values per ancestry
avg_q <- df %>%
  group_by(.data[[final_col]]) %>%
  summarise(across(all_of(q_cols), mean), .groups = "drop")

print(avg_q)

# Determine which ancestry corresponds to each Q column
q_labels <- colnames(avg_q)[2:27]  # Dynamically use num_q to get Q columns
ancestry_labels <- apply(avg_q[, q_labels], 2, function(x) {
  avg_q[[final_col]][which.max(x)]  # Find the ancestry with the highest value for each Q
})

new_colnames <- paste0("Q_", ancestry_labels)
# Combine all labels
colnames(df) <- c(base_cols, new_colnames, final_col)

df <- df %>%
  # Step 1: Find the maximum Q-value for each row
  mutate(SortIndex = apply(select(., starts_with("Q_")), 1, max)) %>%
  
  # Step 2: Assign Ancestry_Inf based on the column associated with the maximum Q-value
  mutate(
    Ancestry_Inf = case_when(
      SortIndex > 0.50 ~ names(select(., starts_with("Q_")))[apply(select(., starts_with("Q_")), 1, which.max)],
      TRUE ~ "MIX"  # If no ancestry exceeds 0.50, mark as "mixed"
    )
  ) %>%
  
  # Step 3: Sort the data by Ancestry_Inf and SortIndex in descending order
  arrange(Ancestry_Inf, desc(SortIndex)) %>%
  
  # Step 4: Preserve the new order of IID for plotting
  mutate(IID = factor(IID, levels = IID))

# Save the annotated output
write.table(df, "annotated_Q_file.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
```

## Step 6: Barplot (optional)
```R
library(ggplot2)
df_long <- df %>%
  select(FID, IID, starts_with("Q_")) %>%
  pivot_longer(cols = starts_with("Q_"), names_to = "Ancestry_Component", values_to = "Proportion")

# Plot
df_long %>% filter(grepl("pash1", FID)) %>% ggplot(aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  geom_bar(stat = "identity", width = 1) +
  labs(title = "Ancestry Proportions per Individual", x = "Individuals", y = "Proportion") -> p1
# Save the plot
ggsave("barplot.jpeg", plot = p1, width = 12, height = 6, dpi = 600, device = "jpeg")

```
## Step 7: Map to principal components (optinal)
```R
library(ggplot2)
pca_data <- read.table("your_pca_results.mds", header = TRUE, stringsAsFactors = FALSE) # e.g. from RICOPILI or EIGENSTRAT
annotated_q_data <- read.table("annotated_Q_file.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge the two datasets on FID and IID
merged_data <- pca_data %>%
  inner_join(annotated_q_data, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)  # Extract ancestry prefix
head(merged_data)
# Create the scatter plot
ggplot(merged_data, aes(x = C1, y = C2, color = Ancestry_Inf, shape = Ancestry_Prefix)) +
  geom_point(size = 2) +
  theme_minimal(base_size = 12) +
  labs(title = "PC1 vs PC2 with inferred Ancestry", x = "PC1 (C1)", y = "PC2 (C2)", color = "Ancestry_Inf", shape = "Ancestry Prefix") +
    scale_shape_manual(values = c("AFR" = 16, "AMR" = 17, "EAS" = 8, "EUR" = 19, "SAS" = 20, "MIX" = 5))  # Shape per prefix
# Save the plot
ggsave("pc1pc2_scatterplot.jpeg", width = 8, height = 6, dpi = 600, device = "jpeg")