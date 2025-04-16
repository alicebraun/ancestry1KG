args <- commandArgs(trailingOnly = TRUE)
submit_dir <- args[1]
prefix <- args[2]
cat("Working directory:", submit_dir, "\n")
cat("Prefix:", prefix, "\n")
setwd(submit_dir)

# Load necessary library
library(dplyr)
library(tidyr)

# Build file names
fam_file <- paste0("1kg_",prefix, ".fam")
q_file   <- paste0("1kg_",prefix, ".26.Q")
pop_file <- paste0("1kg_",prefix, ".pop")

# Read each file
fam <- read.table(fam_file, header = FALSE)
q   <- read.table(q_file, header = FALSE)
pop <- read.table(pop_file, header = FALSE)

# Combine them as a space-separated string (like for a shell command or input list)
df <- cbind(fam, q, pop)

# Define and assign all column names
base_cols <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")
q_cols <- paste0("Q", 1:26)
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

# Get all Q_ columns
q_cols <- grep("^Q_", colnames(df), value = TRUE)

# Get superpopulation prefixes
superpops <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Find the max Q-value (SortIndex) and associated column
df <- df %>%
  mutate(
    SortIndex = apply(select(., all_of(q_cols)), 1, max),
    Max_Ancestry = names(select(., all_of(q_cols)))[
      apply(select(., all_of(q_cols)), 1, which.max)
    ]
  )

# If max ancestry < 0.70, check superpop sums
# First compute superpop-level Qs
for (super in superpops) {
  pattern <- paste0("^Q_", super, "_")
  super_cols <- grep(pattern, colnames(df), value = TRUE)
  df[[paste0("Q_", super)]] <- rowSums(df[, super_cols], na.rm = TRUE)
}

# Assign final Ancestry_Inf (either fine-grained or superpopulation)
df <- df %>%
  mutate(
    Ancestry_Inf = case_when(
      SortIndex > 0.70 ~ gsub("^Q_", "", Max_Ancestry),
      Q_AFR >= 0.70 ~ "AFR",
      Q_AMR >= 0.70 ~ "AMR",
      Q_EAS >= 0.70 ~ "EAS",
      Q_EUR >= 0.70 ~ "EUR",
      Q_SAS >= 0.70 ~ "SAS",
      TRUE ~ "MIX"
    )
  ) %>%
  arrange(Ancestry_Inf, desc(SortIndex)) %>%
  mutate(IID = factor(IID, levels = IID))

write.table(df, paste0("1kg_",prefix,".Q.annotated"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


