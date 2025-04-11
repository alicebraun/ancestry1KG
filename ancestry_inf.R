# Load necessary library
library(dplyr)
library(ggplot2)
library(tidyr)

# Define your prefix (from an environment variable or hardcoded)
prefix <- "pash1"
k <- "26"
# Build base path
base <- paste0("1kg_scz_", prefix, ".geno.05.merged.pruned")

# Build file names
fam_file <- paste0(base, ".fam")
q_file   <- paste0(base, "."k".Q")
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

cat("Annotated file saved as 'annotated_Q_file.txt'\n")


# Reshape for plotting
df_long <- df %>%
  select(FID, IID, starts_with("Q_")) %>%
  pivot_longer(cols = starts_with("Q_"), names_to = "Ancestry_Component", values_to = "Proportion")

palette used in the bar plot
fill_colors <- scale_fill_discrete()$palette(26)

# Plot
df_long %>% filter(grepl("pash1", FID)) %>% ggplot(aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = "white"),  # Set white background
    panel.background = element_rect(fill = "white", color = "white")  # Ensure panel is white as well
  ) + scale_fill_manual(values = fill_colors)  +
  labs(title = "Ancestry Proportions per Individual", x = "Individuals", y = "Proportion") -> p1

# Save the plot
ggsave("barplot.jpeg", plot = p1, width = 12, height = 6, dpi = 600, device = "jpeg")
cat("✔ Barplot saved to:", output_plot, "\n")


# Step 1: Load the datasets
pca_data <- read.table("pca1_scz_pash1_mix_ab-qc1.menv.mds", header = TRUE, stringsAsFactors = FALSE)
annotated_q_data <- read.table("annotated_Q_file.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge the two datasets on FID and IID
merged_data <- pca_data %>%
  inner_join(annotated_q_data, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)  # Extract prefix
head(merged_data)
# Create the scatter plot
ggplot(merged_data, aes(x = C1, y = C2, color = Ancestry_Inf, shape = Ancestry_Prefix)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "PC1 vs PC2 with inferred Ancestry", x = "PC1 (C1)", y = "PC2 (C2)", color = "Ancestry_Inf", shape = "Ancestry Prefix") +
    scale_shape_manual(values = c("AFR" = 16, "AMR" = 17, "EAS" = 8, "EUR" = 19, "SAS" = 20, "MIX" = 5)) +  # Shape per prefix
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white", color = "white"))

# Optional: Save the plot as a JPEG image
ggsave("pc1pc2_scatterplot.jpeg", width = 8, height = 6, dpi = 600, device = "jpeg")
cat("✔ Scatterplot saved as 'PC1_vs_PC2_scatterplot.jpeg'\n")