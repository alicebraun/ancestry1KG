# Load necessary library
library(dplyr)
library(ggplot2)
library(tidyr)

# Define your prefix (from an environment variable or hardcoded)
getwd()
setwd("your_path_here")
prefix <- "your_prefix"

df <- read.table(paste0(prefix,"_annotated_Q_file.txt"),
                 header = TRUE, 
                 sep = "\t", 
                 quote = "", 
                 stringsAsFactors = FALSE)

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
ggsave("barplot_test.jpeg", plot = p1, width = 12, height = 6, dpi = 600, device = "jpeg")
cat("✔ Barplot saved to:", output_plot, "\n")


# Step 1: Load the datasets
pca_data <- read.table("pca1_scz_pash1_mix_ab-qc1.menv.mds", header = TRUE, stringsAsFactors = FALSE)

# Merge the two datasets on FID and IID
merged_data <- pca_data %>%
  inner_join(df, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)  # Extract prefix
head(merged_data)
# Create the scatter plot
ggplot(merged_data, aes(x = C1, y = C2, color = Max_Ancestry, shape = Ancestry_Prefix)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "PC1 vs PC2 with inferred Ancestry", x = "PC1 (C1)", y = "PC2 (C2)", color = "Ancestry_Inf", shape = "Ancestry Prefix") +
    scale_shape_manual(values = c("AFR" = 16, "AMR" = 17, "EAS" = 8, "EUR" = 19, "SAS" = 20, "MIX" = 5)) +  # Shape per prefix
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white", color = "white"))

# Optional: Save the plot as a JPEG image
ggsave("pc1pc2_test.jpeg", width = 8, height = 6, dpi = 600, device = "jpeg")
cat("✔ Scatterplot saved as 'PC1_vs_PC2_scatterplot.jpeg'\n")