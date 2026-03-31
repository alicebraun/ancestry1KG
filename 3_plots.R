# Load necessary library
library(dplyr)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
submit_dir <- args[1]
prefix <- args[2]
cat("Working directory:", submit_dir, "\n")
cat("Prefix:", prefix, "\n")
setwd(submit_dir)


# Load data
df <- read.table(paste0(prefix,".Q.annotated"),
                 header = TRUE, 
                 sep = "\t", 
                 quote = "", 
                 stringsAsFactors = FALSE)


# Define the components you want to include in the plot
included_components <- c("Q_AFR", "Q_AMR", "Q_EAS", "Q_EUR", "Q_SAS")

# Ensure those columns exist (optional but safe)
included_components <- intersect(included_components, names(df))

# Reshape only those ancestry components for plotting
df_long <- df %>%
  select(FID, IID, all_of(included_components)) %>%
  pivot_longer(cols = all_of(included_components),
               names_to = "Ancestry_Component",
               values_to = "Proportion")

#palette used in the bar plot
fill_colors <- scale_fill_discrete()$palette(5)

# Plot
df_long %>% filter(grepl("cas|con", FID)) %>% ggplot(aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
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
ggsave(paste0(prefix,"_barplot.jpeg"), plot = p1, width = 12, height = 6, dpi = 600, device = "jpeg")


# Load the pca covariate file
mds_file <- args[3]

# Check if the file exists
if (!file.exists(mds_file)) {
  stop(paste("The specified PCA file does not exist:", mds_file))
}

# Read in the PCA/PCoA covariate data
pca_data <- read.table(mds_file, header = TRUE, stringsAsFactors = FALSE)

cat("Successfully loaded PCA covariate file with", nrow(pca_data), "rows and", ncol(pca_data), "columns.\n")


# Merge the PCs and annotad Q.file on FID and IID
merged_data <- pca_data %>%
  inner_join(df, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)  # Extract prefix
head(merged_data)
# Create the scatter plot
#merged_data %>% filter(Ancestry_Inf != "MIX") %>% ggplot(aes(x = C1, y = C2, color = Max_Ancestry, shape = Ancestry_Prefix)) +
#merged_data %>% ggplot(aes(x = C1, y = C2, color = Max_Ancestry, shape = Ancestry_Prefix)) +
#  geom_point(size = 2, alpha = 0.8) +
#  theme_minimal(base_size = 12) +
#  labs(title = "PC1 vs PC2 with inferred Ancestry", x = "PC1 (C1)", y = "PC2 (C2)", color = "Ancestry_Inf", shape = "Ancestry Prefix") +
#    scale_shape_manual(values = c("AFR" = 16, "AMR" = 17, "EAS" = 8, "EUR" = 19, "SAS" = 20, "MIX" = 5)) +  # Shape per prefix
#  theme(legend.position = "right",
#        plot.background = element_rect(fill = "white", color = "white"),
#        panel.background = element_rect(fill = "white", color = "white"))
# Optional: Save the plot as a JPEG image
#ggsave(paste0(prefix,"_pcaplot.jpeg"), width = 8, height = 6, dpi = 600, device = "jpeg")


# Merge the PCs and annotated Q.file on FID and IID
merged_data <- pca_data %>%
  inner_join(df, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
#merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)
merged_data$Ancestry_Prefix <- sub("^([A-Z]+)_.*", "\\1", merged_data$Ancestry_Inf)
unique(merged_data$Ancestry_Prefix)
# Create the scatter plot

# Assign shape group
merged_data$shape <- ifelse(grepl("mis", merged_data$FID), "1KG", "inferred")


# PCs to loop over
pcs <- c("C2", "C3", "C4", "C5")

for (pc in pcs) {
  # ----- Non-MIX plot ----- #
  p_nonmix <- ggplot() +
    geom_point(
      data = merged_data %>% filter(shape == "1KG"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1,
      alpha = 0.5
    ) +
    geom_point(
      data = merged_data %>% filter(Ancestry_Inf != "MIX", shape == "inferred"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1
    ) +
    scale_shape_manual(values = c("1KG" = 1, "inferred" = 6)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PC1 vs ", pc, " with inferred Ancestry (No MIX)"),
      x = "PC1 (C1)", y = paste0(pc, " (", pc, ")"),
      color = "Ancestry Prefix",
      shape = "Reference"
    ) +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  ggsave(paste0(prefix, "_pca_plot_PC1_vs_", pc, "_nomix.jpeg"), p_nonmix, width = 8, height = 6, dpi = 600, device = "jpeg")

  # ----- MIX plot ----- #
  p_mix <- ggplot() +
    geom_point(
      data = merged_data %>% filter(shape == "1KG"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1,
      alpha = 0.3
    ) +
    geom_point(
      data = merged_data %>% filter(Ancestry_Inf == "MIX", shape == "inferred"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1
    ) +
    scale_shape_manual(values = c("1KG" = 1, "inferred" = 6)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PC1 vs ", pc, " with inferred Ancestry (MIX only)"),
      x = "PC1 (C1)", y = paste0(pc, " (", pc, ")"),
      color = "Ancestry Prefix",
      shape = "Reference"
    ) +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  ggsave(paste0(prefix, "_pca_plot_PC1_vs_", pc, "_mix.jpeg"), p_mix, width = 8, height = 6, dpi = 600, device = "jpeg")
}


summary(merged_data[, c("C1", "C2", "C3", "C4", "C5")])
merged_data %>% filter(is.na(C1) | is.na(C3) | is.na(C4) | is.na(C5)) %>% nrow()


# Load necessary library
library(dplyr)
library(ggplot2)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
submit_dir <- args[1]
prefix <- args[2]
cat("Working directory:", submit_dir, "\n")
cat("Prefix:", prefix, "\n")
setwd(submit_dir)


# Load data
df <- read.table(paste0(prefix,".Q.annotated"),
                 header = TRUE, 
                 sep = "\t", 
                 quote = "", 
                 stringsAsFactors = FALSE)


# Define the components you want to include in the plot
included_components <- c("Q_AFR", "Q_AMR", "Q_EAS", "Q_EUR", "Q_SAS")

# Ensure those columns exist (optional but safe)
included_components <- intersect(included_components, names(df))

# Reshape only those ancestry components for plotting
df_long <- df %>%
  select(FID, IID, all_of(included_components)) %>%
  pivot_longer(cols = all_of(included_components),
               names_to = "Ancestry_Component",
               values_to = "Proportion")

#palette used in the bar plot
fill_colors <- scale_fill_discrete()$palette(5)

# Plot
df_long %>% filter(grepl("cas|con", FID)) %>% ggplot(aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
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
ggsave(paste0(prefix,"_barplot.jpeg"), plot = p1, width = 12, height = 6, dpi = 600, device = "jpeg")


# Load the pca covariate file
mds_file <- args[3]

# Check if the file exists
if (!file.exists(mds_file)) {
  stop(paste("The specified PCA file does not exist:", mds_file))
}

# Read in the PCA/PCoA covariate data
pca_data <- read.table(mds_file, header = TRUE, stringsAsFactors = FALSE)

cat("Successfully loaded PCA covariate file with", nrow(pca_data), "rows and", ncol(pca_data), "columns.\n")


# Merge the PCs and annotad Q.file on FID and IID
merged_data <- pca_data %>%
  inner_join(df, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)  # Extract prefix
head(merged_data)
# Create the scatter plot
#merged_data %>% filter(Ancestry_Inf != "MIX") %>% ggplot(aes(x = C1, y = C2, color = Max_Ancestry, shape = Ancestry_Prefix)) +
#merged_data %>% ggplot(aes(x = C1, y = C2, color = Max_Ancestry, shape = Ancestry_Prefix)) +
#  geom_point(size = 2, alpha = 0.8) +
#  theme_minimal(base_size = 12) +
#  labs(title = "PC1 vs PC2 with inferred Ancestry", x = "PC1 (C1)", y = "PC2 (C2)", color = "Ancestry_Inf", shape = "Ancestry Prefix") +
#    scale_shape_manual(values = c("AFR" = 16, "AMR" = 17, "EAS" = 8, "EUR" = 19, "SAS" = 20, "MIX" = 5)) +  # Shape per prefix
#  theme(legend.position = "right",
#        plot.background = element_rect(fill = "white", color = "white"),
#        panel.background = element_rect(fill = "white", color = "white"))
# Optional: Save the plot as a JPEG image
#ggsave(paste0(prefix,"_pcaplot.jpeg"), width = 8, height = 6, dpi = 600, device = "jpeg")


# Merge the PCs and annotated Q.file on FID and IID
merged_data <- pca_data %>%
  inner_join(df, by = c("FID", "IID"))

# Create a new column for the prefix (e.g., AFR, EUR)
#merged_data$Ancestry_Prefix <- gsub("Q_([A-Za-z]+)_.*", "\\1", merged_data$Ancestry_Inf)
merged_data$Ancestry_Prefix <- sub("^([A-Z]+)_.*", "\\1", merged_data$Ancestry_Inf)
unique(merged_data$Ancestry_Prefix)
# Create the scatter plot

# Assign shape group
merged_data$shape <- ifelse(grepl("mis", merged_data$FID), "1KG", "inferred")


# PCs to loop over
pcs <- c("C2", "C3", "C4", "C5")

for (pc in pcs) {
  # ----- Non-MIX plot ----- #
  p_nonmix <- ggplot() +
    geom_point(
      data = merged_data %>% filter(shape == "1KG"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1,
      alpha = 0.5
    ) +
    geom_point(
      data = merged_data %>% filter(Ancestry_Inf != "MIX", shape == "inferred"),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix", shape = "shape"),
      size = 1
    ) +
    scale_shape_manual(values = c("1KG" = 1, "inferred" = 6)) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PC1 vs ", pc, " with inferred Ancestry (No MIX)"),
      x = "PC1 (C1)", y = paste0(pc, " (", pc, ")"),
      color = "Ancestry Prefix",
      shape = "Reference"
    ) +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  ggsave(paste0(prefix, "_pca_plot_PC1_vs_", pc, "_nomix.jpeg"), p_nonmix, width = 8, height = 6, dpi = 600, device = "jpeg")
}
  for (pc in pcs) {
  p_mis_colored <- ggplot() +
    # Gray dots for all non-"mis" individuals
    geom_point(
      data = merged_data %>% filter(!grepl("mis", FID)),
      aes_string(x = "C1", y = pc),
      color = "gray70",
      shape = 16,
      size = 1,
      alpha = 0.5
    ) +
    # Colored triangles for "mis" individuals
    geom_point(
      data = merged_data %>% filter(grepl("mis", FID)),
      aes_string(x = "C1", y = pc, color = "Ancestry_Prefix"),
      shape = 17,
      size = 1,
        alpha = 0.5
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste0("PC1 vs ", pc, ": 1KG individuals colored by ancestry"),
      x = "PC1 (C1)", y = paste0(pc, " (", pc, ")"),
      color = "Ancestry Prefix"
    ) +
    theme(
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  ggsave(
    paste0(prefix, "_pca_plot_PC1_vs_", pc, "_mis_colored.jpeg"),
    p_mis_colored,
    width = 8, height = 6,
    dpi = 600,
    device = "jpeg"
  )
}