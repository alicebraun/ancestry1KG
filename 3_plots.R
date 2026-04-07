# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)   # for two-panel PCA plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript 3_plots.R <outname> <submit_dir> <pca_file>\n  e.g. Rscript 3_plots.R clz2a /path/to/anc clz2a.menv.mds")
}

outname    <- args[1]
submit_dir <- args[2]
mds_file   <- args[3]
prefix     <- paste0("1kg_", outname, ".geno.05.pruned")

cat("Working directory:", submit_dir, "\n")
cat("Outname:", outname, "\n")
cat("Q.annotated prefix:", prefix, "\n")
cat("PCA file:", mds_file, "\n")
setwd(submit_dir)


# ============================================================ #
#  Ancestry barplot (inferred samples only)
# ============================================================ #

df <- read.table(paste0(prefix, ".Q.annotated"),
                 header = TRUE, sep = "\t", quote = "",
                 stringsAsFactors = FALSE)

included_components <- intersect(
  c("Q_AFR", "Q_AMR", "Q_EAS", "Q_EUR", "Q_SAS"),
  names(df)
)

df_inferred <- df %>%
  filter(Ancestry_Pheno == "-") %>%
  select(FID, IID, all_of(included_components))

# Sort individuals: primary sort by dominant component, secondary by its proportion
dominant <- apply(df_inferred[, included_components, drop = FALSE], 1, which.max)
dom_prop  <- apply(df_inferred[, included_components, drop = FALSE], 1, max)
iid_order <- df_inferred$IID[order(dominant, -dom_prop)]

df_long <- df %>%
  select(FID, IID, Ancestry_Pheno, all_of(included_components)) %>%
  pivot_longer(cols = all_of(included_components),
               names_to  = "Ancestry_Component",
               values_to = "Proportion") %>%
  mutate(IID = factor(IID, levels = iid_order))

fill_colors <- scales::hue_pal()(length(included_components))

p_bar <- df_long %>%
  filter(Ancestry_Pheno == "-") %>%
  ggplot(aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = fill_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.background  = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white")
  ) +
  labs(title = "Ancestry Proportions per Individual",
       x = "Individuals", y = "Proportion")

ggsave(paste0(prefix, "_barplot.jpeg"),
       plot = p_bar, width = 12, height = 6, dpi = 600, device = "jpeg")


# ============================================================ #
#  Load PCA and merge with annotated Q file
# ============================================================ #

pca_data <- read.table(mds_file, header = TRUE, stringsAsFactors = FALSE)
cat("Loaded PCA file:", nrow(pca_data), "rows,", ncol(pca_data), "columns\n")

# Normalise column names:
#   pcaer/PLINK2 may write "#FID" as first column header
#   pcaer MDS files use MDS1, MDS2, ... — rename to C1, C2, ... for consistency
names(pca_data) <- sub("^X\\.", "", names(pca_data))   # R turns "#FID" into "X.FID"
names(pca_data) <- sub("^#", "", names(pca_data))       # strip any remaining leading #
names(pca_data) <- gsub("^MDS([0-9]+)$", "C\\1", names(pca_data))
names(pca_data) <- gsub("^PC([0-9]+)$",  "C\\1", names(pca_data))
cat("PCA column names:", paste(names(pca_data), collapse = ", "), "\n")
cat("Q.annotated FID examples:", paste(head(df$FID, 3), collapse = ", "), "\n")
cat("PCA FID examples:        ", paste(head(pca_data$FID, 3), collapse = ", "), "\n")

merged_data <- pca_data %>%
  inner_join(df, by = "IID")

if (nrow(merged_data) == 0) {
  stop(paste0(
    "inner_join returned 0 rows — IID values do not match between the ",
    "PCA file and the Q.annotated file.\n",
    "  PCA IID examples: ", paste(head(pca_data$IID, 5), collapse = ", "), "\n",
    "  Q   IID examples: ", paste(head(df$IID, 5), collapse = ", ")
  ))
}

# Ancestry superpopulation prefix (AFR, EUR, EAS, AMR, SAS, MIX)
merged_data$Ancestry_Prefix <- sub("^([A-Z]+)_.*", "\\1", merged_data$Ancestry_Inf)

# Label 1KG reference vs samples being inferred.
# Ancestry_Pheno comes from the .pop file: "-" = sample to infer, anything else = 1KG reference.
merged_data$Group <- ifelse(merged_data$Ancestry_Pheno == "-", "Inferred", "1KG")

cat("Ancestry prefixes found:", paste(unique(merged_data$Ancestry_Prefix), collapse = ", "), "\n")
cat("Group counts:\n")
print(table(merged_data$Group))

# Shared colour scale so both panels use identical colours
ancestry_levels <- sort(unique(merged_data$Ancestry_Prefix))
ancestry_colors <- setNames(scales::hue_pal()(length(ancestry_levels)), ancestry_levels)

pcs <- c("C2", "C3", "C4", "C5")


# ============================================================ #
#  PCA plots — two separate panels (1KG | Inferred)
# ============================================================ #

for (pc in pcs) {

  # Shared axis limits across both groups so panels are directly comparable
  x_range <- range(merged_data$C1,          na.rm = TRUE)
  y_range <- range(merged_data[[pc]], na.rm = TRUE)

  shared_coords <- coord_cartesian(xlim = x_range, ylim = y_range)

  # ---------- Left panel: 1KG reference ---------- #
  p_1kg <- ggplot(
    merged_data %>% filter(Group == "1KG"),
    aes(x = C1, y = .data[[pc]], color = Ancestry_Prefix)
  ) +
    geom_point(shape = 17, size = 1, alpha = 0.7) +
    scale_color_manual(values = ancestry_colors, drop = FALSE) +
    shared_coords +
    theme_minimal(base_size = 12) +
    labs(title = "1KG Reference", x = "PC1", y = pc, color = "Ancestry") +
    theme(
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  # ---------- Right panel: Inferred samples ---------- #
  p_inf <- ggplot(
    merged_data %>% filter(Group == "Inferred"),
    aes(x = C1, y = .data[[pc]], color = Ancestry_Prefix)
  ) +
    geom_point(shape = 16, size = 1, alpha = 0.7) +
    scale_color_manual(values = ancestry_colors, drop = FALSE) +
    shared_coords +
    theme_minimal(base_size = 12) +
    labs(title = "Inferred Ancestry", x = "PC1", y = pc, color = "Ancestry") +
    theme(
      legend.position  = "right",
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white")
    )

  p_panels <- p_1kg + p_inf +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(title = paste0("PC1 vs ", pc))

  ggsave(
    paste0(prefix, "_pca_plot_PC1_vs_", pc, "_panels.jpeg"),
    p_panels, width = 16, height = 6, dpi = 600, device = "jpeg"
  )
}
