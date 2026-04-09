# 5_report.R — ancestry assignment summary report
#
# Usage: Rscript 5_report.R <outname> [submit_dir]
#   e.g. Rscript 5_report.R clz3a
#        Rscript 5_report.R clz3a /path/to/anc

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 5_report.R <outname> [submit_dir]")

outname    <- args[1]
submit_dir <- ifelse(length(args) >= 2, args[2], getwd())
setwd(submit_dir)

prefix   <- paste0("1kg_", outname, ".geno.05.pruned")
ann_file <- paste0(prefix, ".Q.annotated")

if (!file.exists(ann_file)) stop("File not found: ", ann_file)

df <- read.table(ann_file, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")

# Split 1KG reference vs target samples
ref    <- df[df$Ancestry_Pheno != "-", ]
target <- df[df$Ancestry_Pheno == "-", ]

# Superpop from Ancestry_Inf
target$Superpop <- sub("_.*", "", target$Ancestry_Inf)

# Pre-QC count from merged fam (1KG + cohort before --mind/--geno filters)
n_ref  <- nrow(ref)
merged_fam <- paste0("1kg_", outname, ".fam")
if (file.exists(merged_fam)) {
  n_preqc <- nrow(read.table(merged_fam, header = FALSE)) - n_ref
} else {
  n_preqc <- NA
}

n_postqc  <- nrow(target)
n_removed <- if (!is.na(n_preqc)) n_preqc - n_postqc else NA

tbl        <- sort(table(target$Superpop), decreasing = TRUE)
n_mix      <- sum(target$Superpop == "MIX")
n_assigned <- n_postqc - n_mix

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

# ---- Print report ---- #
report_lines <- c(
  "",
  "================================================",
  paste("  Ancestry Assignment Report:", outname),
  "================================================",
  if (!is.na(n_preqc)) sprintf("  Input (pre-QC)        : %d", n_preqc),
  if (!is.na(n_removed)) sprintf("  Removed (QC filters)  : %d", n_removed),
  sprintf("  Retained (post-QC)    : %d", n_postqc),
  sprintf("  Assigned (non-MIX)    : %d  (%.1f%%)", n_assigned, 100 * n_assigned / n_postqc),
  sprintf("  Unassigned (MIX)      : %d  (%.1f%%)", n_mix,      100 * n_mix      / n_postqc),
  "",
  "  Breakdown by superpopulation:",
  sapply(names(tbl), function(pop) {
    sprintf("    %-6s : %6d  (%5.1f%%)", pop, tbl[[pop]], 100 * tbl[[pop]] / n_postqc)
  }),
  "------------------------------------------------",
  sprintf("  1KG reference samples : %d", n_ref),
  "================================================",
  ""
)
cat(paste(report_lines, collapse = "\n"), "\n")

# Write to file
out_file <- paste0(outname, "_ancestry_report.txt")
writeLines(report_lines, out_file)
cat("Report written to:", out_file, "\n")

# ---- Single-line Excel-friendly output ---- #
superpops   <- c("AFR", "AMR", "EAS", "EUR", "SAS", "MIX")
header_vals <- c("abbreviation", "n_preqc", "n_removed", "n_postqc",
                 "n_assigned", "pct_assigned",
                 paste0("n_", superpops), paste0("pct_", superpops))
row_vals <- c(
  outname,
  if (!is.na(n_preqc))   n_preqc   else "NA",
  if (!is.na(n_removed)) n_removed else "NA",
  n_postqc,
  n_assigned,
  round(100 * n_assigned / n_postqc, 1),
  sapply(superpops, function(p) { n <- tbl[p]; if (is.na(n)) 0L else as.integer(n) }),
  sapply(superpops, function(p) { n <- tbl[p]; if (is.na(n)) 0 else round(100 * n / n_postqc, 1) })
)

cat("\n--- Excel row (space-separated) ---\n")
cat(paste(header_vals, collapse = " "), "\n")
cat(paste(row_vals,    collapse = " "), "\n\n")
