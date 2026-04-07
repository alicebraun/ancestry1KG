#!/usr/bin/env Rscript
# 0b_harmonize_alleles.R
#
# Harmonize alleles in an input PLINK BIM against a reference BIM, matching
# variants by chromosome and position.
#
# Outputs:
#   --out-flip    : variant IDs where a complement strand flip resolves the match
#   --out-exclude : variant IDs to remove (not in reference, allele mismatch,
#                   or palindromic A/T / C/G SNP that cannot be strand-resolved)
#
# Usage:
#   Rscript 0b_harmonize_alleles.R \
#     --input-bim input.bim \
#     --ref-bim   reference.bim \
#     --out-flip  to_flip.txt \
#     --out-exclude to_exclude.txt

suppressPackageStartupMessages(library(data.table))

# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    if      (args[i] == "--input-bim")   { out$input_bim   <- args[i+1]; i <- i+2 }
    else if (args[i] == "--ref-bim")     { out$ref_bim     <- args[i+1]; i <- i+2 }
    else if (args[i] == "--out-flip")    { out$out_flip    <- args[i+1]; i <- i+2 }
    else if (args[i] == "--out-exclude") { out$out_exclude <- args[i+1]; i <- i+2 }
    else stop(paste("Unknown argument:", args[i]))
  }
  missing <- c("input_bim","ref_bim","out_flip","out_exclude")[
    !c("input_bim","ref_bim","out_flip","out_exclude") %in% names(out)]
  if (length(missing)) stop(paste("Missing:", paste(missing, collapse=", ")))
  out
}

opt <- parse_args(args)

# --------------------------------------------------------------------------- #
# Read BIM files
# --------------------------------------------------------------------------- #
bim_cols <- c("chr", "id", "cm", "pos", "a1", "a2")

cat("Reading input BIM:", opt$input_bim, "\n")
inp <- fread(opt$input_bim, header = FALSE, col.names = bim_cols)

cat("Reading reference BIM:", opt$ref_bim, "\n")
ref <- fread(opt$ref_bim, header = FALSE, col.names = bim_cols)

# Normalise
inp[, chr := sub("^chr", "", chr)]
ref[, chr := sub("^chr", "", chr)]
inp[, c("a1","a2") := .(toupper(a1), toupper(a2))]
ref[, c("a1","a2") := .(toupper(a1), toupper(a2))]

# --------------------------------------------------------------------------- #
# Vectorised allele harmonisation via data.table join
# --------------------------------------------------------------------------- #

# Keep first occurrence per position in reference
ref <- ref[!duplicated(paste(chr, pos, sep = ":"))]
setnames(ref, c("a1","a2"), c("ra1","ra2"))

# Join on chr + pos
inp <- ref[, .(chr, pos, ra1, ra2)][inp, on = .(chr, pos)]

# Complement lookup (vectorised)
comp <- function(x) chartr("ACGT", "TGCA", x)

inp[, `:=`(
  ca1 = comp(a1),
  ca2 = comp(a2)
)]

# Classify each variant
inp[, status := fcase(
  is.na(ra1),                                          "exclude",   # not in reference
  (a1 == ra1 & a2 == ra2) | (a1 == ra2 & a2 == ra1),  "ok",        # direct match
  comp(a1) == a2,                                      "exclude",   # palindromic
  (ca1 == ra1 & ca2 == ra2) | (ca1 == ra2 & ca2 == ra1), "flip",   # complement match
  default =                                            "exclude"    # no match
)]

# --------------------------------------------------------------------------- #
# Write outputs
# --------------------------------------------------------------------------- #
writeLines(inp[status == "flip",    id], con = opt$out_flip)
writeLines(inp[status == "exclude", id], con = opt$out_exclude)

cat(sprintf("Variants processed : %d\n", nrow(inp)))
cat(sprintf("  OK (direct match): %d\n", inp[status == "ok",      .N]))
cat(sprintf("  To flip          : %d\n", inp[status == "flip",    .N]))
cat(sprintf("  To exclude       : %d\n", inp[status == "exclude", .N]))
