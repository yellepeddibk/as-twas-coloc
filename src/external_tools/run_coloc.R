#!/usr/bin/env Rscript
# run_coloc.R
# -----------
# Companion R script for the Python coloc_r_wrapper.py interface.
#
# Usage (called by Python subprocess):
#   Rscript run_coloc.R <input_json> <output_json>
#
# Input JSON schema:
#   {
#     "dataset1": { "snp": [...], "beta": [...], "varbeta": [...], "N": int, "type": "cc"|"quant", "s": float },
#     "dataset2": { "snp": [...], "beta": [...], "varbeta": [...], "N": int, "type": "quant" },
#     "priors":   { "p1": 1e-4, "p2": 1e-4, "p12": 1e-5 }
#   }
#
# Output JSON schema:
#   {
#     "PP.H0.abf": ..., "PP.H1.abf": ..., "PP.H2.abf": ...,
#     "PP.H3.abf": ..., "PP.H4.abf": ..., "nsnps": int
#   }
#
# Pitfall notes:
#   - dataset1 and dataset2 MUST have "varbeta", NOT "se".
#     varbeta = se^2.  Passing se instead of varbeta produces incorrect posteriors.
#   - Do NOT pre-filter SNPs by p-value before calling coloc.abf.
#     Removing non-significant SNPs inflates PP4 artificially.
#   - coloc.abf assumes each dataset has a single causal variant per locus.
#     For multi-causal loci, use coloc.susie (requires SuSiE fine-mapping first).
#
# Requirements:
#   install.packages("coloc")    # CRAN
#   install.packages("jsonlite") # CRAN

suppressPackageStartupMessages({
  library(coloc)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript run_coloc.R <input_json> <output_json>")
}

input_file  <- args[1]
output_file <- args[2]

# ── Load input ────────────────────────────────────────────────────────────────
payload <- fromJSON(input_file, simplifyVector = FALSE)

ds1     <- payload$dataset1
ds2     <- payload$dataset2
priors  <- payload$priors

# Convert list elements to named vectors as required by coloc
to_named_vec <- function(x, name = "snp") {
  vals <- unlist(x)
  vals
}

dataset1 <- list(
  snp      = unlist(ds1$snp),
  beta     = as.numeric(unlist(ds1$beta)),
  varbeta  = as.numeric(unlist(ds1$varbeta)),
  N        = as.integer(ds1$N),
  type     = ds1$type
)
if (!is.null(ds1$s)) dataset1$s <- as.numeric(ds1$s)

dataset2 <- list(
  snp      = unlist(ds2$snp),
  beta     = as.numeric(unlist(ds2$beta)),
  varbeta  = as.numeric(unlist(ds2$varbeta)),
  N        = as.integer(ds2$N),
  type     = ds2$type
)
if (!is.null(ds2$sdY)) dataset2$sdY <- as.numeric(ds2$sdY)

# ── Run coloc.abf ─────────────────────────────────────────────────────────────
result <- tryCatch(
  coloc.abf(
    dataset1 = dataset1,
    dataset2 = dataset2,
    p1  = as.numeric(priors$p1),
    p2  = as.numeric(priors$p2),
    p12 = as.numeric(priors$p12)
  ),
  error = function(e) {
    message("coloc.abf error: ", conditionMessage(e))
    stop(e)
  }
)

# ── Extract posterior probabilities ───────────────────────────────────────────
summary_row <- as.list(result$summary)

output <- list(
  "PP.H0.abf" = as.numeric(summary_row[["PP.H0.abf"]]),
  "PP.H1.abf" = as.numeric(summary_row[["PP.H1.abf"]]),
  "PP.H2.abf" = as.numeric(summary_row[["PP.H2.abf"]]),
  "PP.H3.abf" = as.numeric(summary_row[["PP.H3.abf"]]),
  "PP.H4.abf" = as.numeric(summary_row[["PP.H4.abf"]]),
  "nsnps"     = as.integer(summary_row[["nsnps"]])
)

# ── Write output JSON ─────────────────────────────────────────────────────────
writeLines(toJSON(output, auto_unbox = TRUE, digits = 8), output_file)
message("coloc.abf complete. PP4 = ", round(output[["PP.H4.abf"]], 4))
