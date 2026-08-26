#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readxl))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
input <- file.path(
  project_root, "data", "metadata", "source", "elife-98366-supp1-v1.xlsx"
)
output <- file.path(project_root, "data", "metadata", "clinical_metadata.csv")

clinical <- read_excel(input, sheet = "Clinical_info", skip = 2)
clinical <- clinical[!is.na(clinical$Sample), ]
clinical <- clinical[!grepl("^\\*", clinical$Sample), ]
write.csv(clinical, output, row.names = FALSE, na = "")
cat("Wrote", nrow(clinical), "clinical sample records to", output, "\n")
