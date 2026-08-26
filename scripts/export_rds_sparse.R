#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
input <- file.path(
  project_root, "data", "processed", "GSE205335_Lung_IO_UMI_matrix.rds"
)
output_dir <- file.path(project_root, "data", "processed", "sparse_export")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matrix <- readRDS(input)
stopifnot(inherits(matrix, "dgCMatrix"))
stopifnot(all(matrix@x >= 0), all(matrix@x == round(matrix@x)))

writeBin(as.integer(matrix@x), file.path(output_dir, "data.int32.bin"), size = 4)
writeBin(as.integer(matrix@i), file.path(output_dir, "indices.int32.bin"), size = 4)
writeBin(as.integer(matrix@p), file.path(output_dir, "indptr.int32.bin"), size = 4)
writeLines(rownames(matrix), file.path(output_dir, "genes.txt"), useBytes = TRUE)
writeLines(colnames(matrix), file.path(output_dir, "cells.txt"), useBytes = TRUE)

metadata <- data.frame(
  n_genes = nrow(matrix),
  n_cells = ncol(matrix),
  nnz = length(matrix@x)
)
write.csv(metadata, file.path(output_dir, "matrix_shape.csv"), row.names = FALSE)
cat("Exported", nrow(matrix), "genes x", ncol(matrix), "cells with",
    length(matrix@x), "nonzero values\n")

