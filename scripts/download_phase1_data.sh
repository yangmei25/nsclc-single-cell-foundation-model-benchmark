#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_dir="${project_root}/data/raw"
source_dir="${project_root}/data/metadata/source"

mkdir -p "${raw_dir}" "${source_dir}"

curl -L --fail --show-error \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  -o "${raw_dir}/GSE205335_Lung_IO_UMI_matrix.rds.gz"

curl -L --fail --show-error \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz" \
  -o "${source_dir}/GSE205335_Lung_IO_CellIdentity.txt.gz"

curl -L --fail --show-error \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/miniml/GSE205335_family.xml.tgz" \
  -o "${source_dir}/GSE205335_family.xml.tgz"

curl -L --fail --show-error \
  "https://cdn.elifesciences.org/articles/98366/elife-98366-supp1-v1.xlsx" \
  -o "${source_dir}/elife-98366-supp1-v1.xlsx"

# GEO wrapped the already-compressed RDS in an additional gzip layer.
# Keep the download unchanged and remove only that outer layer.
gzip -dc "${raw_dir}/GSE205335_Lung_IO_UMI_matrix.rds.gz" \
  > "${project_root}/data/processed/GSE205335_Lung_IO_UMI_matrix.rds"

echo "Phase I source downloads and matrix decompression complete."
