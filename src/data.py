"""Phase I data loading and AnnData construction utilities."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def normalized_sample_key(sample: str, platform: str | None = None) -> str:
    """Create a join key for paper sample IDs and annotation `orig.ident` IDs."""
    value = str(sample)
    if platform is not None and str(platform).strip():
        value += str(platform)
    return re.sub(r"[^a-z0-9]", "", value.lower())


def load_cell_metadata(project_root: Path) -> pd.DataFrame:
    path = (
        project_root
        / "data/metadata/source/GSE205335_Lung_IO_CellIdentity.txt.gz"
    )
    cells = pd.read_csv(path, sep="\t", keep_default_na=False)
    if not cells["barcode"].is_unique:
        raise ValueError("Cell barcodes are not unique")
    cells["sample_key"] = cells["orig.ident"].map(normalized_sample_key)
    return cells.set_index("barcode", drop=False)


def load_clinical_metadata(project_root: Path) -> pd.DataFrame:
    path = project_root / "data/metadata/clinical_metadata.csv"
    clinical = pd.read_csv(path, keep_default_na=False)
    clinical["sample_key"] = [
        normalized_sample_key(sample, platform)
        for sample, platform in zip(
            clinical["Sample"], clinical["Platform"], strict=True
        )
    ]
    if not clinical["sample_key"].is_unique:
        raise ValueError("Clinical sample join keys are not unique")
    sampling_date = pd.to_datetime(clinical["Sampling date"], errors="coerce")
    treatment_date = pd.to_datetime(
        clinical["Immunotherapy start date"], errors="coerce"
    )
    clinical["timepoint_derived"] = "unknown_or_same_day"
    clinical.loc[sampling_date < treatment_date, "timepoint_derived"] = "pre"
    clinical.loc[sampling_date > treatment_date, "timepoint_derived"] = "post"
    return clinical


def build_anndata(project_root: Path):
    import anndata as ad

    export = project_root / "data/processed/sparse_export"
    shape = pd.read_csv(export / "matrix_shape.csv").iloc[0]
    n_genes = int(shape["n_genes"])
    n_cells = int(shape["n_cells"])
    nnz = int(shape["nnz"])

    data = np.memmap(export / "data.int32.bin", dtype=np.int32, mode="r", shape=(nnz,))
    indices = np.memmap(
        export / "indices.int32.bin", dtype=np.int32, mode="r", shape=(nnz,)
    )
    indptr = np.memmap(
        export / "indptr.int32.bin", dtype=np.int32, mode="r", shape=(n_cells + 1,)
    )
    genes = pd.Index((export / "genes.txt").read_text().splitlines(), name="gene")
    cell_ids = pd.Index((export / "cells.txt").read_text().splitlines(), name="cell_id")
    if len(genes) != n_genes or len(cell_ids) != n_cells:
        raise ValueError("Exported names do not match matrix dimensions")

    gene_by_cell = sparse.csc_matrix(
        (data, indices, indptr), shape=(n_genes, n_cells), copy=False
    )
    cell_by_gene = gene_by_cell.T.tocsr()

    cells = load_cell_metadata(project_root)
    missing_annotation = cell_ids.difference(cells.index)
    extra_annotation = cells.index.difference(cell_ids)
    if len(missing_annotation) or len(extra_annotation):
        raise ValueError(
            f"Expression/annotation mismatch: {len(missing_annotation)} expression-only, "
            f"{len(extra_annotation)} annotation-only cells"
        )
    cells = cells.loc[cell_ids].copy()

    clinical = load_clinical_metadata(project_root).set_index("sample_key")
    missing_samples = sorted(set(cells["sample_key"]) - set(clinical.index))
    if missing_samples:
        raise ValueError(f"Cell samples absent from clinical metadata: {missing_samples}")
    cells = cells.join(clinical, on="sample_key", rsuffix="_clinical")

    var = pd.DataFrame(index=genes)
    var["mt"] = var.index.str.startswith("MT-")
    adata = ad.AnnData(X=cell_by_gene, obs=cells, var=var)
    adata.uns["source_geo"] = "GSE205335"
    adata.uns["source_publication_doi"] = "10.7554/eLife.98366"
    adata.uns["counts_layer"] = "X contains raw UMI counts"
    return adata
