#!/usr/bin/env python3
"""Generate frozen scGPT embeddings in an isolated Python environment."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ["MPLBACKEND"] = "Agg"

import numpy as np
import scanpy as sc
import torch
import scgpt as scg
from scgpt.tasks import cell_emb

LABEL_COLUMNS = ["cluster.total", "lineage.total", "lineage.sub", "cluster.sub", "celltype"]
OBS_TO_SAVE = ["Sample", "Patient", "Platform", "Tissue origin"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--device", choices=["cuda", "cpu"], default="cuda")
    return parser.parse_args()


def force_single_process_dataloader() -> None:
    """Prevent scGPT 0.2.4 from duplicating its dense matrix across workers."""
    original = cell_emb.DataLoader

    def safe_dataloader(*args, **kwargs):
        kwargs["num_workers"] = 0
        kwargs["pin_memory"] = False
        return original(*args, **kwargs)

    cell_emb.DataLoader = safe_dataloader


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    source_path = root / "data/processed/GSE205335_phase1_raw_counts.h5ad"
    model_dir = root / "models/scGPT_human"
    output_path = root / "results/embeddings/GSE205335_scgpt_frozen_v1.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable in the isolated scGPT process")
    np.random.seed(args.random_state)
    torch.manual_seed(args.random_state)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.random_state)
    torch.set_grad_enabled(False)

    adata = sc.read_h5ad(source_path)
    if adata.shape != (96_505, 33_714):
        raise ValueError(f"Unexpected source shape: {adata.shape}")
    cohort_mask = (
        adata.obs["core.patient"].astype(str).eq("Core")
        & adata.obs["Tissue origin"].astype(str).eq("Metastatic LN")
    )
    work = adata[cohort_mask].copy()
    if work.n_obs != 29_614 or work.obs["Sample"].nunique() != 10:
        raise ValueError(
            f"Expected 29,614 cells from 10 site-matched samples; got "
            f"{work.n_obs:,} cells from {work.obs['Sample'].nunique()} samples"
        )
    work.obs.drop(columns=LABEL_COLUMNS, inplace=True)
    work.var["gene_name"] = work.var_names.astype(str)
    del adata

    device_name = torch.cuda.get_device_name(0) if args.device == "cuda" else "CPU"
    print(f"Embedding {work.n_obs:,} cells on {device_name}")
    force_single_process_dataloader()
    embedding = scg.tasks.embed_data(
        work,
        model_dir,
        gene_col="gene_name",
        obs_to_save=OBS_TO_SAVE,
        batch_size=args.batch_size,
        device=args.device,
        use_fast_transformer=False,
        return_new_adata=True,
    )
    if embedding.n_obs != work.n_obs or not embedding.obs_names.equals(work.obs_names):
        raise ValueError("Embedding cells do not match the selected input")
    if not np.isfinite(embedding.X).all():
        raise ValueError("Embedding contains non-finite values")
    embedding.write_h5ad(output_path, compression="gzip")
    print(f"Wrote {output_path} with shape {embedding.shape}")


if __name__ == "__main__":
    main()
