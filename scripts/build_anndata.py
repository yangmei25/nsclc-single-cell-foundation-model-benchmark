#!/usr/bin/env python3
"""Construct the initial GSE205335 AnnData object from the sparse R export."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.data import build_anndata  # noqa: E402


def main() -> int:
    adata = build_anndata(ROOT)
    output = ROOT / "data/processed/GSE205335_phase1_raw_counts.h5ad"
    adata.write_h5ad(output, compression="gzip")
    print(f"Wrote {adata.n_obs:,} cells x {adata.n_vars:,} genes to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

