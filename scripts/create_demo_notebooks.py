#!/usr/bin/env python3
"""Generate the five-notebook scRNA/scGPT demo workflow."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKER_SOURCE = (ROOT / "scripts/run_scgpt_embedding.py").read_text()


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}


def py(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": text.splitlines(True)}


def write(name, cells, display="Python (scrna-oncology-fm)"):
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": display, "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = ROOT / "notebooks" / name
    path.write_text(json.dumps(notebook, indent=1) + "\n")
    print("Generated", path)


setup = '''from pathlib import Path
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path.cwd().resolve()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures/demo"
TABLES = RESULTS / "tables"
EMBEDDINGS = RESULTS / "embeddings"
for directory in (FIGURES, TABLES, EMBEDDINGS):
    directory.mkdir(parents=True, exist_ok=True)
RAW_PATH = DATA / "processed/GSE205335_phase1_raw_counts.h5ad"
RANDOM_STATE = 0
'''


# ---------------------------------------------------------------------------
# Notebook 1
# ---------------------------------------------------------------------------
write("01_data_download_and_preparation.ipynb", [
    md('''# 1 — Download, prepare, and audit the complete dataset

## Question

Can the complete GSE205335 dataset be reconstructed as one reusable AnnData
object with raw counts, sample metadata, clinical metadata, and published cell
annotations aligned correctly?

This notebook keeps **all 96,505 cells and 33,714 genes** so the project can be
extended later. It does not perform response analysis or model evaluation.'''),
    md('''## 0. Reproducible source preparation

For a fresh checkout, run these commands once from a terminal:

```bash
bash scripts/download_phase1_data.sh
Rscript scripts/export_rds_sparse.R
Rscript scripts/prepare_clinical_metadata.R
python scripts/build_anndata.py
```

The cells below audit the resulting file rather than redownloading 3 GB every
time the notebook is opened.'''),
    py(setup + '''
assert RAW_PATH.exists(), (
    f"Missing {RAW_PATH}. Run the four preparation commands shown above."
)
adata = ad.read_h5ad(RAW_PATH, backed="r")
print(adata)'''),
    md('''## 1. Structural and count-matrix checks'''),
    py('''assert adata.shape == (96_505, 33_714)
assert adata.obs_names.is_unique
assert adata.var_names.is_unique
for start in range(0, adata.n_obs, 10_000):
    assert adata.X[start:min(start + 10_000, adata.n_obs)].min() >= 0
required = ["Sample", "Patient", "Tissue origin", "Platform", "core.patient",
            "Response", "lineage.total", "lineage.sub", "celltype"]
assert not set(required) - set(adata.obs.columns)
print("Shape:", adata.shape)
print("Samples:", adata.obs["Sample"].nunique())
print("Patients:", adata.obs["Patient"].nunique())
print("Sparse/raw non-negative counts verified across all cells")'''),
    md('''## 2. Dataset inventory and fixed demo cohort'''),
    py('''inventory = pd.DataFrame({
    "cells": adata.obs.groupby("Sample", observed=True).size(),
    "Patient": adata.obs.groupby("Sample", observed=True)["Patient"].first(),
    "Tissue": adata.obs.groupby("Sample", observed=True)["Tissue origin"].first(),
    "Platform": adata.obs.groupby("Sample", observed=True)["Platform"].first(),
    "Core": adata.obs.groupby("Sample", observed=True)["core.patient"].first(),
})
inventory.to_csv(TABLES / "dataset_sample_inventory.csv")
display(inventory)

demo_mask = adata.obs["core.patient"].astype(str).eq("Core") & adata.obs["Tissue origin"].astype(str).eq("Metastatic LN")
demo = adata.obs.loc[demo_mask].groupby("Sample", observed=True).agg(
    cells=("Sample", "size"), Patient=("Patient", "first"),
    Response=("Response", "first"), Platform=("Platform", "first"),
).sort_index()
assert demo["cells"].sum() == 29_614 and len(demo) == 10
demo.to_csv(TABLES / "demo_10sample_cohort.csv")
display(demo)
print("Demo cohort: 29,614 cells, 10 metastatic-LN samples")
adata.file.close()'''),
    md('''## Output

- Complete reusable dataset: `data/processed/GSE205335_phase1_raw_counts.h5ad`
- Complete sample inventory: `results/tables/dataset_sample_inventory.csv`
- Fixed demo cohort: `results/tables/demo_10sample_cohort.csv`'''),
])


# ---------------------------------------------------------------------------
# Notebook 2
# ---------------------------------------------------------------------------
write("02_scanpy_pca_10_samples.ipynb", [
    md('''# 2 — Scanpy PCA representations for the 10 samples

## Questions

1. Does conventional Scanpy PCA separate broad cell types in all 29,614 cells?
2. After subsetting, does CD4-only PCA resolve CD4 states?
3. Does CD8-only PCA resolve CD8 states?

Published annotations are used only to color and evaluate the unsupervised PCA
spaces; they are not inputs to normalization, HVG selection, PCA, or Leiden.'''),
    py(setup + '''
source = sc.read_h5ad(RAW_PATH)
mask = source.obs["core.patient"].astype(str).eq("Core") & source.obs["Tissue origin"].astype(str).eq("Metastatic LN")
demo = source[mask].copy()
del source
assert demo.n_obs == 29_614 and demo.obs["Sample"].nunique() == 10
print(demo)'''),
    md('''## 1. Shared Scanpy pipeline'''),
    py('''def build_pca(raw, label, resolution=0.6):
    work = raw.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=2_000, flavor="seurat_v3", batch_key="Platform")
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.pca(work, n_comps=50, mask_var="highly_variable", random_state=RANDOM_STATE)
    sc.pp.neighbors(work, use_rep="X_pca", n_neighbors=15, n_pcs=50, random_state=RANDOM_STATE)
    sc.tl.umap(work, random_state=RANDOM_STATE)
    sc.tl.leiden(work, resolution=resolution, key_added=f"leiden_{label}",
                 random_state=RANDOM_STATE, flavor="igraph", n_iterations=2, directed=False)
    return work

all_cells = build_pca(demo, "all")
cd4 = build_pca(demo[demo.obs["lineage.sub"].astype(str).eq("CD4+ T cells")], "cd4", 0.8)
cd8 = build_pca(demo[demo.obs["lineage.sub"].astype(str).eq("CD8+ T cells")], "cd8", 0.8)
print("All/CD4/CD8:", all_cells.n_obs, cd4.n_obs, cd8.n_obs)'''),
    md('''## 2. Three PCA/UMAP views'''),
    py('''plots = [
    (all_cells, "lineage.total", "All cells — broad lineages", "scanpy_all_cells_pca_umap.png"),
    (cd4, "celltype", "CD4-only — fine states", "scanpy_cd4_pca_umap.png"),
    (cd8, "celltype", "CD8-only — fine states", "scanpy_cd8_pca_umap.png"),
]
for obj, color, title, filename in plots:
    with plt.rc_context({"figure.figsize": (9, 6)}):
        sc.pl.umap(obj, color=color, title=title, show=False)
        plt.savefig(FIGURES / filename, dpi=220, bbox_inches="tight")
        plt.show()'''),
    md('''## 3. Save compact PCA representations'''),
    py('''def save_representation(obj, path, leiden_key):
    out = ad.AnnData(X=obj.obsm["X_pca"].astype(np.float32), obs=obj.obs.copy())
    out.obsm["X_umap"] = obj.obsm["X_umap"].astype(np.float32)
    out.uns["representation"] = "Scanpy log-normalization, 2,000 batch-aware HVGs, 50-PC PCA"
    out.uns["leiden_key"] = leiden_key
    out.write_h5ad(path, compression="gzip")

save_representation(all_cells, EMBEDDINGS / "GSE205335_10sample_scanpy_pca_v1.h5ad", "leiden_all")
save_representation(cd4, EMBEDDINGS / "GSE205335_10sample_cd4_scanpy_pca_v1.h5ad", "leiden_cd4")
save_representation(cd8, EMBEDDINGS / "GSE205335_10sample_cd8_scanpy_pca_v1.h5ad", "leiden_cd8")
print("Saved three PCA representations")'''),
])


# ---------------------------------------------------------------------------
# Notebook 3
# ---------------------------------------------------------------------------
write("03_scgpt_cell_embeddings_10_samples.ipynb", [
    md('''# 3 — Frozen scGPT embeddings for all cells in the 10 samples

## Question

Can frozen scGPT generate a reusable 512-dimensional embedding for every one
of the same 29,614 cells used by the Scanpy baseline?

Run this notebook on Colab GPU. No cell-type or response label enters scGPT.'''),
    md('''## 0. Install one lightweight file-reading dependency

The default Colab kernel needs only `anndata` to verify input/output files.
Scanpy and scGPT are installed later inside the isolated Python 3.11
environment. Let this cell finish; it may take about a minute on a new runtime.'''),
    py('''import importlib.util, subprocess, sys
if importlib.util.find_spec("anndata") is None:
    print("Installing anndata in the default Colab kernel; please wait...")
    subprocess.run([sys.executable, "-m", "pip", "install", "anndata"], check=True)
print("anndata is available")'''),
    py('''from pathlib import Path
import subprocess, sys
import numpy as np
import anndata as ad
import torch

IN_COLAB = Path("/content").exists()
if IN_COLAB:
    from google.colab import drive
    drive.mount("/content/drive")
    ROOT = Path("/content/drive/MyDrive/scrna_oncology_fm")
else:
    ROOT = Path.cwd().resolve()
    if ROOT.name == "notebooks": ROOT = ROOT.parent

MODEL_DIR = ROOT / "models/scGPT_human"
OUTPUT = ROOT / "results/embeddings/GSE205335_scgpt_frozen_v1.h5ad"
SCGPT_PYTHON = Path("/content/scgpt-venv/bin/python") if IN_COLAB else Path(sys.executable)
for name in ["best_model.pt", "args.json", "vocab.json"]:
    assert (MODEL_DIR / name).exists(), f"Missing {MODEL_DIR / name}"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
if IN_COLAB and DEVICE != "cuda": raise RuntimeError("Select a Colab GPU runtime")
print("Device:", DEVICE, torch.cuda.get_device_name(0) if DEVICE == "cuda" else "CPU")'''),
    md('''## 1. Create/reuse the isolated scGPT environment'''),
    py('''needs_install = IN_COLAB and (not SCGPT_PYTHON.exists() or subprocess.run(
    [str(SCGPT_PYTHON), "-c", "import torch, torchtext, scgpt"], capture_output=True).returncode != 0)
if needs_install:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "uv"], check=True)
    if not SCGPT_PYTHON.exists():
        subprocess.run(["uv", "python", "install", "3.11"], check=True)
        subprocess.run(["uv", "venv", str(SCGPT_PYTHON.parent.parent), "--python", "3.11"], check=True)
    subprocess.run(["uv", "pip", "install", "--python", str(SCGPT_PYTHON), "torch==2.3.0",
                    "--index-url", "https://download.pytorch.org/whl/cu121"], check=True)
    subprocess.run(["uv", "pip", "install", "--python", str(SCGPT_PYTHON), "numpy<2",
                    "torchtext==0.18.0", "scgpt==0.2.4", "fsspec[http]<=2024.6.1",
                    "scanpy", "anndata", "pandas", "scipy", "ipython"], check=True)

# scGPT imports IPython even in the non-interactive worker. Repair an existing
# environment created by an earlier notebook version without reinstalling it.
if subprocess.run([str(SCGPT_PYTHON), "-c", "import IPython"], capture_output=True).returncode != 0:
    subprocess.run(["uv", "pip", "install", "--python", str(SCGPT_PYTHON), "ipython"], check=True)
probe = subprocess.run([
    str(SCGPT_PYTHON), "-c",
    "import os; os.environ['MPLBACKEND']='Agg'; "
    "import torch; print('torch', torch.__version__); "
    "import torchtext; print('torchtext', torchtext.__version__); "
    "import scgpt; print('scGPT', scgpt.__version__); "
    "print('CUDA', torch.cuda.is_available())",
], capture_output=True, text=True)
print(probe.stdout)
if probe.returncode != 0:
    print("scGPT environment probe failed. Full error:")
    print(probe.stderr)
    raise RuntimeError("scGPT environment is incomplete; use the full error printed above.")
print("scGPT environment is ready")'''),
    md('''## 2. Verify cohort and embed all cells'''),
    py(f'''source = ad.read_h5ad(ROOT / "data/processed/GSE205335_phase1_raw_counts.h5ad", backed="r")
mask = source.obs["core.patient"].astype(str).eq("Core") & source.obs["Tissue origin"].astype(str).eq("Metastatic LN")
assert int(mask.sum()) == 29_614 and source.obs.loc[mask, "Sample"].nunique() == 10
source.file.close()

worker = Path("/content/run_scgpt_embedding.py") if IN_COLAB else ROOT / "scripts/run_scgpt_embedding.py"
worker.write_text({WORKER_SOURCE!r})
RUN_EMBEDDING = True
if RUN_EMBEDDING:
    subprocess.run([str(SCGPT_PYTHON), str(worker), "--root", str(ROOT), "--batch-size", "16",
                    "--random-state", "0", "--device", DEVICE], check=True)

embedding = ad.read_h5ad(OUTPUT)
assert embedding.shape == (29_614, 512)
assert embedding.obs["Sample"].nunique() == 10
assert np.isfinite(embedding.X).all()
print("Complete scGPT embedding:", embedding.shape)
print("Saved:", OUTPUT)'''),
])


# ---------------------------------------------------------------------------
# Notebook 4
# ---------------------------------------------------------------------------
write("04_cell_level_pca_vs_scgpt.ipynb", [
    md('''# 4 — Cell-level benchmark: Scanpy PCA versus frozen scGPT

## Question

On exactly the same cells, which representation better preserves known broad
cell identities and fine CD4/CD8 states?

Evidence includes UMAPs, label silhouette, kNN label accuracy, Leiden ARI/NMI,
and platform silhouette. UMAP appearance alone is not treated as proof.'''),
    py(setup + '''
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.neighbors import KNeighborsClassifier

PCA_PATH = EMBEDDINGS / "GSE205335_10sample_scanpy_pca_v1.h5ad"
SCGPT_PATH = EMBEDDINGS / "GSE205335_scgpt_frozen_v1.h5ad"
assert PCA_PATH.exists() and SCGPT_PATH.exists(), "Run Notebooks 2 and 3 first"
pca = sc.read_h5ad(PCA_PATH)
scgpt = sc.read_h5ad(SCGPT_PATH)
assert pca.obs_names.equals(scgpt.obs_names)
label_columns = ["lineage.total", "lineage.sub", "celltype"]
scgpt.obs[label_columns] = pca.obs.loc[scgpt.obs_names, label_columns].copy()
print("Matched cells:", pca.n_obs)'''),
    md('''## 1. Build comparable graphs and UMAPs'''),
    py('''views = {"Scanpy_PCA": pca.copy(), "frozen_scGPT": scgpt.copy()}
for name, obj in views.items():
    sc.pp.neighbors(obj, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
    sc.tl.umap(obj, random_state=RANDOM_STATE)
    sc.tl.leiden(obj, resolution=0.6, key_added="cluster", random_state=RANDOM_STATE,
                 flavor="igraph", n_iterations=2, directed=False)

fig, axes = plt.subplots(1, 2, figsize=(15, 6))
for ax, (name, obj) in zip(axes, views.items()):
    sc.pl.umap(obj, color="lineage.total", title=name, ax=ax, show=False)
fig.tight_layout()
fig.savefig(FIGURES / "cell_level_broad_pca_vs_scgpt.png", dpi=220, bbox_inches="tight")
plt.show()'''),
    md('''## 2. Quantitative broad-cell benchmark

Silhouette uses a fixed 5,000-cell subsample to avoid constructing a full
29,614 × 29,614 distance matrix. kNN is evaluated by sample-held-out folds, so
cells from the test sample never appear in training.'''),
    py('''def held_out_sample_knn(matrix, labels, samples, k=15):
    predictions = np.empty(len(labels), dtype=object)
    for sample in pd.unique(samples):
        test = samples == sample
        model = KNeighborsClassifier(n_neighbors=k, weights="distance", metric="cosine", n_jobs=-1)
        model.fit(matrix[~test], labels[~test])
        predictions[test] = model.predict(matrix[test])
    return float(np.mean(predictions == labels))

def benchmark(name, obj, labels, scope):
    matrix = np.asarray(obj.X, dtype=np.float32)
    labels = np.asarray(labels).astype(str)
    samples = obj.obs["Sample"].astype(str).to_numpy()
    platform = obj.obs["Platform"].astype(str).to_numpy()
    return {
        "scope": scope, "representation": name, "n_cells": len(labels),
        "label_silhouette": silhouette_score(matrix, labels, metric="cosine", sample_size=min(5000, len(labels)), random_state=0),
        "sample_held_out_knn_accuracy": held_out_sample_knn(matrix, labels, samples),
        "leiden_ARI": adjusted_rand_score(labels, obj.obs["cluster"].astype(str)),
        "leiden_NMI": normalized_mutual_info_score(labels, obj.obs["cluster"].astype(str)),
        "platform_silhouette": silhouette_score(matrix, platform, metric="cosine", sample_size=min(5000, len(labels)), random_state=0),
    }

rows = [benchmark(name, obj, obj.obs["lineage.total"], "broad_all_cells") for name, obj in views.items()]'''),
    md('''## 3. Fine-state CD4 and CD8 benchmark'''),
    py('''fig, axes = plt.subplots(2, 2, figsize=(16, 12))
for row_index, lineage in enumerate(["CD4+ T cells", "CD8+ T cells"]):
    for col_index, (name, obj) in enumerate(views.items()):
        mask = obj.obs["lineage.sub"].astype(str).eq(lineage)
        subset = obj[mask].copy()
        sc.pp.neighbors(subset, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
        sc.tl.umap(subset, random_state=RANDOM_STATE)
        sc.tl.leiden(subset, resolution=0.8, key_added="cluster", random_state=RANDOM_STATE,
                     flavor="igraph", n_iterations=2, directed=False)
        sc.pl.umap(subset, color="celltype", title=f"{lineage}: {name}", ax=axes[row_index, col_index], show=False)
        rows.append(benchmark(name, subset, subset.obs["celltype"], lineage))
fig.tight_layout()
fig.savefig(FIGURES / "cell_level_CD4_CD8_pca_vs_scgpt.png", dpi=220, bbox_inches="tight")
plt.show()

metrics = pd.DataFrame(rows)
display(metrics)
metrics.to_csv(TABLES / "cell_level_pca_vs_scgpt_metrics.csv", index=False)'''),
    md('''## Reading the evidence

- Higher label silhouette, held-out kNN accuracy, ARI and NMI are better.
- Platform silhouette should be near zero; a large positive value suggests a
  technical platform effect.
- A credible conclusion reports where scGPT wins, ties, or loses rather than
  assuming the foundation model must be superior.'''),
])


# ---------------------------------------------------------------------------
# Notebook 5
# ---------------------------------------------------------------------------
write("05_sample_level_embeddings.ipynb", [
    md('''# 5 — Sample-level response representations

## Questions

1. What does conventional cell-state composition show?
2. What do conventional CD4/CD8 functional gene programs show?
3. Do pooled scGPT embeddings add a complementary response pattern?

All representations are frozen before response is revealed. Response is shown
by color; the two samples from P1076 use triangles while all others use circles.'''),
    py(setup + '''
from itertools import combinations
from sklearn.metrics import pairwise_distances, silhouette_score

SCGPT_PATH = EMBEDDINGS / "GSE205335_scgpt_frozen_v1.h5ad"
assert SCGPT_PATH.exists(), "Run Notebook 3 first"
embedding = sc.read_h5ad(SCGPT_PATH)
source = sc.read_h5ad(RAW_PATH, backed="r")
meta = source.obs.loc[embedding.obs_names, ["Sample", "Patient", "Tissue origin", "Platform", "lineage.sub", "celltype"]].copy()
response_locked = source.obs.loc[embedding.obs_names, ["Sample", "Response"]].copy()
source.file.close()
assert embedding.n_obs == 29_614 and meta["Sample"].nunique() == 10
sample_key = meta.groupby("Sample", observed=True)[["Patient", "Tissue origin", "Platform"]].first().sort_index()
print("All embedded cells:", embedding.n_obs)
display(meta["Sample"].value_counts().reindex(sample_key.index).to_frame("whole_cells"))'''),
    md('''## 1. Conventional line A — cell-state composition

Each sample is represented by published cell-state proportions. A centered
log-ratio (CLR) transform handles the fact that proportions sum to one. Response
is not used to select states or construct the representation.'''),
    py('''sample_ids = sample_key.index.tolist()
valid_state = meta["celltype"].notna() & ~meta["celltype"].astype(str).isin(["AMB cells", "NA", "nan"])
counts = pd.crosstab(
    meta.loc[valid_state, "Sample"].astype(str),
    meta.loc[valid_state, "celltype"].astype(str),
).reindex(sample_ids, fill_value=0)
proportions = counts.div(counts.sum(axis=1), axis=0)
pseudocount = 0.5 / counts.sum(axis=1)
adjusted = proportions.add(pseudocount, axis=0)
adjusted = adjusted.div(adjusted.sum(axis=1), axis=0)
composition_clr = np.log(adjusted).sub(np.log(adjusted).mean(axis=1), axis=0)
composition_matrix = composition_clr.to_numpy(dtype=np.float32)
composition_matrix /= np.maximum(np.linalg.norm(composition_matrix, axis=1, keepdims=True), 1e-12)
display(proportions.round(3))
print("Composition features:", composition_matrix.shape)'''),
    md('''## 2. Conventional line B — CD4/CD8 functional programs

Predefined marker programs summarize Treg suppression, helper/memory, TH17,
cytotoxicity, exhaustion, interferon response and proliferation. Scores are
mean log-normalized expression within CD4/CD8 cells, aggregated per sample,
standardized across samples, and never selected using response.'''),
    py('''programs = {
    "Treg_suppression": ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
    "CD4_helper_memory": ["IL7R", "CCR7", "TCF7", "LTB"],
    "TH17": ["KLRB1", "CCR6", "RORA"],
    "cytotoxicity": ["NKG7", "CCL5", "PRF1", "GZMB", "GNLY"],
    "exhaustion": ["PDCD1", "TOX", "LAG3", "HAVCR2", "CXCL13"],
    "IFN_response": ["ISG15", "IFIT1", "IFIT3", "MX1"],
    "proliferation": ["MKI67", "TOP2A", "STMN1"],
}
cd48_ids = meta.index[meta["lineage.sub"].astype(str).isin(["CD4+ T cells", "CD8+ T cells"])]
raw = sc.read_h5ad(RAW_PATH)
available_genes = sorted({gene for genes in programs.values() for gene in genes if gene in raw.var_names})
functional_cells = raw[cd48_ids, available_genes].copy()
del raw
sc.pp.normalize_total(functional_cells, target_sum=1e4)
sc.pp.log1p(functional_cells)

cell_scores = pd.DataFrame(index=functional_cells.obs_names)
for program, genes in programs.items():
    genes = [gene for gene in genes if gene in functional_cells.var_names]
    cell_scores[program] = np.asarray(functional_cells[:, genes].X.mean(axis=1)).ravel()
cell_scores["Sample"] = meta.loc[cell_scores.index, "Sample"].astype(str)
functional_scores = cell_scores.groupby("Sample", observed=True)[list(programs)].mean().reindex(sample_ids)
functional_z = (functional_scores - functional_scores.mean()) / functional_scores.std(ddof=0).replace(0, 1)
functional_matrix = functional_z.to_numpy(dtype=np.float32)
functional_matrix /= np.maximum(np.linalg.norm(functional_matrix, axis=1, keepdims=True), 1e-12)
display(functional_scores.round(3))
print("Functional features:", functional_matrix.shape)'''),
    md('''## 3. Foundation-model line — pooled scGPT embeddings

Mean pooling treats cells equally. Cluster-aware pooling gives globally rarer
unsupervised scGPT states up to 3× weight. Both are calculated for all cells
and for the combined CD4/CD8 compartment.'''),
    py('''def pool(obj, metadata, sample_ids):
    work = obj.copy()
    sc.pp.neighbors(work, use_rep="X", n_neighbors=15, random_state=RANDOM_STATE)
    sc.tl.leiden(work, resolution=0.3, key_added="cluster", random_state=RANDOM_STATE,
                 flavor="igraph", n_iterations=2, directed=False)
    matrix = np.asarray(work.X, dtype=np.float32)
    clusters = work.obs["cluster"].astype(str)
    frequency = clusters.value_counts(normalize=True)
    cell_weight = clusters.map((1 / np.sqrt(frequency)).clip(upper=3.0)).to_numpy()
    means, weighted = [], []
    for sample in sample_ids:
        selected = metadata["Sample"].astype(str).to_numpy() == sample
        assert selected.sum() > 0
        means.append(matrix[selected].mean(axis=0))
        weighted.append(np.average(matrix[selected], axis=0, weights=cell_weight[selected]))
    def norm(x):
        x = np.asarray(x, dtype=np.float32)
        return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)
    return norm(means), norm(weighted)

whole_mean, whole_weighted = pool(embedding, meta, sample_ids)
cd48_mask = meta["lineage.sub"].astype(str).isin(["CD4+ T cells", "CD8+ T cells"]).to_numpy()
cd48 = embedding[cd48_mask].copy()
cd48_meta = meta.iloc[np.flatnonzero(cd48_mask)].copy()
coverage = cd48_meta["Sample"].value_counts().reindex(sample_ids).fillna(0).astype(int)
assert coverage.gt(0).all()
display(coverage.to_frame("CD4_CD8_cells"))
cd48_mean, cd48_weighted = pool(cd48, cd48_meta, sample_ids)

representations = {
    "composition_CLR": composition_matrix,
    "CD4_CD8_programs": functional_matrix,
    "whole_mean": whole_mean,
    "whole_cluster_aware": whole_weighted,
    "CD4_CD8_mean": cd48_mean,
    "CD4_CD8_cluster_aware": cd48_weighted,
}
print({key: value.shape for key, value in representations.items()})'''),
    md('''## 4. Reveal response and visualize all three evidence lines'''),
    py('''response = response_locked.groupby("Sample", observed=True)["Response"].first().reindex(sample_ids)
sample_key["Response"] = response.astype(str)
assert sample_key["Response"].value_counts().to_dict() == {"Non-responder": 7, "Responder": 3}
colors = {"Responder": "#2878B5", "Non-responder": "#D9534F"}
special_samples = set(sample_key.index[sample_key["Patient"].astype(str).eq("P1076")])
assert special_samples == {"EBUS_76", "NECK_05"}
titles = {
    "composition_CLR": "Conventional: cell-state composition (CLR)",
    "CD4_CD8_programs": "Conventional: CD4/CD8 functional programs",
    "whole_mean": "scGPT: whole cells — mean pooling",
    "whole_cluster_aware": "scGPT: whole cells — cluster-aware pooling",
    "CD4_CD8_mean": "scGPT: CD4/CD8 — mean pooling",
    "CD4_CD8_cluster_aware": "scGPT: CD4/CD8 — cluster-aware pooling",
}
umaps = {}
for name, matrix in representations.items():
    view = ad.AnnData(X=matrix, obs=sample_key.copy())
    sc.pp.neighbors(
        view, n_neighbors=3, use_rep="X", metric="cosine",
        random_state=RANDOM_STATE,
    )
    sc.tl.umap(view, random_state=RANDOM_STATE)
    umaps[name] = view.obsm["X_umap"].copy()
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for i, sample in enumerate(view.obs_names):
        response_group = view.obs.loc[sample, "Response"]
        marker = "^" if sample in special_samples else "o"
        ax.scatter(view.obsm["X_umap"][i, 0], view.obsm["X_umap"][i, 1],
                   s=105, color=colors[response_group], marker=marker,
                   edgecolor="white", linewidth=0.8)
        ax.annotate(sample, view.obsm["X_umap"][i], xytext=(4, 4), textcoords="offset points", fontsize=8)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["Responder"],
               markeredgecolor="white", markersize=9, label="Responder"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["Non-responder"],
               markeredgecolor="white", markersize=9, label="Non-responder"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="#777777",
               markeredgecolor="white", markersize=9, label="P1076 samples"),
    ]
    ax.set(title=titles[name], xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(handles=handles, frameon=False); fig.tight_layout()
    fig.savefig(FIGURES / f"sample_{name}_response_umap.png", dpi=240, bbox_inches="tight")
    plt.show()'''),
    md('''## 5. Distance evidence across the three lines'''),
    py('''def pair_mean(distance, left, right=None):
    left = np.flatnonzero(left)
    pairs = list(combinations(left, 2)) if right is None else [(i, j) for i in left for j in np.flatnonzero(right)]
    return float(np.mean([distance[i, j] for i, j in pairs]))

labels = sample_key["Response"].to_numpy()
rows = []
for name, matrix in representations.items():
    d = pairwise_distances(matrix, metric="cosine")
    r = labels == "Responder"; n = ~r
    centroids = np.vstack([matrix[r].mean(0), matrix[n].mean(0)])
    observed = silhouette_score(matrix, labels, metric="cosine")
    null = []
    for chosen in combinations(range(10), 3):
        perm = np.array(["Non-responder"] * 10, dtype=object); perm[list(chosen)] = "Responder"
        null.append(silhouette_score(matrix, perm, metric="cosine"))
    rows.append({"representation": name, "silhouette_cosine": observed,
                 "within_responder": pair_mean(d, r), "within_nonresponder": pair_mean(d, n),
                 "between_groups": pair_mean(d, r, n),
                 "centroid_distance": pairwise_distances(centroids, metric="cosine")[0, 1],
                 "exact_10sample_permutation_p": np.mean(np.asarray(null) >= observed),
                 "n_permutations": len(null)})
metrics = pd.DataFrame(rows).sort_values("silhouette_cosine", ascending=False)
display(metrics)
metrics.to_csv(TABLES / "sample_level_three_line_metrics.csv", index=False)

direction_rows = []
response_mask = sample_key["Response"].eq("Responder").to_numpy()
for line, frame in {
    "composition_CLR": composition_clr,
    "CD4_CD8_programs": functional_z,
}.items():
    for feature in frame.columns:
        responder_mean = float(frame.loc[response_mask, feature].mean())
        nonresponder_mean = float(frame.loc[~response_mask, feature].mean())
        direction_rows.append({
            "line": line, "feature": feature,
            "responder_mean": responder_mean,
            "nonresponder_mean": nonresponder_mean,
            "responder_minus_nonresponder": responder_mean - nonresponder_mean,
        })
feature_directions = pd.DataFrame(direction_rows)
display(feature_directions.loc[feature_directions["line"].eq("CD4_CD8_programs")]
        .sort_values("responder_minus_nonresponder", ascending=False))
display(feature_directions.loc[feature_directions["line"].eq("composition_CLR")]
        .assign(abs_difference=lambda x: x["responder_minus_nonresponder"].abs())
        .sort_values("abs_difference", ascending=False).head(12).drop(columns="abs_difference"))
feature_directions.to_csv(TABLES / "sample_level_feature_directions.csv", index=False)'''),
    md('''## 6. Save reusable sample representations'''),
    py('''out = ad.AnnData(X=whole_mean, obs=sample_key.copy())
out.var_names = embedding.var_names.copy()
out.obsm["X_whole_cluster_aware"] = whole_weighted
out.obsm["X_CD4_CD8_mean"] = cd48_mean
out.obsm["X_CD4_CD8_cluster_aware"] = cd48_weighted
out.obsm["X_composition_CLR"] = composition_matrix
out.obsm["X_CD4_CD8_programs"] = functional_matrix
out.uns["composition_features"] = composition_clr.columns.astype(str).tolist()
out.uns["functional_programs"] = list(programs)
for name, coords in umaps.items(): out.obsm[f"X_umap_{name}"] = coords
path = EMBEDDINGS / "GSE205335_10sample_sample_embeddings_v1.h5ad"
out.write_h5ad(path, compression="gzip")
print("Saved", path)'''),
    md('''## 7. Demo conclusion

Compare whether response is visible through cell composition, through
interpretable CD4/CD8 functions, and through frozen scGPT geometry. Agreement
across lines is stronger evidence than one attractive UMAP. With ten samples
this remains an explanatory demo, not a validated predictor.'''),
])
