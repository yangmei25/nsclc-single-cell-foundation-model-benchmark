# Final report — scGPT–Scanpy NSCLC single-cell demo

## Project status

The five-notebook demo is complete. It starts from the complete GSE205335
dataset, fixes a site-matched 10-sample cohort, builds conventional and scGPT
cell representations, benchmarks them at cell level, and compares three
sample-level evidence lines.

This is an independent reanalysis of a selected public-data subset. It was
designed to demonstrate a Scanpy-versus-scGPT cell-to-sample workflow, not to
reproduce the source paper's figures, full-cohort statistics, or complete set
of biological claims.

## Data

- Full reusable dataset: 96,505 cells × 33,714 genes
- Demo cohort: 10 core metastatic-LN samples
- Demo cells: 29,614
- Response groups: 3 responder and 7 non-responder samples
- Frozen scGPT output: 29,614 cells × 512 dimensions

## Cell-level question

**Question:** Does frozen scGPT preserve known cell identity better than a
conventional Scanpy PCA baseline on exactly the same cells?

| Scope | Representation | Label silhouette | Sample-held-out kNN accuracy | Leiden ARI | Leiden NMI | Platform silhouette |
|---|---|---:|---:|---:|---:|---:|
| Broad cells | Scanpy PCA | 0.107 | 0.832 | 0.234 | 0.503 | 0.107 |
| Broad cells | frozen scGPT | 0.064 | **0.923** | 0.234 | **0.510** | 0.147 |
| CD4 | Scanpy PCA | 0.013 | 0.445 | 0.176 | 0.273 | 0.238 |
| CD4 | frozen scGPT | -0.036 | **0.476** | **0.184** | **0.313** | 0.368 |
| CD8 | Scanpy PCA | **0.067** | **0.620** | **0.275** | **0.372** | 0.117 |
| CD8 | frozen scGPT | 0.027 | 0.580 | 0.233 | 0.356 | 0.241 |

### Cell-level conclusion

Frozen scGPT improves cross-sample recognition of broad cell identity and gives
a small CD4 kNN/NMI advantage. It is not universally superior: Scanpy PCA is
better for CD8 fine states and has lower platform separation. The foundation
model therefore adds transferable identity information, but it also retains
more 3p/5p technical structure.

## Sample-level question

**Question:** Which response-blind sample representation best separates the
three responder and seven non-responder samples?

| Evidence line | Representation | Cosine silhouette | Exact 3-vs-7 permutation p |
|---|---|---:|---:|
| Conventional function | CD4/CD8 programs | **0.258** | **0.033** |
| Foundation model | CD4/CD8 mean | 0.157 | 0.183 |
| Foundation model | CD4/CD8 cluster-aware | 0.152 | 0.208 |
| Foundation model | Whole-cell mean | 0.033 | 0.358 |
| Foundation model | Whole-cell cluster-aware | 0.029 | 0.383 |
| Conventional composition | Cell-state CLR | 0.025 | 0.375 |

### Biological direction

Responder samples show higher combined proliferation, exhaustion and
cytotoxicity programs. Non-responder samples show higher TH17 and CD4
helper/memory programs. Composition alone does not separate the groups well,
although non-responders trend toward higher follicular B, CD4 TRM and CD4 TH17,
while responders trend toward higher CD8 TEX, CD8 IFN-response and CD8 TRM.

### Sample-level conclusion

The strongest response structure is an interpretable CD4/CD8 functional-state
profile. CD4/CD8 scGPT pooling captures a weaker version of the pattern, while
whole-cell pooling dilutes it. Cluster-aware weighting does not improve on the
simple mean in this dataset.

The two samples from P1076 are retained as separate points and shown as
triangles in every sample UMAP; response is encoded by color.

## Main figures

- `results/figures/demo/cell_level_broad_pca_vs_scgpt.png`
- `results/figures/demo/cell_level_CD4_CD8_pca_vs_scgpt.png`
- `results/figures/demo/sample_composition_CLR_response_umap.png`
- `results/figures/demo/sample_CD4_CD8_programs_response_umap.png`
- `results/figures/demo/sample_whole_mean_response_umap.png`
- `results/figures/demo/sample_whole_cluster_aware_response_umap.png`
- `results/figures/demo/sample_CD4_CD8_mean_response_umap.png`
- `results/figures/demo/sample_CD4_CD8_cluster_aware_response_umap.png`

## Reproducible workflow

1. `notebooks/01_data_download_and_preparation.ipynb`
2. `notebooks/02_scanpy_pca_10_samples.ipynb`
3. `notebooks/03_scgpt_cell_embeddings_10_samples.ipynb`
4. `notebooks/04_cell_level_pca_vs_scgpt.ipynb`
5. `notebooks/05_sample_level_embeddings.ipynb`

Notebook 3 is the only GPU step. The complete T4 run took approximately 7
minutes 15 seconds. Notebooks 4 and 5 use saved embeddings and run locally.

## Interpretation boundary

This is a runnable explanatory demo, not a clinical predictor. Ten samples can
demonstrate a coherent analysis and an exact label-permutation result, but they
cannot establish external validity. The defensible claim is that functional
CD4/CD8 summaries are most informative here, with frozen scGPT providing
complementary but weaker sample-level structure.
