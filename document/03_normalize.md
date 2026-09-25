# 03_normalize — Normalization, cell-cycle regression, batch correction (Harmony substitute for fastMNN)

**Question** Reproduce the paper's Seurat normalize -> HVG -> cell-cycle regression -> batch-correction stage on the QC-filtered cells.

**Script** `scripts/03_normalize.py`

## Inputs
- outputs/objects/02_qc_filtered.h5ad

## What was done
- sc.pp.normalize_total(target_sum=1e4) + log1p, equivalent to Seurat's LogNormalize(scale.factor=10000).
- Cell-cycle scoring via sc.tl.score_genes_cell_cycle using Seurat's own cc.genes S/G2M lists (43 S genes, 54 G2M genes, extracted directly from the R Seurat package and mapped to Ensembl IDs; all present after the gene filter).
- 2000 HVGs selected (Seurat flavor), matching Seurat's FindVariableFeatures default n=2000.
- sc.pp.regress_out on cc_diff = G2M_score - S_score ONLY (per paper's stated 'regressing out the difference between the G2M and S phase scores'), then sc.pp.scale(max_value=10) (Seurat ScaleData default clip), then PCA (50 PCs, arpack).
- DEVIATION: paper used scran::fastMNN for batch correction. bioconductor-batchelor has no osx-arm64 conda build and the pip fallback mnnpy fails to compile (Apple clang lacks -fopenmp) on this machine. Substituted Harmony (harmonypy 2.0.2, run directly -- the scanpy.external wrapper mis-transposes Z_corr in this harmonypy version) with batch key = biosd_sample (22 batches, converged in 4 iterations), producing a 50-dim corrected embedding analogous to fastMNN's corrected PCs.

## Outputs
- outputs/objects/03_normalize_normalized_corrected.h5ad (40,732 cells x 2,000 HVGs, obsm['X_pca'] and obsm['X_pca_harmony'])

## Findings
- Harmony batch correction converged after 4 of a max 20 iterations across the 22 sample batches.
- All 43 S-phase and 54 G2M Seurat cell-cycle genes survived the earlier gene filter, so cell-cycle scoring used the full canonical Tirosh et al. gene lists with no substitutions.

## Still open
- Batch correction method (Harmony) differs from the paper's fastMNN -- both are MNN/graph-based batch-mixing methods but are not numerically identical; downstream UMAP topology should be compared qualitatively to Fig. 2A, not expected to match exactly.

