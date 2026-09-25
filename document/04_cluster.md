# 04_cluster — Louvain clustering + UMAP (full dataset and primary-only)

**Question** Reproduce the paper's Louvain-clustering fidelity check (ARI/AMI of clusters vs origin/region) and compute UMAP embeddings for Fig.1B (primary-only) and Fig.2A (full dataset) analogs.

**Script** `scripts/04_cluster.py`

## Inputs
- outputs/objects/03_normalize_normalized_corrected.h5ad (full dataset, 40,732 cells)
- outputs/objects/02_qc_filtered.h5ad (subset to PRI for the primary-only pipeline)

## What was done
- Primary-only pipeline (independent of the full-dataset one, to match Fig.1B's own reported scope): QC-passed PRI cells (7,316) -> normalize_total+log1p -> cell-cycle score/regress cc_diff -> 2000 HVGs -> scale -> PCA(50) -> Harmony(batch=individual, 10 batches) -> neighbors(k=15) -> UMAP.
- Full-dataset pipeline (all origins, 40,732 cells, from step 03) -> neighbors(k=15) on X_pca_harmony -> UMAP.
- Louvain resolution scanned on the full dataset (sc.tl.louvain, 'louvain' package via igraph) over [0.005-1.0]; resolution=0.05 is the first value producing exactly 3 clusters, matching the paper's target of 3 clusters (paper implicitly targets the 3 origins: PRI/ORG/BTO).
- ARI/AMI between the 3 Louvain clusters and (a) origin and (b) region labels, averaged over 10 random seeds at resolution=0.05.

## Outputs
- outputs/objects/04_cluster_full_clustered.h5ad
- outputs/objects/04_cluster_primary_only_umap.h5ad
- outputs/tables/04_cluster_louvain_ari_ami.csv

## Findings
- Mean ARI(origin)=0.652, AMI(origin)=0.586 vs mean ARI(region)=0.014, AMI(region)=0.016 across 10 seeds -- same QUALITATIVE conclusion as the paper (clusters track origin, not region), but the paper reports origin ARI/AMI>0.95 and region<0.3, i.e. our origin correspondence is considerably weaker in absolute terms.
- The gap is attributable to using Harmony instead of fastMNN for batch correction, the Atlas-requantified (vs authors' own CellRanger) counts, and the resulting different cell composition after QC (see document/02_qc.md).
- Primary-only subset: 7,316 cells across 10 individuals, close to the paper's reported 7,295 cells / n=10 individuals for Fig.1B.

## Still open
- Paper's post-hoc removal of small mito/gene-count-outlier Louvain clusters was not applied here as a separate iteration -- the isOutlier-based cell QC in step 02 already removes most such cells before clustering.

