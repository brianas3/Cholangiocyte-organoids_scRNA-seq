"""
Step 04b -- Primary-only UMAP WITHOUT Harmony, plus small-cluster removal.

Why: in the primary tissue each donor contributes exactly one region
(CBD = Donor 1-3, GB = 4-6, IHD = 7-10; biosd_sample is 1:1 with individual), so
region is fully nested in the batch key. Harmony on `individual` (04_cluster.py)
therefore removes the region signal (kNN same-region purity at 30 PCs:
raw PCA 0.99 vs Harmony 0.89; IHD 0.96 vs 0.52). The paper's fastMNN only
corrects mutual nearest neighbours and keeps most regional structure.
Caveat: with this design donor effect and region effect cannot be separated
statistically, uncorrected or corrected.

To keep one region as ONE cluster (donor split inside CBD: Donor 3 vs Donor 1+2), donor
effects are removed WITHIN each region only: every donor's PCA centroid is shifted to its
region's centroid (X_pca_wr = X_pca - donor_mean + region_mean). Between-region differences
are untouched because the shift never crosses regions.

Also does the paper's post-hoc removal of small clusters (fig. S3-S5): Leiden at a
fine resolution, drop clusters < MIN_FRAC of cells, then recompute neighbours/UMAP.

Input:  outputs/objects/04_cluster_primary_only_umap.h5ad (uses its X_pca)
Output: outputs/objects/04b_primary_nocorrect_umap.h5ad
        outputs/tables/04b_primary_knn_region_purity.csv
"""
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
N_PCS = 30
LEIDEN_RES = 1.0
MIN_FRAC = 0.01  # drop clusters smaller than 1% of primary cells


def knn_purity(adata, rep, n_pcs=N_PCS, k=15):
    """Mean fraction of each cell's kNN-graph neighbours sharing its region."""
    codes = pd.Categorical(adata.obs["region"]).codes
    onehot = np.eye(codes.max() + 1)[codes]
    tmp = sc.AnnData(np.zeros((adata.n_obs, 1)))
    tmp.obsm["X"] = adata.obsm[rep][:, :n_pcs]
    sc.pp.neighbors(tmp, use_rep="X", n_neighbors=k)
    nb = tmp.obsp["connectivities"].tocsr() @ onehot
    nb /= nb.sum(1, keepdims=True)
    own = nb[np.arange(adata.n_obs), codes]
    return pd.Series(own).groupby(adata.obs["region"].values, observed=True).mean()


def center_donor_within_region(adata, rep="X_pca"):
    """Shift each donor's centroid onto its region centroid (donors nested in region)."""
    X = np.asarray(adata.obsm[rep], dtype=np.float64).copy()
    reg, ind = adata.obs["region"].astype(str).values, adata.obs["individual"].astype(str).values
    out = X.copy()
    for r in np.unique(reg):
        mr = reg == r
        region_mean = X[mr].mean(0)
        for d in np.unique(ind[mr]):
            md = mr & (ind == d)
            out[md] = X[md] - X[md].mean(0) + region_mean
    return out


def main():
    a = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_primary_only_umap.h5ad")
    n0 = a.n_obs
    a.obsm["X_pca_wr"] = center_donor_within_region(a)
    REP = "X_pca_wr"

    sc.pp.neighbors(a, use_rep=REP, n_pcs=N_PCS, n_neighbors=15, random_state=0)
    sc.tl.leiden(a, resolution=LEIDEN_RES, random_state=0, key_added="leiden_fine",
                 flavor="igraph", n_iterations=2, directed=False)
    sizes = a.obs["leiden_fine"].value_counts()
    small = sizes[sizes < MIN_FRAC * n0].index
    print(f"{len(sizes)} clusters, dropping {len(small)} small: "
          f"{sizes[small].to_dict()} ({sizes[small].sum()} cells)")
    a = a[~a.obs["leiden_fine"].isin(small)].copy()

    sc.pp.neighbors(a, use_rep=REP, n_pcs=N_PCS, n_neighbors=15, random_state=0)
    sc.tl.umap(a, random_state=0)

    purity = pd.DataFrame({"pca_no_correction": knn_purity(a, "X_pca"),
                           "pca_donor_centered_within_region": knn_purity(a, REP),
                           "harmony_individual": knn_purity(a, "X_pca_harmony")})
    purity.to_csv(f"{ROOT}/outputs/tables/04b_primary_knn_region_purity.csv")
    print(purity.round(3))

    a.write_h5ad(f"{ROOT}/outputs/objects/04b_primary_nocorrect_umap.h5ad")
    print(f"{n0} -> {a.n_obs} cells")


if __name__ == "__main__":
    main()
