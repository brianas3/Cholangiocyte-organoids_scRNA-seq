"""
Step 04c -- Full-dataset UMAP WITHOUT Harmony (for Fig.2A).

Every biosd_sample belongs to exactly ONE origin x region category (22 samples, 9
categories), so the sample batch key is nested in the biology. Harmony on
`biosd_sample` (03_normalize.py) therefore erases category structure: kNN same-category
purity at 30 PCs, uncorrected PCA vs Harmony = PRI CBD 0.99 vs 0.75, PRI GB 0.99 vs 0.58,
PRI IHD 0.97 vs 0.38 (the paper's Fig.2A shows the three PRI regions as distinct lobes).

Instead, sample offsets are removed WITHIN each origin x region category only (same
approach as 04b): X_pca_wc = X_pca - sample_mean + category_mean. Origin and region
differences are untouched. Caveat: sample and category effects cannot be separated
statistically in this design.

Input:  outputs/objects/04_cluster_full_clustered.h5ad (uses its X_pca)
Output: outputs/objects/04c_full_nocorrect_umap.h5ad
        outputs/tables/04c_full_knn_category_purity.csv
"""
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
N_PCS = 30


def center_sample_within_category(adata, rep="X_pca"):
    X = np.asarray(adata.obsm[rep], dtype=np.float64)
    cat, smp = adata.obs["cat"].astype(str).values, adata.obs["biosd_sample"].astype(str).values
    out = X.copy()
    for c in np.unique(cat):
        mc = cat == c
        cat_mean = X[mc].mean(0)
        for s in np.unique(smp[mc]):
            ms = mc & (smp == s)
            out[ms] = X[ms] - X[ms].mean(0) + cat_mean
    return out


def knn_purity(adata, rep, n_pcs=N_PCS, k=15):
    codes = pd.Categorical(adata.obs["cat"]).codes
    onehot = np.eye(codes.max() + 1)[codes]
    tmp = sc.AnnData(np.zeros((adata.n_obs, 1)))
    tmp.obsm["X"] = adata.obsm[rep][:, :n_pcs]
    sc.pp.neighbors(tmp, use_rep="X", n_neighbors=k)
    nb = tmp.obsp["connectivities"].tocsr() @ onehot
    nb /= nb.sum(1, keepdims=True)
    own = nb[np.arange(adata.n_obs), codes]
    return pd.Series(own).groupby(adata.obs["cat"].values, observed=True).mean()


def main():
    a = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_full_clustered.h5ad")
    a.obs["cat"] = a.obs["origin"].astype(str) + " " + a.obs["region"].astype(str)
    a.obsm["X_pca_wc"] = center_sample_within_category(a)

    sc.pp.neighbors(a, use_rep="X_pca_wc", n_pcs=N_PCS, n_neighbors=15, random_state=0)
    sc.tl.umap(a, random_state=0)

    purity = pd.DataFrame({"pca_no_correction": knn_purity(a, "X_pca"),
                           "pca_sample_centered_within_category": knn_purity(a, "X_pca_wc"),
                           "harmony_sample": knn_purity(a, "X_pca_harmony")})
    purity.to_csv(f"{ROOT}/outputs/tables/04c_full_knn_category_purity.csv")
    print(purity.round(3))
    a.write_h5ad(f"{ROOT}/outputs/objects/04c_full_nocorrect_umap.h5ad")


if __name__ == "__main__":
    main()
