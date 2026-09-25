"""
Step 04 -- Louvain clustering + UMAP, on both the full dataset and an
independent primary-tissue-only pipeline (to match Fig.1B's own reported
scope: 7,295 cells / 10 individuals in the paper).

Inputs:
  outputs/objects/02_qc_filtered.h5ad   (subset to origin == "PRI")
  outputs/objects/03_normalize_normalized_corrected.h5ad  (full dataset)

Outputs:
  outputs/objects/04_cluster_primary_only_umap.h5ad
  outputs/objects/04_cluster_full_clustered.h5ad
  outputs/tables/04_cluster_louvain_ari_ami.csv
"""
import json
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"


def load_cc_genes(adata_var_names):
    """Seurat cc.genes S/G2M lists, Ensembl-mapped (see scripts/03_normalize.py); restricted
    to genes present in this object's var_names (all pass in practice -- see document/03_normalize.md)."""
    with open(f"{ROOT}/data/metadata/cc_genes_ensembl.json") as f:
        cc = json.load(f)
    s_ids, g2m_ids = list(cc["s_genes"].values()), list(cc["g2m_genes"].values())
    return ([g for g in s_ids if g in adata_var_names],
            [g for g in g2m_ids if g in adata_var_names])


def harmony_pca(pca, batch_series, batch_key):
    meta = batch_series.to_frame(batch_key)
    ho = harmonypy.run_harmony(np.asarray(pca, dtype=np.float64), meta, [batch_key],
                                max_iter_harmony=20, random_state=0)
    return ho.Z_corr  # (n_cells, n_pcs) in harmonypy 2.0.2 -- no transpose needed


def build_primary_only(s_present, g2m_present):
    qc = sc.read_h5ad(f"{ROOT}/outputs/objects/02_qc_filtered.h5ad")
    pri = qc[qc.obs["origin"] == "PRI"].copy()
    pri.layers["counts"] = pri.X.copy()
    sc.pp.normalize_total(pri, target_sum=1e4)
    sc.pp.log1p(pri)
    sc.tl.score_genes_cell_cycle(pri, s_genes=s_present, g2m_genes=g2m_present)
    pri.obs["cc_diff"] = pri.obs["G2M_score"] - pri.obs["S_score"]
    sc.pp.highly_variable_genes(pri, n_top_genes=2000, flavor="seurat")
    pri_hvg = pri[:, pri.var["highly_variable"]].copy()
    sc.pp.regress_out(pri_hvg, ["cc_diff"], n_jobs=8)
    sc.pp.scale(pri_hvg, max_value=10)
    sc.tl.pca(pri_hvg, n_comps=50, svd_solver="arpack")
    pri_hvg.obsm["X_pca_harmony"] = harmony_pca(pri_hvg.obsm["X_pca"], pri_hvg.obs["individual"], "individual")
    sc.pp.neighbors(pri_hvg, use_rep="X_pca_harmony", n_neighbors=15, random_state=0)
    sc.tl.umap(pri_hvg, random_state=0)
    return pri_hvg


def scan_louvain_resolution(adata, target_k=3, grid=(0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.0)):
    for res in grid:
        sc.tl.louvain(adata, resolution=res, random_state=0, key_added="louvain_tmp")
        if adata.obs["louvain_tmp"].nunique() == target_k:
            return res
    raise RuntimeError(f"No resolution in grid gave exactly {target_k} clusters")


def main():
    qc_genes_probe = sc.read_h5ad(f"{ROOT}/outputs/objects/02_qc_filtered.h5ad").var_names
    s_present, g2m_present = load_cc_genes(qc_genes_probe)

    pri_hvg = build_primary_only(s_present, g2m_present)
    pri_hvg.write_h5ad(f"{ROOT}/outputs/objects/04_cluster_primary_only_umap.h5ad")

    adata_hvg = sc.read_h5ad(f"{ROOT}/outputs/objects/03_normalize_normalized_corrected.h5ad")
    sc.pp.neighbors(adata_hvg, use_rep="X_pca_harmony", n_neighbors=15, random_state=0)
    sc.tl.umap(adata_hvg, random_state=0)

    res = scan_louvain_resolution(adata_hvg, target_k=3)
    rows = []
    for seed in range(10):
        sc.tl.louvain(adata_hvg, resolution=res, random_state=seed, key_added=f"louvain_s{seed}")
        labels = adata_hvg.obs[f"louvain_s{seed}"]
        rows.append({
            "seed": seed, "n_clusters": labels.nunique(),
            "ARI_origin": adjusted_rand_score(adata_hvg.obs["origin"], labels),
            "AMI_origin": adjusted_mutual_info_score(adata_hvg.obs["origin"], labels),
            "ARI_region": adjusted_rand_score(adata_hvg.obs["region"], labels),
            "AMI_region": adjusted_mutual_info_score(adata_hvg.obs["region"], labels),
        })
        del adata_hvg.obs[f"louvain_s{seed}"]
    pd.DataFrame(rows).to_csv(f"{ROOT}/outputs/tables/04_cluster_louvain_ari_ami.csv", index=False)

    adata_hvg.write_h5ad(f"{ROOT}/outputs/objects/04_cluster_full_clustered.h5ad")
    print(f"Louvain resolution for {3} clusters: {res}")


if __name__ == "__main__":
    main()
