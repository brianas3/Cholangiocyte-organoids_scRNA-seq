"""
Step 03 -- Normalize, score/regress cell cycle, batch-correct.

Reproduces Seurat LogNormalize -> FindVariableFeatures -> CellCycleScoring
(regress G2M-S diff only) -> ScaleData -> PCA, per the paper's Methods.

DEVIATION: the paper used scran::fastMNN for batch correction.
bioconductor-batchelor has no osx-arm64 conda build, and the pip fallback
'mnnpy' fails to compile (Apple clang lacks -fopenmp) on Apple Silicon.
Harmony (harmonypy) is substituted here -- see document/03_normalize.md.

Inputs:  outputs/objects/02_qc_filtered.h5ad
Outputs: outputs/objects/03_normalize_normalized_corrected.h5ad
"""
import json
import numpy as np
import scanpy as sc
import harmonypy

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"


def get_seurat_cc_genes():
    """Seurat's cc.genes$s.genes / $g2m.genes (Tirosh et al.), mapped to Ensembl IDs.
    Extracted once via: R -e 'library(Seurat); cc.genes$s.genes; cc.genes$g2m.genes'
    then bulk-mapped through https://rest.ensembl.org/lookup/symbol/homo_sapiens.
    Cached here as cc_genes_ensembl.json (see workspace) -- regenerate if missing.
    """
    with open(f"{ROOT}/data/metadata/cc_genes_ensembl.json") as f:
        cc = json.load(f)
    return list(cc["s_genes"].values()), list(cc["g2m_genes"].values())


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/02_qc_filtered.h5ad")
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    s_ids, g2m_ids = get_seurat_cc_genes()
    s_present = [g for g in s_ids if g in adata.var_names]
    g2m_present = [g for g in g2m_ids if g in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_present, g2m_genes=g2m_present)
    adata.obs["cc_diff"] = adata.obs["G2M_score"] - adata.obs["S_score"]

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.regress_out(adata_hvg, ["cc_diff"], n_jobs=8)
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

    data_mat = np.asarray(adata_hvg.obsm["X_pca"], dtype=np.float64)
    meta = adata_hvg.obs[["biosd_sample"]].copy()
    ho = harmonypy.run_harmony(data_mat, meta, ["biosd_sample"], max_iter_harmony=20, random_state=0)
    # NOTE: harmonypy 2.0.2 returns Z_corr already as (n_cells, n_pcs) -- do NOT
    # transpose (scanpy.external.pp.harmony_integrate's `.Z_corr.T` mis-shapes it
    # for this version; call harmonypy directly as done here).
    adata_hvg.obsm["X_pca_harmony"] = ho.Z_corr

    adata_hvg.write_h5ad(f"{ROOT}/outputs/objects/03_normalize_normalized_corrected.h5ad")
    print(f"Normalized+corrected: {adata_hvg.shape}")


if __name__ == "__main__":
    main()
