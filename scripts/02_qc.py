"""
Step 02 -- QC filtering of the Atlas-hosted E-MTAB-8495 raw count matrix,
reproducing the gene/cell QC criteria described in Sampaziotis et al. 2021
(Science 371:839-846) Methods.

Inputs:
  data/Result files/E-MTAB-8495-quantification-raw-files/*.mtx(+cols/rows)
  data/E-MTAB-8495_biosample_to_sample.csv

Outputs:
  outputs/tables/02_qc_cells_per_sample.csv
  outputs/figures/02_qc_diagnostics.png / .pdf
  outputs/objects/02_qc_filtered.h5ad
"""
import json, requests, time
import numpy as np
import pandas as pd
import scipy.io as sio
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
DATA = f"{ROOT}/data/Result files/E-MTAB-8495-quantification-raw-files"
MANIFEST = f"{ROOT}/data/E-MTAB-8495_biosample_to_sample.csv"


def load_raw_adata():
    raw = sio.mmread(f"{DATA}/E-MTAB-8495.aggregated_filtered_counts.mtx").tocsr()
    rows = pd.read_csv(f"{DATA}/E-MTAB-8495.aggregated_filtered_counts.mtx_rows",
                        sep="\t", header=None, names=["ensembl_id", "ensembl_id2"])
    cols = pd.read_csv(f"{DATA}/E-MTAB-8495.aggregated_filtered_counts.mtx_cols",
                        header=None, names=["barcode"])
    import anndata as ad
    adata = ad.AnnData(X=raw.T.tocsr())
    adata.obs_names = cols["barcode"].values.astype(str)
    adata.var_names = rows["ensembl_id"].values.astype(str)
    adata.obs["biosd_sample"] = adata.obs_names.str.split("-").str[0]

    manifest = pd.read_csv(MANIFEST)
    adata.obs = (adata.obs.reset_index().rename(columns={"index": "barcode_full"})
                 .merge(manifest, on="biosd_sample", how="left"))
    adata.obs.index = adata.obs_names = cols["barcode"].values.astype(str)

    region_map = {"intrahepatic bile duct": "IHD", "common bile duct": "CBD", "gall bladder": "GB"}
    adata.obs["region"] = adata.obs["organism_part"].map(region_map)
    adata.obs["origin"] = np.select_ = None  # placeholder, replaced below
    def origin(row):
        if row["growth_condition"] == "primary tissue":
            return "PRI"
        return "BTO" if row["stimulus"] == "bile" else "ORG"
    adata.obs["origin"] = adata.obs.apply(origin, axis=1)
    return adata


def fetch_mito_and_marker_ids():
    """Ensembl chrMT gene ids + EPCAM/KRT7/KRT19 ids (network calls to rest.ensembl.org)."""
    r = requests.get("https://rest.ensembl.org/overlap/region/human/MT:1-16569",
                      params={"feature": "gene", "content-type": "application/json"}, timeout=60)
    mt_ids = [g["gene_id"] for g in r.json()]
    marker_ids = {}
    for sym in ["EPCAM", "KRT7", "KRT19"]:
        for _ in range(5):
            rr = requests.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}",
                               params={"content-type": "application/json"}, timeout=30)
            if rr.ok:
                marker_ids[sym] = rr.json()["id"]
                break
            time.sleep(3)
    return mt_ids, marker_ids


def is_outlier_mad(x, nmads=3, kind="lower"):
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    lo, hi = med - nmads * mad, med + nmads * mad
    if kind == "lower":
        return x < lo
    if kind == "higher":
        return x > hi
    return (x < lo) | (x > hi)


def main():
    adata = load_raw_adata()
    mt_ids, marker_ids = fetch_mito_and_marker_ids()

    adata.var["mt"] = adata.var_names.isin(mt_ids)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

    for sym, gid in marker_ids.items():
        adata.obs[f"counts_{sym}"] = np.asarray(adata[:, gid].X.todense()).ravel() if gid in adata.var_names else 0
    adata.obs["biliary_marker_pos"] = (
        (adata.obs["counts_EPCAM"] > 3) | (adata.obs["counts_KRT7"] > 3) | (adata.obs["counts_KRT19"] > 3))

    adata.obs["log10_total_counts"] = np.log10(adata.obs["total_counts"])
    adata.obs["log10_n_genes"] = np.log10(adata.obs["n_genes_by_counts"])

    low_lib = np.zeros(adata.n_obs, dtype=bool)
    low_genes = np.zeros(adata.n_obs, dtype=bool)
    high_mito = np.zeros(adata.n_obs, dtype=bool)
    for _, idx in adata.obs.groupby("biosd_sample").groups.items():
        mask = adata.obs_names.isin(idx)
        low_lib[mask] = is_outlier_mad(adata.obs.loc[mask, "log10_total_counts"], kind="lower")
        low_genes[mask] = is_outlier_mad(adata.obs.loc[mask, "log10_n_genes"], kind="lower")
        high_mito[mask] = is_outlier_mad(adata.obs.loc[mask, "pct_counts_mt"], kind="higher")
    adata.obs["qc_outlier"] = low_lib | low_genes | high_mito
    adata.obs["pass_qc"] = (~adata.obs["qc_outlier"]) & adata.obs["biliary_marker_pos"]

    # Gene filter: >0 counts in >=3 cells per batch, in every batch of >=1 origin
    X = adata.X.tocsc()
    gt0 = X > 0
    origins, samples = adata.obs["origin"].values, adata.obs["biosd_sample"].values
    gene_pass = np.zeros(adata.n_vars, dtype=bool)
    for org in ["PRI", "ORG", "BTO"]:
        org_mask = origins == org
        pass_all = np.ones(adata.n_vars, dtype=bool)
        for b in np.unique(samples[org_mask]):
            bmask = org_mask & (samples == b)
            n_pos = np.asarray(gt0[bmask, :].sum(axis=0)).ravel()
            pass_all &= (n_pos >= 3)
        gene_pass |= pass_all
    adata.var["gene_filter_pass"] = gene_pass

    adata_qc = adata[adata.obs["pass_qc"].values, adata.var["gene_filter_pass"].values].copy()
    adata_qc.write_h5ad(f"{ROOT}/outputs/objects/02_qc_filtered.h5ad")

    paper_counts = {("PRI", "IHD"): 587, ("PRI", "CBD"): 3006, ("PRI", "GB"): 3702,
                    ("ORG", "IHD"): 6641, ("ORG", "CBD"): 5321, ("ORG", "GB"): 5859,
                    ("BTO", "IHD"): 3653, ("BTO", "CBD"): 3224, ("BTO", "GB"): 3815}
    qc_counts = adata.obs[adata.obs["pass_qc"]].groupby(["origin", "region"]).size()
    rows = [{"origin": o, "region": r, "n_pass_qc_this_pipeline": int(qc_counts.get((o, r), 0)),
             "n_paper_reported": n} for (o, r), n in paper_counts.items()]
    pd.DataFrame(rows).to_csv(f"{ROOT}/outputs/tables/02_qc_cells_per_sample.csv", index=False)

    print(f"QC-passed: {adata_qc.n_obs} cells x {adata_qc.n_vars} genes")


if __name__ == "__main__":
    main()
