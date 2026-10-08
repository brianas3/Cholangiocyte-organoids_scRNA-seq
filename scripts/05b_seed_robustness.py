"""
Step 05b -- UMAP seed robustness for Fig.1B.

Re-runs sc.tl.umap on the SAME kNN graph (04b, donor-centered-within-region PCA) with
seeds 0-4. Seed only changes the 2-D layout, not the graph or clusters. Reports the
kNN same-region purity measured IN THE UMAP EMBEDDING for each seed (mean +- sd), so
the claim "regions separate regardless of seed" is quantified, not eyeballed.

Input:  outputs/objects/04b_primary_nocorrect_umap.h5ad
Output: outputs/figures/05b_seed_robustness.png / .pdf
        outputs/tables/05b_seed_robustness_purity.csv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
SEEDS = [0, 1, 2, 3, 4]
COLORS = {"IHD": "#D4A017", "CBD": "#7C8C3C", "GB": "#A13D2B"}


def umap_region_purity(emb, region, k=15):
    codes = pd.Categorical(region).codes
    onehot = np.eye(codes.max() + 1)[codes]
    tmp = sc.AnnData(np.zeros((len(region), 1)))
    tmp.obsm["X"] = emb
    sc.pp.neighbors(tmp, use_rep="X", n_neighbors=k)
    nb = tmp.obsp["connectivities"].tocsr() @ onehot
    nb /= nb.sum(1, keepdims=True)
    own = nb[np.arange(len(region)), codes]
    return pd.Series(own).groupby(np.asarray(region), observed=True).mean()


def main():
    a = sc.read_h5ad(f"{ROOT}/outputs/objects/04b_primary_nocorrect_umap.h5ad")
    region = a.obs["region"].astype(str).values
    rows, embs = [], {}
    for s in SEEDS:
        sc.tl.umap(a, random_state=s)
        embs[s] = a.obsm["X_umap"].copy()
        rows.append(umap_region_purity(embs[s], region).rename(s))
    tab = pd.DataFrame(rows)
    tab.index.name = "seed"
    tab.loc["mean"], tab.loc["sd"] = tab.iloc[:len(SEEDS)].mean(), tab.iloc[:len(SEEDS)].std()
    tab.round(3).to_csv(f"{ROOT}/outputs/tables/05b_seed_robustness_purity.csv")
    print(tab.round(3))

    fig, axes = plt.subplots(2, 3, figsize=(10, 6.4))
    rng = np.random.default_rng(0)
    order = rng.permutation(len(region))
    for ax, s in zip(axes.ravel(), SEEDS):
        e = embs[s][order]
        ax.scatter(e[:, 0], e[:, 1], s=2.5, c=[COLORS[r] for r in region[order]],
                   edgecolors="black", linewidths=0.05, alpha=0.9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"seed {s}", fontsize=9)
    ax = axes.ravel()[-1]
    ax.axis("off")
    for i, r in enumerate(["IHD", "CBD", "GB"]):
        ax.scatter([0.1], [0.8 - i * 0.12], s=40, c=COLORS[r], edgecolors="black",
                   linewidths=0.4, transform=ax.transAxes)
        ax.text(0.18, 0.8 - i * 0.12, f"PRI {r}", transform=ax.transAxes, va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{ROOT}/outputs/figures/05b_seed_robustness.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/05b_seed_robustness.pdf")


if __name__ == "__main__":
    main()
