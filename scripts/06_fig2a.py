"""
Step 06 -- Fig.2A analog: full-dataset (primary + organoid + bile-treated)
UMAP, colored by the 9 origin x region categories in a single panel, matching
Fig.2A of Sampaziotis et al. 2021 (Science 371:839-846).

Colors and layout approximated BY EYE from the published figure: boxed panel
(all 4 spines visible), black-edged filled points, one legend below the axes
laid out as a 3x3 grid (columns = region: IHD/CBD/GB, rows = origin:
PRI/ORG/BTO), plain "UMAP 1"/"UMAP 2" axis-label text, no panel title. Exact
hex values are not published by the authors -- this is a visual
approximation, not a machine-verified palette match.

Input:  outputs/objects/04_cluster_full_clustered.h5ad
Output: outputs/figures/06_fig2a_full_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"

CAT_COLORS = {
    "PRI IHD": "#D4A017", "ORG IHD": "#2166AC", "BTO IHD": "#E67E22",
    "PRI CBD": "#7C8C3C", "ORG CBD": "#2E7D32", "BTO CBD": "#C2A66B",
    "PRI GB":  "#A13D2B", "ORG GB":  "#7EC8E3", "BTO GB":  "#6E4B1D",
}
# legend layout matches the paper: 3 columns (IHD, CBD, GB) x 3 rows (PRI, ORG, BTO)
LEGEND_ORDER = ["PRI IHD", "ORG IHD", "BTO IHD",
                "PRI CBD", "ORG CBD", "BTO CBD",
                "PRI GB",  "ORG GB",  "BTO GB"]


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_full_clustered.h5ad")
    emb, obs = adata.obsm["X_umap"], adata.obs
    obs = obs.assign(cat=obs["origin"].astype(str) + " " + obs["region"].astype(str))

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    for cat in LEGEND_ORDER:
        m = (obs["cat"] == cat).values
        ax.scatter(emb[m, 0], emb[m, 1], s=3, c=CAT_COLORS[cat], edgecolors="black",
                   linewidths=0.05, alpha=0.85, label=cat, rasterized=True)

    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.margins(0.04)
    ax.set_xlabel("UMAP 1", fontsize=9, labelpad=4)
    ax.set_ylabel("UMAP 2", fontsize=9, labelpad=4)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False,
              markerscale=2.2, fontsize=7.5, handletextpad=0.4, labelspacing=0.6,
              columnspacing=1.2)
    fig.tight_layout()

    fig.savefig(f"{ROOT}/outputs/figures/06_fig2a_full_umap.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/06_fig2a_full_umap.pdf")


if __name__ == "__main__":
    main()
