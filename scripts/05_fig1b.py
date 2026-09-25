"""
Step 05 -- Fig.1B analog: primary cholangiocyte UMAP colored by region.

Colors and layout approximated BY EYE from Fig.1B of Sampaziotis et al. 2021
(Science 371:839-846): boxed panel (all 4 spines visible), black-edged filled
points, "PRI <region>" legend in the lower-right, plain "UMAP 1"/"UMAP 2"
axis-label text (no title on the panel itself, no corner arrows). Exact hex
values are not published by the authors -- these are a visual approximation,
not a machine-verified palette match.

Input:  outputs/objects/04_cluster_primary_only_umap.h5ad
Output: outputs/figures/05_fig1b_primary_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
REGION_COLORS = {
    "IHD": "#D4A017",  # golden/mustard yellow
    "CBD": "#7C8C3C",  # olive / yellow-green
    "GB":  "#A13D2B",  # dark red / maroon
}


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_primary_only_umap.h5ad")
    emb, obs = adata.obsm["X_umap"], adata.obs

    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    for reg in ["IHD", "CBD", "GB"]:
        m = (obs["region"] == reg).values
        ax.scatter(emb[m, 0], emb[m, 1], s=10, c=REGION_COLORS[reg], edgecolors="black",
                   linewidths=0.15, alpha=0.9, label=f"PRI {reg}", rasterized=True)

    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.margins(0.04)
    ax.set_xlabel("UMAP 1", fontsize=9, labelpad=4)
    ax.set_ylabel("UMAP 2", fontsize=9, labelpad=4)

    ax.legend(loc="lower right", frameon=False, markerscale=1.4, fontsize=8,
              handletextpad=0.4, labelspacing=0.5, borderaxespad=0.8)
    fig.tight_layout()

    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.pdf")


if __name__ == "__main__":
    main()
