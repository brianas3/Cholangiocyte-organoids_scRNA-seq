"""
Step 06 -- Fig.2A analog: full-dataset (primary + organoid + bile-treated)
UMAP, colored by origin and by region side by side.

Input:  outputs/objects/04_cluster_full_clustered.h5ad
Output: outputs/figures/06_fig2a_full_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
REGION_COLORS = {"IHD": "#1b9e77", "CBD": "#d95f02", "GB": "#7570b3"}
ORIGIN_COLORS = {"PRI": "#e41a1c", "ORG": "#999999", "BTO": "#4daf4a"}


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_full_clustered.h5ad")
    emb, obs = adata.obsm["X_umap"], adata.obs

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6))

    ax = axes[0]
    for org, c in ORIGIN_COLORS.items():
        m = (obs["origin"] == org).values
        ax.scatter(emb[m, 0], emb[m, 1], s=2, alpha=0.35, color=c, linewidths=0, label=org, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(loc="upper left", frameon=False, markerscale=3, fontsize=8, bbox_to_anchor=(0.0, 1.0))
    ax.set_title("Colored by origin", fontsize=9)
    ax.margins(0.05)
    ax.annotate("", xy=(0.16, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=1.1, color="black"))
    ax.annotate("", xy=(0.02, 0.16), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=1.1, color="black"))
    ax.text(0.09, -0.015, "UMAP1", transform=ax.transAxes, fontsize=7, ha="center", va="top")
    ax.text(-0.015, 0.09, "UMAP2", transform=ax.transAxes, fontsize=7, ha="right", va="center", rotation=90)

    ax = axes[1]
    for reg, c in REGION_COLORS.items():
        m = (obs["region"] == reg).values
        ax.scatter(emb[m, 0], emb[m, 1], s=2, alpha=0.35, color=c, linewidths=0, label=reg, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(loc="upper left", frameon=False, markerscale=3, fontsize=8, bbox_to_anchor=(0.0, 1.0))
    ax.set_title("Colored by region", fontsize=9)
    ax.margins(0.05)

    fig.suptitle(f"Organoids converge across regions and partially re-diverge after bile treatment "
                 f"(n={adata.n_obs:,} cells)", fontsize=10, y=1.03)
    fig.tight_layout()

    fig.savefig(f"{ROOT}/outputs/figures/06_fig2a_full_umap.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/06_fig2a_full_umap.pdf")


if __name__ == "__main__":
    main()
