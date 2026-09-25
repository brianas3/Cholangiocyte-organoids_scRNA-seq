"""
Step 05 -- Fig.1B analog: primary cholangiocyte UMAP colored by region.

Input:  outputs/objects/04_cluster_primary_only_umap.h5ad
Output: outputs/figures/05_fig1b_primary_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
REGION_COLORS = {"IHD": "#1b9e77", "CBD": "#d95f02", "GB": "#7570b3"}


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/04_cluster_primary_only_umap.h5ad")
    emb, obs = adata.obsm["X_umap"], adata.obs

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    for reg, c in REGION_COLORS.items():
        m = (obs["region"] == reg).values
        ax.scatter(emb[m, 0], emb[m, 1], s=4, alpha=0.6, color=c, linewidths=0, label=reg, rasterized=True)

    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.margins(0.06)

    ax.annotate("", xy=(0.16, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=1.1, color="black"))
    ax.annotate("", xy=(0.02, 0.16), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=1.1, color="black"))
    ax.text(0.09, -0.015, "UMAP1", transform=ax.transAxes, fontsize=7, ha="center", va="top")
    ax.text(-0.015, 0.09, "UMAP2", transform=ax.transAxes, fontsize=7, ha="right", va="center", rotation=90)

    ax.legend(loc="upper left", frameon=False, markerscale=2, fontsize=8, labelspacing=0.6,
              bbox_to_anchor=(0.02, 0.98))
    ax.set_title(f"Primary cholangiocytes retain region-specific transcriptional identity\n"
                 f"(n={adata.n_obs:,} cells, {obs['individual'].nunique()} individuals)", fontsize=9)
    fig.tight_layout()

    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.pdf")


if __name__ == "__main__":
    main()
