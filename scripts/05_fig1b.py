"""
Step 05 -- Fig.1B analog: primary cholangiocyte UMAP colored by region.

Colors and layout approximated BY EYE from Fig.1B of Sampaziotis et al. 2021
(Science 371:839-846): boxed panel (all 4 spines visible), black-edged filled
points, per-region donor-shade legend (upper-left), plain "UMAP 1"/"UMAP 2"
axis-label text (no title on the panel itself, no corner arrows). Exact hex
values are not published by the authors -- these are a visual approximation,
not a machine-verified palette match.

Input:  outputs/objects/04b_primary_nocorrect_umap.h5ad
Output: outputs/figures/05_fig1b_primary_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
# One hue per region, one shade per donor (dark -> light by donor number), as in the
# paper's legend: PRI IHD 4 dots (Donor 7-10), PRI CBD 3 (Donor 1-3), PRI GB 3 (Donor 4-6).
# The paper's caption does not state that shades are donors -- inferred from the 4/3/3
# dot counts matching the 10 individuals. Shades are a visual approximation.
DONOR_SHADES = {
    "IHD": ["#8A6A0C", "#B8860B", "#DDAA1F", "#F5CF3E"],
    "CBD": ["#4F5F22", "#7C8C3C", "#A7B85A"],
    "GB":  ["#7A1F14", "#A13D2B", "#D2574A"],
}


def donor_colors(obs, reg):
    donors = sorted(obs.loc[obs["region"] == reg, "individual"].astype(str).unique(),
                    key=lambda d: int(d.split()[-1]))
    return dict(zip(donors, DONOR_SHADES[reg]))


def main():
    adata = sc.read_h5ad(f"{ROOT}/outputs/objects/04b_primary_nocorrect_umap.h5ad")
    emb, obs = adata.obsm["X_umap"], adata.obs

    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    rng = np.random.default_rng(0)
    legend_rows = []
    for reg in ["IHD", "CBD", "GB"]:
        cmap = donor_colors(obs, reg)
        legend_rows.append((reg, list(cmap.values())))
        idx = np.where((obs["region"] == reg).values)[0]
        idx = rng.permutation(idx)  # shuffle donors so no donor is drawn on top
        cols = obs["individual"].astype(str).iloc[idx].map(cmap).values
        ax.scatter(emb[idx, 0], emb[idx, 1], s=10, c=cols, edgecolors="black",
                   linewidths=0.15, alpha=0.9)

    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.margins(0.04)
    ax.set_xlabel("UMAP 1", fontsize=9, labelpad=4)
    ax.set_ylabel("UMAP 2", fontsize=9, labelpad=4)

    # legend (upper-left): per-region row of donor-shade dots, then "PRI <region>" label.
    # Dots left-aligned, label placed after the longest row (4 dots) so labels line up.
    x0, dx, y0, dy = 0.04, 0.045, 0.95, 0.065
    for i, (reg, cols) in enumerate(legend_rows):
        y = y0 - i * dy
        for j, c in enumerate(cols):
            ax.scatter([x0 + j * dx], [y], s=22, c=c, edgecolors="black", linewidths=0.4,
                       transform=ax.transAxes, zorder=10, clip_on=False)
        ax.text(x0 + 4 * dx, y, f"PRI {reg}", transform=ax.transAxes, fontsize=8,
                va="center", ha="left", zorder=10)
    fig.tight_layout()

    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.png", dpi=600)
    fig.savefig(f"{ROOT}/outputs/figures/05_fig1b_primary_umap.pdf")


if __name__ == "__main__":
    main()
