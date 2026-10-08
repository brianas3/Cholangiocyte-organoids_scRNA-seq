"""
Step 09 -- Supplementary Fig.S1A analog: full-dataset UMAP coloured by sample.

Paper (Sampaziotis 2021, Fig.S1A, supplement p.10): UMAP of all sequenced samples; each
patient / organoid line has a unique colour + marker. Marker shape by group: PRI CBD =
square, PRI GB = triangle-down, PRI IHD = circle, ORG = triangle-left, BTO = triangle-up.

Differences from the paper:
- "PRI IHD 5" (MacParland et al. 2018, cluster 17) is not in E-MTAB-8495 -> not drawn.
- Legend labels "<group> 1/2/3" are assigned by ascending biosd_sample ID within each
  origin x region group; the paper does not publish its own sample-to-label key, so the
  numbering is an assumption (colours are visual approximations, not the authors' values).
- Layout: UMAP from 04c (sample offsets removed WITHIN each origin x region category), so
  mixing of samples inside a category is partly by construction. Use for layout comparison,
  not as an independent batch-effect check.

Input:  outputs/objects/04c_full_nocorrect_umap.h5ad
Output: outputs/figures/09_figS1a_sample_umap.png / .pdf
"""
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from matplotlib.lines import Line2D

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"

# (origin, region) -> marker; colours per sample index, approximated by eye from Fig.S1A
MARKER = {"PRI CBD": "s", "PRI GB": "v", "PRI IHD": "o", "ORG": "<", "BTO": "^"}
COLORS = {
    ("PRI", "CBD"): ["#B4F065", "#0B5C5C", "#4F5B0E"],
    ("PRI", "GB"):  ["#F4645F", "#F28A1E", "#E81CE8"],
    ("PRI", "IHD"): ["#6FA8DC", "#F2D21A", "#8B3A12", "#B79CE0"],
    ("ORG", "CBD"): ["#6FE3C8", "#F4A0C8"],
    ("ORG", "GB"):  ["#C8F0B0", "#C8601A"],
    ("ORG", "IHD"): ["#C8C814", "#D0D0D0"],
    ("BTO", "CBD"): ["#DCC8F0", "#F4A0A8"],
    ("BTO", "GB"):  ["#FFFBC0", "#A050C8"],
    ("BTO", "IHD"): ["#C0DCF8", "#50E8F0"],
}
# legend order as in the paper
ORDER = [("PRI", "CBD"), ("PRI", "GB"), ("PRI", "IHD"),
         ("ORG", "CBD"), ("ORG", "GB"), ("ORG", "IHD"),
         ("BTO", "CBD"), ("BTO", "GB"), ("BTO", "IHD")]


def main():
    a = sc.read_h5ad(f"{ROOT}/outputs/objects/04c_full_nocorrect_umap.h5ad")
    emb, obs = a.obsm["X_umap"], a.obs
    rng = np.random.default_rng(0)

    # one row per sample: (label, marker, colour, cell indices)
    samples = []
    for o, r in ORDER:
        ids = sorted(obs.loc[(obs.origin == o) & (obs.region == r), "biosd_sample"].astype(str).unique())
        mk = MARKER[f"{o} {r}"] if o == "PRI" else MARKER[o]
        for i, sid in enumerate(ids):
            idx = np.where((obs.biosd_sample.astype(str) == sid).values)[0]
            samples.append((f"{o} {r} {i + 1}", mk, COLORS[(o, r)][i], idx))

    allidx = np.concatenate([s[3] for s in samples])
    lab = np.concatenate([[k] * len(s[3]) for k, s in enumerate(samples)])
    perm = rng.permutation(len(allidx))  # shuffle so no sample is systematically on top
    allidx, lab = allidx[perm], lab[perm]

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    for k, (name, mk, col, _) in enumerate(samples):
        m = lab == k
        ax.scatter(emb[allidx[m], 0], emb[allidx[m], 1], s=2.2, c=col, marker=mk,
                   linewidths=0, alpha=0.9)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
    handles = [Line2D([], [], marker=mk, ls="", color=col, markersize=5.5, label=name)
               for name, mk, col, _ in samples]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=6.8, labelspacing=0.45, handletextpad=0.3)
    fig.tight_layout()
    fig.savefig(f"{ROOT}/outputs/figures/09_figS1a_sample_umap.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{ROOT}/outputs/figures/09_figS1a_sample_umap.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
