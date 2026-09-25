# 06_fig2a — Fig. 2A analog: full-dataset UMAP (primary + organoid + bile-treated)

**Question** Reproduce Fig.2A's depiction of organoids converging across regions and bile-treated organoids shifting back toward a region-specific (gallbladder-associated) identity.

**Script** `scripts/06_fig2a.py`

## Inputs
- outputs/objects/04_cluster_full_clustered.h5ad (40,732 cells, Harmony-corrected UMAP)

## What was done
- Two-panel UMAP of all QC-passed cells (PRI+ORG+BTO): left colored by origin, right by region, using the same region color key as Fig.1B analog for cross-figure consistency.
- figure-style embedding conventions applied: no tick marks, small corner UMAP1/UMAP2 arrow indicator, frameless legend.

## Outputs
- outputs/figures/06_fig2a_full_umap.png (+ .pdf)

## Findings
- Origin panel: PRI, ORG and BTO occupy three largely distinct zones (PRI bottom, ORG large cluster at right, BTO upper-left, adjoining PRI) -- reproducing the paper's central claim that organoid culture shifts cells into a shared transcriptional state distinct from primary tissue, and that bile treatment shifts organoids back toward the primary-tissue neighborhood.
- Region panel: within the PRI zone, IHD/CBD/GB remain visibly separated (as in Fig.1B); within the ORG zone, the three regions are thoroughly mixed (no visible separation) -- reproducing the paper's claim that organoid culture erases regional identity; the BTO zone shows partial region structure re-emerging, consistent with bile-induced partial reacquisition of regional identity.

## Still open
- This is a qualitative reproduction of the reported topology, not a quantitative match to Fig.2A's exact cell positions (different batch-correction method and requantified counts, per document/03_normalize.md and document/02_qc.md).

## Revision (style match to the published figure)

Redrawn as a **single panel** colored by the full 9-way origin x region
category (PRI/ORG/BTO x IHD/CBD/GB), replacing the earlier two-panel
(origin-only, region-only) layout, to match Fig.2A's actual design. Colors,
black point outlines, boxed axes (all 4 spines), plain "UMAP 1"/"UMAP 2"
axis-label text, and the 3x3 legend grid (columns = region, rows = origin)
are approximated **by eye** from the published PDF -- the authors did not
release a machine-readable palette, so these hex values are a visual match,
not a verified one. See `scripts/06_fig2a.py` for the exact color table.
