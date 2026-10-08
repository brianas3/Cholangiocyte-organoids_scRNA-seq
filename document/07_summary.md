# 07_summary — UMAP reproduction summary vs Sampaziotis et al. 2021

**Question** How closely does this pipeline, run on the Atlas-hosted re-quantification of E-MTAB-8495, reproduce the paper's Fig.1B and Fig.2A UMAPs and the underlying cluster-fidelity statistics?

**Script** `scripts/07_summary.py`

## Inputs
- document/02_qc.md
- document/03_normalize.md
- document/04_cluster.md
- document/06_fig2a.md
- scripts/04b_primary_nocorrect.py, 04c_full_nocorrect.py, 05b_seed_robustness.py (2026-10-08 fix; see Findings)

## What was done
- Compared reproduced cell counts per origin/region (outputs/tables/02_qc_cells_per_sample.csv) against the paper's reported numbers.
- Compared reproduced Louvain-cluster ARI/AMI vs origin/region (outputs/tables/04_cluster_louvain_ari_ami.csv) against the paper's reported >0.95 (origin) / <0.3 (region).
- Visually compared outputs/figures/05_fig1b_primary_umap.png and outputs/figures/06_fig2a_full_umap.png against Fig.1B and Fig.2A of the paper.

## Outputs
- This document (document/07_summary.md) -- no new data outputs.

## Findings
- CELL COUNTS: organoid (ORG) and bile-treated organoid (BTO) counts land within ~1.06-1.21x of the paper across all 3 regions. Primary tissue (PRI) diverges more: IHD 357 vs paper 587 (0.61x), GB 2985 vs 3702 (0.81x), CBD 3974 vs 3006 (1.32x).
- CLUSTER FIDELITY: Louvain resolution 0.05 gives exactly 3 clusters (as in the paper). Mean ARI(origin)=0.65/AMI=0.59 vs ARI(region)=0.01/AMI=0.02 across 10 seeds -- the SAME qualitative conclusion as the paper (clusters track origin, not region) but weaker in absolute terms than the paper's >0.95/<0.3.
- BATCH-CORRECTION FIX (2026-10-08, supersedes the earlier UMAP findings): every donor/sample belongs to exactly ONE region (PRI: CBD = Donor 1-3, GB = 4-6, IHD = 7-10; all 22 samples map to a single origin x region category), so region is fully nested in the Harmony batch key. Harmony (03_normalize.py, 04_cluster.py) therefore erased region structure: kNN same-category purity at 30 PCs, Harmony vs uncorrected PCA = PRI CBD 0.75 vs 0.99, PRI GB 0.58 vs 0.99, PRI IHD 0.38 vs 0.97 (full dataset); primary-only IHD 0.52 vs 0.96. Fig.1B and Fig.2A are now rebuilt WITHOUT Harmony: scripts/04b_primary_nocorrect.py and scripts/04c_full_nocorrect.py shift each donor/sample centroid onto its region (or origin x region) centroid, never across regions (X_pca_wr / X_pca_wc). Same-category purity after correction: PRI 0.98-0.995, ORG 0.88-0.92, BTO CBD 0.94, BTO GB/IHD 0.51-0.53 (outputs/tables/04b_primary_knn_region_purity.csv, 04c_full_knn_category_purity.csv). Without correction, Donor 3 vs Donor 1+2 split PRI CBD into two clusters (donor effect, not biology); the within-region centering merges them.
- UMAP TOPOLOGY (Fig.1B, 7,205 primary cells after removing 2 small Leiden clusters of 111 cells; paper 7,295 / 10 individuals): IHD, CBD and GB each form one main cluster; CBD is a single cluster. Donor shades inside each region are intermixed (the paper's legend shows 4/3/3 shades = 10 individuals; the caption does not say shades are donors -- inferred). Small stray fragments remain in GB and CBD in every seed and were not characterised. Seed robustness (scripts/05b_seed_robustness.py, seeds 0-4): UMAP-embedding same-region purity CBD 0.995, GB 0.994, IHD 0.973, sd <= 0.003 -- seed only changes layout (rotation, relative position, shape), not the conclusion; the paper's tidy three-blob layout will not be reproduced exactly (Seurat vs scanpy UMAP, different input counts).
- UMAP TOPOLOGY (Fig.2A, 40,732 cells vs paper's 35,603): PRI IHD/CBD/GB are three distinct groups as in the paper; ORG IHD and CBD largely overlap with ORG GB adjacent (not fully mixed); BTO CBD forms its own cluster while BTO GB and BTO IHD overlap (paper: partial overlap). BTO CBD separation may reflect its only 2 samples (sample effect retained) -- hypothesis, not tested.
- NOT RE-RUN: the Louvain ARI/AMI above still uses the Harmony embedding (04_cluster.py), so it inherits the nesting problem; rerun on X_pca_wc before quoting it.

## Still open
- ROOT CAUSE OF REMAINING GAPS (in order of likely impact; Harmony/region nesting, now fixed for the figures, was the largest cause of the unseparated UMAPs): (1) The paper's own QC'd/batch-corrected object was never deposited -- only raw fastq are on ArrayExpress (E-MTAB-8495) -- so this pipeline necessarily starts from EBI Single Cell Expression Atlas's independent re-quantification of those fastq (different CellRanger version/reference/cell-calling), not the authors' original counts. (2) Batch correction: Harmony was replaced by within-region centering (a mean shift only, NOT fastMNN and not the paper's method) after Harmony proved unusable for this nested design; fastMNN was not used because bioconductor-batchelor has no osx-arm64 conda build and the pip fallback mnnpy fails to compile (Apple clang lacks -fopenmp) on this Mac -- see document/03_normalize.md and document/01... (env build notes in step-status log). (3) scater::isOutlier's exact log-transform/direction convention is not fully specified in the paper's Methods text; this pipeline follows the standard scran/OSCA convention. (4) The paper's IHD primary group of 587 cells across '5 patients' includes a 5th sample that is the external MacParland et al. 2018 dataset (cluster 17), which is NOT part of the E-MTAB-8495 accession and was therefore excluded here -- this alone likely explains most of the IHD undercount.
- Given (1)-(4), this is a faithful re-implementation of the DESCRIBED pipeline on the SAME underlying fastq, reproducing the paper's qualitative conclusions (regional diversity in primary tissue; convergence in organoid culture; partial reacquisition after bile treatment) -- it is not, and cannot be, a cell-for-cell reproduction of the original figures.
- Not attempted in this pass: PAGA connectivity, diffusion pseudotime/Monocle2, SC3+clustree subpopulation search, and the post-hoc small-cluster removal step -- all of which the paper also applies specifically to the primary-tissue analysis (fig. S3-S5) and could be added as further steps.

