# 07_summary — UMAP reproduction summary vs Sampaziotis et al. 2021

**Question** How closely does this pipeline, run on the Atlas-hosted re-quantification of E-MTAB-8495, reproduce the paper's Fig.1B and Fig.2A UMAPs and the underlying cluster-fidelity statistics?

**Script** `scripts/07_summary.py`

## Inputs
- document/02_qc.md
- document/03_normalize.md
- document/04_cluster.md
- document/06_fig2a.md

## What was done
- Compared reproduced cell counts per origin/region (outputs/tables/02_qc_cells_per_sample.csv) against the paper's reported numbers.
- Compared reproduced Louvain-cluster ARI/AMI vs origin/region (outputs/tables/04_cluster_louvain_ari_ami.csv) against the paper's reported >0.95 (origin) / <0.3 (region).
- Visually compared outputs/figures/05_fig1b_primary_umap.png and outputs/figures/06_fig2a_full_umap.png against Fig.1B and Fig.2A of the paper.

## Outputs
- This document (document/07_summary.md) -- no new data outputs.

## Findings
- CELL COUNTS: organoid (ORG) and bile-treated organoid (BTO) counts land within ~1.06-1.21x of the paper across all 3 regions. Primary tissue (PRI) diverges more: IHD 357 vs paper 587 (0.61x), GB 2985 vs 3702 (0.81x), CBD 3974 vs 3006 (1.32x).
- CLUSTER FIDELITY: Louvain resolution 0.05 gives exactly 3 clusters (as in the paper). Mean ARI(origin)=0.65/AMI=0.59 vs ARI(region)=0.01/AMI=0.02 across 10 seeds -- the SAME qualitative conclusion as the paper (clusters track origin, not region) but weaker in absolute terms than the paper's >0.95/<0.3.
- UMAP TOPOLOGY: Fig.1B analog (7,316 primary cells, 10 individuals, vs paper's 7,295/10) shows GB as a clearly separated island with CBD/IHD forming an adjacent, overlapping group -- matching the paper's finding that adjacent regions are more transcriptionally similar than distal ones.
- Fig.2A analog (40,732 cells vs paper's 35,603) shows PRI, ORG and BTO occupying three largely distinct zones, with region separation clearly visible within PRI, erased within ORG, and partially restored within BTO -- reproducing the paper's central plasticity claim.

## Still open
- ROOT CAUSE OF REMAINING GAPS (in order of likely impact): (1) The paper's own QC'd/batch-corrected object was never deposited -- only raw fastq are on ArrayExpress (E-MTAB-8495) -- so this pipeline necessarily starts from EBI Single Cell Expression Atlas's independent re-quantification of those fastq (different CellRanger version/reference/cell-calling), not the authors' original counts. (2) Batch correction used Harmony instead of scran::fastMNN, because bioconductor-batchelor has no osx-arm64 conda build and the pip fallback mnnpy fails to compile (Apple clang lacks -fopenmp) on this Mac -- see document/03_normalize.md and document/01... (env build notes in step-status log). (3) scater::isOutlier's exact log-transform/direction convention is not fully specified in the paper's Methods text; this pipeline follows the standard scran/OSCA convention. (4) The paper's IHD primary group of 587 cells across '5 patients' includes a 5th sample that is the external MacParland et al. 2018 dataset (cluster 17), which is NOT part of the E-MTAB-8495 accession and was therefore excluded here -- this alone likely explains most of the IHD undercount.
- Given (1)-(4), this is a faithful re-implementation of the DESCRIBED pipeline on the SAME underlying fastq, reproducing the paper's qualitative conclusions (regional diversity in primary tissue; convergence in organoid culture; partial reacquisition after bile treatment) -- it is not, and cannot be, a cell-for-cell reproduction of the original figures.
- Not attempted in this pass: PAGA connectivity, diffusion pseudotime/Monocle2, SC3+clustree subpopulation search, and the post-hoc small-cluster removal step -- all of which the paper also applies specifically to the primary-tissue analysis (fig. S3-S5) and could be added as further steps.

