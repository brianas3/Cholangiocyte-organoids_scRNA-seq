# Cholangiocyte organoids scRNA-seq -- reanalysis

Reanalysis of the single-cell RNA-seq data underlying:

> Sampaziotis F, et al. "Cholangiocyte organoids can repair bile ducts after
> transplantation in the human liver." *Science* 371:839-846 (2021).
> https://doi.org/10.1126/science.aaz6964

Goal: reproduce, as closely as the deposited data allow, the paper's
Fig. 1B (region-specific diversity of primary cholangiocytes) and Fig. 2A
(convergence of organoids across regions, partial re-divergence after
gallbladder bile treatment) UMAPs, plus the Louvain-cluster ARI/AMI fidelity
check reported in the Methods.

## Data source

Raw fastq were deposited by the authors on ArrayExpress/BioStudies as
**E-MTAB-8495**, but the authors' own processed/QC'd object (the exact cells
and coordinates behind the published figures) was never deposited.

This repo works from **EBI Single Cell Expression Atlas's own re-quantification**
of those fastq (Alevin-fry against an Ensembl reference + emptyDrops cell
calling), NOT the authors' original CellRanger output -- see
`document/07_summary.md` for what that implies for the numbers below.
Download the raw quantification matrix yourself from:
https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-8495/downloads
(`quantification-raw` file type) and place it under
`data/Result files/E-MTAB-8495-quantification-raw-files/` to rerun the
pipeline -- the matrix itself (>2 GB) is not committed to this repo.

## Pipeline (`scripts/`, run in order)

1. `02_qc.py` -- gene/cell QC reproducing the paper's Methods (3-MAD outlier
   filtering per `scater::isOutlier` convention, EPCAM/KRT7/KRT19 marker-positive
   filter, per-batch gene detection filter).
2. `03_normalize.py` -- Seurat-style LogNormalize, cell-cycle scoring and
   regression of the G2M-S score difference, 2000 HVGs, PCA, Harmony batch
   correction (substituting for `scran::fastMNN` -- see Deviations below).
3. `04_cluster.py` -- Louvain clustering (resolution scanned to exactly 3
   clusters) + UMAP, run once on the primary-tissue-only subset (Fig. 1B scope)
   and once on the full dataset (Fig. 2A scope); ARI/AMI of clusters vs.
   origin/region across 10 seeds.
4. `05_fig1b.py` / `06_fig2a.py` -- render the two UMAP figures.

## Key results (see `document/07_summary.md` for full detail)

| | This reproduction | Paper |
|---|---|---|
| Primary cells / individuals (Fig. 1B) | 7,316 / 10 | 7,295 / 10 |
| Full dataset cells (Fig. 2A) | 40,732 | 35,603 |
| Louvain ARI vs. origin (mean, 10 seeds) | 0.65 | >0.95 |
| Louvain ARI vs. region (mean, 10 seeds) | 0.01 | <0.30 |

Qualitative conclusions reproduce: primary cholangiocytes retain
region-specific identity (Fig. 1B analog), organoid culture erases that
identity while primary tissue keeps it, and bile treatment partially
restores it (Fig. 2A analog). Absolute cluster-fidelity numbers are weaker
than the paper's, attributable to the deviations below.

## Deviations from the paper's Methods

- **Batch correction**: paper used `scran::fastMNN`; this repo uses **Harmony**
  (`harmonypy`), because `bioconductor-batchelor` has no osx-arm64 conda build
  and the pip fallback `mnnpy` fails to compile (Apple clang lacks `-fopenmp`).
- **Underlying counts**: EBI Atlas re-quantification (Alevin-fry + Ensembl
  reference + emptyDrops), not the authors' original CellRanger run (never
  deposited).
- **IHD primary group**: the paper's reported 5-patient/587-cell IHD primary
  group includes a 5th sample from an external dataset (MacParland et al. 2018,
  cluster 17) that is not part of the E-MTAB-8495 accession and is therefore
  absent here.
- Not yet reproduced: PAGA connectivity, diffusion pseudotime / Monocle2,
  SC3+clustree subpopulation search, and the paper's post-hoc removal of
  small mito/gene-count-outlier Louvain clusters.

## Repo layout

```
data/       small sample-metadata tables (manifest, marker/cell-cycle gene IDs)
document/   per-step method notes + final comparison summary
outputs/    figures (PNG+PDF, 600dpi) and small result tables
scripts/    the 5 pipeline scripts, in run order
```

Large intermediate objects (raw count matrix, `.h5ad` checkpoints, ~1-4 GB
each) are intentionally excluded from version control (see `.gitignore`).
