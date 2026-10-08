# 08_region_classifier — regional identity classifier for primary cholangiocytes

**Question** Can IHD / CBD / GB identity of primary cholangiocytes be predicted from expression on a NEW donor, and does a primary-trained classifier still recognise region in organoids (ORG) and bile-treated organoids (BTO)?

**Script** `scripts/08_region_classifier.py`

## Inputs
- outputs/objects/02_qc_filtered.h5ad (log-normalized in the script, target_sum 1e4)
- data/metadata/region_markers_ensembl.json (18 region markers from paper Fig.1D/1E, symbol -> Ensembl ID via Ensembl REST; CYR61 queried as its current symbol CCN1)

## What was done
- Features: 2000 HVGs selected on primary (PRI) cells (Seurat flavor), z-scored inside each training fold.
- Model: multinomial logistic regression (L2, C = 0.01, class_weight = balanced).
- Validation: leave-one-donor-out (LODO) over the 10 primary donors. Region is nested in donor (CBD = Donor 1-3, GB = 4-6, IHD = 7-10), so a random cell split would leak donor identity.
- Baseline: same model on the 18 marker genes only.
- Final model trained on all PRI cells, applied to ORG and BTO cells; mean predicted probability per origin x region x `biosd_sample`.

## Outputs
- outputs/figures/08_region_classifier.png / .pdf (A: LODO confusion matrix, B: marker mean expression z-scored across 9 groups, C: mean predicted probability for PRI (LODO) / ORG / BTO)
- outputs/tables/08_lodo_per_donor.csv, 08_lodo_confusion_hvg2000.csv, 08_top_genes_per_region_ensembl.csv (top 15 coefficients per region, Ensembl IDs), 08_culture_predicted_prob.csv

## Findings
- LODO balanced accuracy: 0.976 (2000 HVG) vs 0.937 (18 markers).
- LODO recall (HVG model): IHD 0.95, CBD 0.98, GB 1.00.
- Per donor (HVG model): Donor 1-7 and 10 reach 0.96-1.00. Donor 8 (12 cells) = 0.58 and Donor 9 (13 cells) = 0.62 -- too few cells to interpret. Marker model is weaker for Donor 3 (0.82), Donor 4 (0.93).
- ORG: PRI-trained model gives mixed CBD/GB probabilities for all three regions (ORG CBD p(CBD) 0.72-0.82; ORG IHD p(CBD) 0.52-0.67; ORG GB p(CBD) 0.42-0.66, p(GB) 0.30-0.54) -- regional identity is not evident, consistent with the paper's loss of regional differences in culture.
- BTO: every region, including BTO CBD and BTO IHD, is predicted GB (p(GB) 0.86-0.94) -- consistent with the paper's shift toward gallbladder identity after bile treatment.

## Still open
- Domain shift: the model never saw culture-related genes, so probabilities on ORG/BTO are relative, not calibrated. A model trained jointly with ORG would be needed to claim "loss" quantitatively.
- ORG/BTO have only 2 samples per region, so sample effect and region effect cannot be separated.
- HVGs are chosen once on all PRI cells (unsupervised, minor leakage into LODO folds).
- Top genes are Ensembl IDs only; they have not been converted to symbols or tested for pathway enrichment.
- The `individual` field for ORG/BTO samples in the source metadata is unreliable (e.g. a BTO CBD sample labelled "Donor 8", the label of a primary IHD donor); organoid analyses here use `biosd_sample`.
