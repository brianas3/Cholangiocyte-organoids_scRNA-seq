"""
Step 08 -- Cholangiocyte REGION-identity classifier (IHD / CBD / GB).

1. Train a multinomial logistic regression on primary (PRI) cells only, features =
   2000 HVGs (log-normalized, z-scored inside each training fold).
2. Validate with leave-one-DONOR-out CV (LODO). Donors are nested in region, so a
   random cell split would leak donor identity; LODO asks "does a region signature
   learned on other donors recognise a NEW donor of that region?". Reported per donor.
3. Compare with a 18-gene marker-only classifier (paper Fig.1D/1E markers, Ensembl IDs
   in data/metadata/region_markers_ensembl.json).
4. Apply the PRI-trained classifier to ORG and BTO cells: mean predicted probability per
   origin x region x sample (does regional identity vanish in organoids / return after bile?).

Caveats (also in document/08_region_classifier.md):
- HVGs are selected once on all PRI cells (unsupervised, minor leakage into LODO folds).
- IHD has 4 donors but Donor 8/9 have 12/13 cells -> their LODO accuracy is very noisy.
- ORG/BTO are a domain shift (culture genes); predictions there are relative, not calibrated.
- With 1-2 donors per organoid line, 'sample' effects cannot be separated from region.

Input:  outputs/objects/02_qc_filtered.h5ad, data/metadata/region_markers_ensembl.json
Output: outputs/figures/08_region_classifier.png / .pdf
        outputs/tables/08_*.csv
"""
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = "/Users/brian/Desktop/RDM/scRNA-seq"
REGIONS = ["IHD", "CBD", "GB"]
COLORS = {"IHD": "#D4A017", "CBD": "#7C8C3C", "GB": "#A13D2B"}


def make_clf():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(C=0.01, max_iter=2000, class_weight="balanced"))


def lodo(X, y, groups):
    pred = np.empty(len(y), dtype=object)
    proba = np.zeros((len(y), len(REGIONS)))
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        clf = make_clf().fit(X[tr], y[tr])
        pred[te] = clf.predict(X[te])
        # columns follow clf.classes_ (sorted); reorder to REGIONS
        idx = [list(clf.classes_).index(r) for r in REGIONS]
        proba[te] = clf.predict_proba(X[te])[:, idx]
    return pred, proba


def main():
    a = sc.read_h5ad(f"{ROOT}/outputs/objects/02_qc_filtered.h5ad")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    pri = a[a.obs["origin"] == "PRI"].copy()
    sc.pp.highly_variable_genes(pri, n_top_genes=2000, flavor="seurat")
    hvg = pri.var_names[pri.var["highly_variable"]]

    def dense(ad, genes):
        return np.asarray(ad[:, genes].X.todense()) if hasattr(ad.X, "todense") else np.asarray(ad[:, genes].X)

    y = pri.obs["region"].astype(str).values
    donors = pri.obs["individual"].astype(str).values
    Xh = dense(pri, hvg)

    with open(f"{ROOT}/data/metadata/region_markers_ensembl.json") as f:
        mk = json.load(f)
    mk_ids = [g for r in REGIONS for g in mk[r].values() if g in a.var_names]
    Xm = dense(pri, mk_ids)

    res = {}
    for name, X in [("hvg2000", Xh), ("markers", Xm)]:
        pred, proba = lodo(X, y, donors)
        per = (pd.DataFrame({"donor": donors, "region": y, "correct": pred == y})
               .groupby(["donor", "region"]).agg(n=("correct", "size"), acc=("correct", "mean")))
        per["model"] = name
        res[name] = (pred, proba, per.reset_index())
        print(name, "balanced acc (LODO):", round(balanced_accuracy_score(y, pred), 3))
    pd.concat([r[2] for r in res.values()]).to_csv(f"{ROOT}/outputs/tables/08_lodo_per_donor.csv", index=False)

    cm = confusion_matrix(y, res["hvg2000"][0], labels=REGIONS)
    cmn = cm / cm.sum(1, keepdims=True)
    pd.DataFrame(cmn, index=REGIONS, columns=REGIONS).round(3).to_csv(
        f"{ROOT}/outputs/tables/08_lodo_confusion_hvg2000.csv")

    # final model on all PRI; top genes (Ensembl IDs) per region
    final = make_clf().fit(Xh, y)
    coef = pd.DataFrame(final[-1].coef_.T, index=hvg, columns=final[-1].classes_)
    top = {r: coef[r].sort_values(ascending=False).head(15).index.tolist() for r in REGIONS}
    pd.DataFrame(top).to_csv(f"{ROOT}/outputs/tables/08_top_genes_per_region_ensembl.csv", index=False)

    # apply to ORG / BTO
    cult = a[a.obs["origin"].isin(["ORG", "BTO"])]
    P = final.predict_proba(dense(cult, hvg))
    P = pd.DataFrame(P[:, [list(final[-1].classes_).index(r) for r in REGIONS]], columns=REGIONS,
                     index=cult.obs_names)
    P["origin"], P["region"] = cult.obs["origin"].values, cult.obs["region"].astype(str).values
    P["sample"] = cult.obs["biosd_sample"].astype(str).values
    summ = P.groupby(["origin", "region", "sample"], observed=True)[REGIONS].mean().round(3)
    summ.to_csv(f"{ROOT}/outputs/tables/08_culture_predicted_prob.csv")
    print(summ)
    pri_P = pd.DataFrame(res["hvg2000"][1], columns=REGIONS)
    pri_P["region"] = y
    pri_summ = pri_P.groupby("region")[REGIONS].mean().round(3)

    # ---- figure ----
    cats = [f"{o} {r}" for o in ["PRI", "ORG", "BTO"] for r in REGIONS]
    sym = {v: k for r in REGIONS for k, v in mk[r].items()}
    Mx = []
    for o in ["PRI", "ORG", "BTO"]:
        for r in REGIONS:
            m = ((a.obs["origin"] == o) & (a.obs["region"] == r)).values
            Mx.append(np.asarray(a[m][:, mk_ids].X.mean(0)).ravel())
    Mx = pd.DataFrame(np.array(Mx), index=cats, columns=[sym[g] for g in mk_ids])
    Mz = (Mx - Mx.mean()) / Mx.std().replace(0, 1)

    fig = plt.figure(figsize=(12, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.9, 1.5], wspace=0.55)
    ax = fig.add_subplot(gs[0])
    ax.imshow(cmn, cmap="Blues", vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{cmn[i, j]:.2f}", ha="center", va="center",
                    color="white" if cmn[i, j] > 0.5 else "black", fontsize=8)
    ax.set_xticks(range(3)); ax.set_xticklabels(REGIONS, fontsize=8)
    ax.set_yticks(range(3)); ax.set_yticklabels(REGIONS, fontsize=8)
    ax.set_xlabel("Predicted", fontsize=8); ax.set_ylabel("True (PRI)", fontsize=8)
    ax.set_title("A  Leave-one-donor-out", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[1])
    im = ax.imshow(Mz.values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(Mz.shape[1])); ax.set_xticklabels(Mz.columns, rotation=60, ha="right",
                                                          fontsize=7, style="italic")
    ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7)
    ax.set_title("B  Marker mean expression (z across groups)", fontsize=9, loc="left")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)

    ax = fig.add_subplot(gs[2])
    groups = [("PRI", pri_summ.reindex(REGIONS))]
    for o in ["ORG", "BTO"]:
        groups.append((o, P[P.origin == o].groupby("region")[REGIONS].mean().reindex(REGIONS)))
    xs, labels, k = [], [], 0
    for o, df in groups:
        for r in REGIONS:
            bottom = 0
            for pr in REGIONS:
                ax.bar(k, df.loc[r, pr], bottom=bottom, color=COLORS[pr], edgecolor="black",
                       linewidth=0.3, width=0.8)
                bottom += df.loc[r, pr]
            labels.append(f"{o}\n{r}"); xs.append(k); k += 1
        k += 0.6
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Mean predicted probability", fontsize=8)
    ax.set_title("C  PRI-trained classifier on PRI (LODO) / ORG / BTO", fontsize=9, loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=COLORS[r], ec="black", lw=0.3) for r in REGIONS]
    ax.legend(handles, [f"p({r})" for r in REGIONS], fontsize=6, frameon=False,
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
    fig.savefig(f"{ROOT}/outputs/figures/08_region_classifier.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{ROOT}/outputs/figures/08_region_classifier.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
