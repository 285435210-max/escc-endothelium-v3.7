# -*- coding: utf-8 -*-
"""Nested-CV sensitivity: re-run DEG + tip intersection + z-scoring INSIDE each inner
training split, then compare out-of-fold AUC with the original (inner folds only chose C).

Same data, seeds and outer folds as 32_panel_nested_cv.py. No plotting, no external cohorts.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
DATA = r"D:\龙虾文档\2026-10-21-20-33-54\data"
OUT = r"D:\escc_sc\out_v43"
os.makedirs(OUT, exist_ok=True)
SEED = 20261105
C_GRID = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
FC_CUT, FDR_CUT = 0.585, 0.05

expr = pd.read_csv(os.path.join(DATA, "GSE53625", "expr_gene.csv"), index_col=0)
expr.columns = [c.strip() for c in expr.columns]
expr = expr.groupby(expr.index).mean()
grp = pd.read_csv(os.path.join(DATA, "GSE53625", "group.csv"))
grp["sample"] = grp["sample"].str.strip()
grp = grp[grp["sample"].isin(expr.columns)].reset_index(drop=True)
expr = expr[grp["sample"].tolist()]
clin = pd.read_csv(os.path.join(ROOT, "data_manifest", "GSE53625_clinical_all.csv"),
                   encoding="utf-8-sig")
clin.columns = [c.strip() for c in clin.columns]
pmap = dict(zip(clin["sample"].astype(str).str.strip(),
                clin["patient id"].astype(str).str.strip()))
patients = np.array([pmap.get(s, s) for s in grp["sample"]])
y = (grp["group"] == "Tumor").astype(int).to_numpy()
cols = grp["sample"].tolist()
tip_set = set(pd.read_csv(os.path.join(ROOT, "09_de", "Tip_signature_shared.csv"))
              ["gene"].astype(str).tolist())
print("[data] genes=%d samples=%d patients=%d tip=%d" %
      (expr.shape[0], expr.shape[1], len(set(patients)), len(tip_set)), flush=True)


def fit_transform_z(df, fit_cols, all_cols=None):
    if all_cols is None:
        all_cols = list(df.columns)
    mu = df[fit_cols].mean(axis=1)
    sd = df[fit_cols].std(axis=1, ddof=1).replace(0, np.nan).fillna(1.0)
    return df[all_cols].sub(mu, axis=0).div(sd, axis=0)


def train_deg(df, tr_cols, y_tr):
    tr_cols = np.asarray(tr_cols)
    T = df[list(tr_cols[y_tr == 1])].to_numpy(dtype=float)
    N = df[list(tr_cols[y_tr == 0])].to_numpy(dtype=float)
    with np.errstate(all="ignore"):
        _, p = stats.ttest_ind(T, N, axis=1, equal_var=False, nan_policy="omit")
        lfc = np.nanmean(T, axis=1) - np.nanmean(N, axis=1)
    padj = multipletests(np.nan_to_num(p, nan=1.0), method="fdr_bh")[1]
    ok = (padj < FDR_CUT) & (np.abs(lfc) >= FC_CUT)
    return set(g for g, k in zip(df.index, ok) if k)


def clf_for(C):
    return LogisticRegression(l1_ratio=1.0, solver="liblinear",
                              random_state=SEED, C=C, max_iter=2000)


def inner_fixed(Xtr, y_tr, p_tr, seed):
    best_c, best_auc = C_GRID[0], -1.0
    cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed)
    splits = list(cv.split(np.zeros(len(y_tr)), y_tr, groups=p_tr))
    for C in C_GRID:
        oof = np.zeros(len(y_tr))
        for i, j in splits:
            m = clf_for(C).fit(Xtr[i], y_tr[i])
            oof[j] = m.predict_proba(Xtr[j])[:, 1]
        a = roc_auc_score(y_tr, oof)
        if a > best_auc:
            best_auc, best_c = a, C
    return best_c, best_auc


def inner_reselect(tr_cols, y_tr, p_tr, seed, outer_cand):
    best_c, best_auc, cand_sizes = C_GRID[0], -1.0, []
    cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed)
    splits = list(cv.split(np.zeros(len(y_tr)), y_tr, groups=p_tr))
    for C in C_GRID:
        oof = np.zeros(len(y_tr))
        for i, j in splits:
            inner_fit = [tr_cols[t] for t in i]
            deg_i = train_deg(expr, inner_fit, y_tr[i])
            cand_i = [g for g in expr.index if g in deg_i and g in tip_set]
            if len(cand_i) < 2:
                cand_i = outer_cand
            cand_sizes.append(len(cand_i))
            z = fit_transform_z(expr.loc[cand_i], inner_fit, tr_cols)
            X = z.T.to_numpy(dtype=float)
            m = clf_for(C).fit(X[i], y_tr[i])
            oof[j] = m.predict_proba(X[j])[:, 1]
        a = roc_auc_score(y_tr, oof)
        if a > best_auc:
            best_auc, best_c = a, C
    return best_c, best_auc, float(np.mean(cand_sizes))


def bootstrap_ci(y, score, patients_arr, n_boot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    uniq = np.unique(patients_arr)
    idx = {p: np.where(patients_arr == p)[0] for p in uniq}
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        ix = np.concatenate([idx[p] for p in pick])
        if len(np.unique(y[ix])) < 2:
            continue
        vals.append(roc_auc_score(y[ix], score[ix]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


outer = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
folds = list(outer.split(np.zeros(len(y)), y, groups=patients))
oof = {"original": np.zeros(len(y)), "variant": np.zeros(len(y))}
rows = []
for k, (i, j) in enumerate(folds, 1):
    tr_cols = [cols[t] for t in i]
    ytr, yte = y[i], y[j]
    deg_tr = train_deg(expr, tr_cols, ytr)
    cand = [g for g in expr.index if g in deg_tr and g in tip_set]
    z = fit_transform_z(expr.loc[cand], tr_cols, cols)
    X = z.T.to_numpy(dtype=float)
    Xtr, Xte = X[i], X[j]

    C_o, A_o = inner_fixed(Xtr, ytr, patients[i], SEED + k)
    C_v, A_v, ncand = inner_reselect(tr_cols, ytr, patients[i], SEED + k, cand)
    aucs = {}
    for tag, C in [("original", C_o), ("variant", C_v)]:
        m = clf_for(C).fit(Xtr, ytr)
        oof[tag][j] = m.predict_proba(Xte)[:, 1]
        aucs[tag] = roc_auc_score(yte, oof[tag][j])
    rows.append(dict(fold=k, n_train=len(i), n_test=len(j),
                     candidates_outer=len(cand), mean_candidates_inner=round(ncand, 1),
                     C_original=C_o, C_variant=C_v,
                     inner_auc_original=round(A_o, 4), inner_auc_variant=round(A_v, 4),
                     test_auc_original=round(aucs["original"], 4),
                     test_auc_variant=round(aucs["variant"], 4)))
    print("[fold %d] cand=%d innerCand=%.1f C=%s/%s innerAUC=%.3f/%.3f testAUC=%.4f/%.4f" %
          (k, len(cand), ncand, C_o, C_v, A_o, A_v,
           aucs["original"], aucs["variant"]), flush=True)

summary = []
for tag in ("original", "variant"):
    a = roc_auc_score(y, oof[tag])
    lo, hi = bootstrap_ci(y, oof[tag], patients)
    summary.append(dict(method=tag, oof_auc=round(a, 4),
                        ci95_low=round(lo, 4), ci95_high=round(hi, 4),
                        mean_C=round(float(np.mean([r["C_" + tag] for r in rows])), 3),
                        mean_inner_auc=round(float(np.mean(
                            [r["inner_auc_" + tag] for r in rows])), 4)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "nestedCV_inner_reselect_per_fold.csv"),
                          index=False, encoding="utf-8-sig")
pd.DataFrame(summary).to_csv(os.path.join(OUT, "nestedCV_inner_reselect_summary.csv"),
                             index=False, encoding="utf-8-sig")
print("\n=== SUMMARY ===")
print(pd.DataFrame(summary).to_string(index=False))
d = summary[1]["oof_auc"] - summary[0]["oof_auc"]
print("\nOOF AUC difference (variant - original) = %+.4f" % d)
print("DONE")
