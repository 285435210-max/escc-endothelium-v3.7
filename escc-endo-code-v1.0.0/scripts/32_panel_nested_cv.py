# -*- coding: utf-8 -*-
"""
B1 修订（审稿意见 M1）：诊断 panel 的三重数据泄漏修复。

旧流程 scripts/11_tip_panel.py 的问题：
  1) 候选基因 = tip 487 ∩ GSE53625 全队列 7,190 DEG  -> 特征选择用到了全部标签
  2) z-score 在折外用全队列均值/方差                  -> 预处理泄漏
  3) StratifiedKFold(5) 未按患者分组                  -> 同一患者 tumor/normal 跨折
  4) 全量重拟合后在同一队列算 AUC = 1.000             -> 重代入(resubstitution)，非验证
  5) 200 组随机基因对照无生成代码、抽样池未说明

本脚本：
  - 外层 StratifiedGroupKFold(5, groups=patient_id)，同一个患者的配对样本永不跨折
  - 每个训练折内重新做 DEG（Welch t + BH + |log2FC|>=0.585），再与 tip 签名取交集
  - z-score 的均值/方差只用训练折
  - 内层 StratifiedGroupKFold(3) 选 L1 惩罚参数 C（嵌套 CV）
  - 只报 out-of-fold AUC；GSE53625 明确标为训练/发现队列
  - 外部队列 GSE20347 / GSE45670 用全量训练模型做独立验证（合法：未参与特征选择）
  - 随机对照从两个池各抽 200 组：7,190 个全队列 DEG 池、221 个 tip∩DEG 候选池
  - 所有随机过程固定种子

输出：09_de/panel/revised/
"""
import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from escc_lib import roc_points

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
DATA = r"D:\龙虾文档\2026-10-21-20-33-54\data"
OUT = os.path.join(ROOT, "09_de", "panel", "revised")
os.makedirs(OUT, exist_ok=True)
SEED = 20261105
C_GRID = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
FC_CUT, FDR_CUT = 0.585, 0.05

# ------------------------------------------------------------------ 数据
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
pmap = dict(zip(clin["sample"].astype(str).str.strip(), clin["patient id"].astype(str).str.strip()))
patients = np.array([pmap.get(s, s) for s in grp["sample"]])
y = (grp["group"] == "Tumor").astype(int).to_numpy()
cols = grp["sample"].tolist()

tip = pd.read_csv(os.path.join(ROOT, "09_de", "Tip_signature_shared.csv"))["gene"].astype(str).tolist()
tip_set = set(tip)
old_coef = pd.read_csv(os.path.join(ROOT, "09_de", "panel", "panel_coefficients.csv"),
                       encoding="utf-8-sig", index_col=0)["coef"]
sel_genes = old_coef[old_coef != 0].index.tolist()
sel_beta = old_coef[old_coef != 0].to_numpy()
print(f"[data] genes={expr.shape[0]} samples={expr.shape[1]} patients={len(set(patients))} "
      f"tip={len(tip_set)} full-model={len(sel_genes)}", flush=True)


# ------------------------------------------------------------------ 工具
def fit_transform_z(df, fit_cols, all_cols=None):
    """逐基因按样本标准化；均值/方差只来自 fit_cols。"""
    if all_cols is None:
        all_cols = list(df.columns)
    mu = df[fit_cols].mean(axis=1)
    sd = df[fit_cols].std(axis=1, ddof=1).replace(0, np.nan).fillna(1.0)
    return df[all_cols].sub(mu, axis=0).div(sd, axis=0)


def train_deg(df, tr_cols, y_tr):
    """训练折内重算肿瘤-癌旁 DEG，返回显著基因集合。口径与 01_deg_gse53625.py 一致。"""
    tr_cols = np.asarray(tr_cols)
    T = df[list(tr_cols[y_tr == 1])].to_numpy(dtype=float)
    N = df[list(tr_cols[y_tr == 0])].to_numpy(dtype=float)
    with np.errstate(all="ignore"):
        _, p = stats.ttest_ind(T, N, axis=1, equal_var=False, nan_policy="omit")
        lfc = np.nanmean(T, axis=1) - np.nanmean(N, axis=1)
    padj = multipletests(np.nan_to_num(p, nan=1.0), method="fdr_bh")[1]
    ok = (padj < FDR_CUT) & (np.abs(lfc) >= FC_CUT)
    return set(g for g, k in zip(df.index, ok) if k)


def inner_choose_C(X, y_tr, p_tr, seed):
    """训练折内 3 折分组 CV 选 C，按内层 OOF AUC。"""
    best_c, best_auc = C_GRID[0], -1.0
    cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=seed)
    splits = list(cv.split(X, y_tr, groups=p_tr))
    for C in C_GRID:
        oof = np.zeros(len(y_tr))
        for i, j in splits:
            clf = LogisticRegression(l1_ratio=1.0, solver="liblinear", random_state=SEED, C=C, max_iter=2000)
            clf.fit(X[i], y_tr[i])
            oof[j] = clf.predict_proba(X[j])[:, 1]
        a = roc_auc_score(y_tr, oof)
        if a > best_auc:
            best_auc, best_c = a, C
    return best_c, best_auc


def bootstrap_auc_by_patient(y, score, patients, n_boot=2000, seed=SEED):
    """按患者整簇重采样（保持配对样本同进同出）。"""
    rng = np.random.default_rng(seed)
    uniq = np.unique(patients)
    idx_by_p = {p: np.where(patients == p)[0] for p in uniq}
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        ix = np.concatenate([idx_by_p[p] for p in pick])
        yy, ss = y[ix], score[ix]
        if len(np.unique(yy)) < 2:
            continue
        vals.append(roc_auc_score(yy, ss))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def cohort_auc(name, beta, genes):
    """外部队列：cohort-wise 标准化（无标签泄漏），用全量训练模型系数打分。"""
    e = pd.read_csv(os.path.join(DATA, "validation", f"{name}_expr_gene.csv"), index_col=0)
    e.columns = [c.strip() for c in e.columns]
    e = e.groupby(e.index).mean()
    g = pd.read_csv(os.path.join(DATA, "validation", f"{name}_group.csv"))
    g["sample"] = g["sample"].str.strip()
    g = g[g["sample"].isin(e.columns)].reset_index(drop=True)
    e = e[g["sample"].tolist()]
    use = [x for x in genes if x in e.index]
    z = fit_transform_z(e.loc[use], list(e.columns)).T.to_numpy()
    b = np.array([beta[genes.index(x)] for x in use])
    s = z @ b
    yy = (g["group"] == "Tumor").astype(int).to_numpy()
    a = roc_auc_score(yy, s)
    lo, hi = bootstrap_auc_by_patient(yy, s, np.arange(len(yy)), n_boot=2000)
    return dict(cohort=name, n_tumor=int(yy.sum()), n_normal=int((1 - yy).sum()),
                n_genes=len(use), AUC=round(a, 4), lower95=round(lo, 4), upper95=round(hi, 4)), (yy, s)


# ------------------------------------------------------------------ 外层划分
outer = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
folds = list(outer.split(np.zeros(len(y)), y, groups=patients))
for k, (i, j) in enumerate(folds, 1):
    assert len(set(patients[i]) & set(patients[j])) == 0, "患者跨折！"
print("[split] 5 折按患者分组，无患者跨折 ✓", flush=True)


# ------------------------------------------------------------------ 主分析：嵌套 CV
oof = np.zeros(len(y), dtype=float)
fold_rows, fold_C = [], []
for k, (i, j) in enumerate(folds, 1):
    tr_cols = [cols[t] for t in i]
    te_cols = [cols[t] for t in j]
    ytr, yte = y[i], y[j]

    deg_tr = train_deg(expr, tr_cols, ytr)
    cand = [g for g in expr.index if g in deg_tr and g in tip_set]
    z = fit_transform_z(expr.loc[cand], tr_cols, cols)
    X = z.T.to_numpy(dtype=float)
    Xtr, Xte = X[i], X[j]

    C, inner_auc = inner_choose_C(Xtr, ytr, patients[i], seed=SEED + k)
    fold_C.append(C)
    clf = LogisticRegression(l1_ratio=1.0, solver="liblinear", random_state=SEED, C=C, max_iter=2000).fit(Xtr, ytr)
    oof[j] = clf.predict_proba(Xte)[:, 1]
    nz = [g for g, c in zip(cand, clf.coef_[0]) if c != 0]
    fold_rows.append(dict(fold=k, n_train=len(i), n_test=len(j),
                          deg_in_fold=len(deg_tr), candidates=len(cand),
                          C=C, inner_cv_auc=round(inner_auc, 4),
                          n_selected=len(nz), test_auc=round(roc_auc_score(yte, oof[j]), 4),
                          selected=";".join(nz)))
    print(f"[fold {k}] DEG={len(deg_tr)} cand={len(cand)} C={C} "
          f"innerAUC={inner_auc:.3f} selected={len(nz)} testAUC={fold_rows[-1]['test_auc']}", flush=True)

auc_oof = roc_auc_score(y, oof)
lo, hi = bootstrap_auc_by_patient(y, oof, patients)
print(f"\n[MAIN] nested/grouped CV out-of-fold AUC = {auc_oof:.4f} (95% CI {lo:.4f}-{hi:.4f})", flush=True)
# 折外预测供 Figure 7 的 ROC 曲线使用（避免把重代入曲线画进图里）
oof_df = pd.DataFrame({"sample": cols, "patient": patients, "y": y, "oof_score": oof})
oof_df.to_csv(os.path.join(OUT, "oof_predictions.csv"), index=False, encoding="utf-8-sig")
oof_df.to_csv(r"D:\escc_sc\data\panel_oof_scores.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(fold_rows).to_csv(os.path.join(OUT, "nested_cv_per_fold.csv"), index=False, encoding="utf-8-sig")

# 基因稳定性
from collections import Counter
cnt = Counter(g for r in fold_rows for g in r["selected"].split(";") if g)
stab = pd.DataFrame({"gene": list(cnt), "folds_selected": list(cnt.values())}).sort_values(
    "folds_selected", ascending=False)
stab.to_csv(os.path.join(OUT, "gene_selection_stability.csv"), index=False, encoding="utf-8-sig")


# ------------------------------------------------------------------ 外部队列（全量训练模型）
ext_rows, curves = [], {"GSE53625 (nested CV OOF)": (y, oof)}
for nm in ["GSE20347", "GSE45670"]:
    r, (yy, ss) = cohort_auc(nm, sel_beta, sel_genes)
    ext_rows.append(r)
    curves[nm] = (yy, ss)
    print(f"[external] {nm}: AUC={r['AUC']} ({r['lower95']}-{r['upper95']})", flush=True)
pd.DataFrame(ext_rows).to_csv(os.path.join(OUT, "external_validation_revised.csv"),
                              index=False, encoding="utf-8-sig")


# ------------------------------------------------------------------ 随机基因对照（四个抽样池）
# 旧稿写 "200 random 23-gene sets gave median AUC 0.531/0.554/0.579"，但全仓没有生成代码，
# 抽样池也未说明。这里显式给出四个池，其中 foldDEG 与主分析条件严格一致。
full_deg = pd.read_csv(os.path.join(ROOT, "01_deg", "GSE53625_DEG_sig.csv"), encoding="utf-8-sig")
pool_a = [g for g in full_deg["gene"].astype(str) if g in expr.index]          # 全队列 7,190 DEG
pool_b = [g for g in old_coef.index.astype(str) if g in expr.index]           # 221 tip∩DEG 候选
pool_g = list(expr.index)                                                      # 全基因组（旧稿口径推测）
fold_deg_sets = [sorted(train_deg(expr, [cols[t] for t in i], y[i])) for i, _ in folds]
POOL_LABEL = {
    "DEG7190": "all-cohort DEGs (feature-selection leakage, permissive control)",
    "candidate221": "221 tip/DEG candidate genes",
    "foldDEG": "fold-internal DEGs (strictly matched to main analysis)",
    "genome": "all measured genes (pool used by the earlier unreported control)",
}
pools = {"DEG7190": pool_a, "candidate221": pool_b, "foldDEG": None, "genome": pool_g}
rng = np.random.default_rng(SEED)
n_rand, k_size = 200, len(sel_genes)
null_rows = []
for pname, pool in pools.items():
    for b in range(n_rand):
        fixed = None if pname == "foldDEG" else [pool[t] for t in
                                                 rng.choice(len(pool), size=k_size, replace=False)]
        po = np.zeros(len(y))
        for k, (i, j) in enumerate(folds):
            gset = ([fold_deg_sets[k][t] for t in
                     rng.choice(len(fold_deg_sets[k]), size=k_size, replace=False)]
                    if pname == "foldDEG" else fixed)
            tr_cols = [cols[t] for t in i]
            zz = fit_transform_z(expr.loc[gset], tr_cols, cols)
            X = zz.T.to_numpy(dtype=float)
            clf = LogisticRegression(l1_ratio=1.0, solver="liblinear",
                                     C=fold_C[k], max_iter=2000).fit(X[i], y[i])
            po[j] = clf.predict_proba(X[j])[:, 1]
        null_rows.append(dict(pool=pname, rep=b + 1, auc=roc_auc_score(y, po)))
    sub = [r["auc"] for r in null_rows if r["pool"] == pname]
    emp_p = (1 + sum(a >= auc_oof for a in sub)) / (n_rand + 1)
    print(f"[null] {pname:<13} median={np.median(sub):.4f} "
          f"range=[{min(sub):.4f},{max(sub):.4f}] empirical P={emp_p:.4f}", flush=True)

null_df = pd.DataFrame(null_rows)
null_df.to_csv(os.path.join(OUT, "random_gene_control_auc.csv"), index=False, encoding="utf-8-sig")
null_summary = []
for pname in pools:
    sub = null_df.loc[null_df.pool == pname, "auc"].to_numpy()
    null_summary.append(dict(pool=pname, description=POOL_LABEL[pname],
                             n_sets=len(sub), size=k_size,
                             median=round(float(np.median(sub)), 4),
                             q025=round(float(np.percentile(sub, 2.5)), 4),
                             q975=round(float(np.percentile(sub, 97.5)), 4),
                             min=round(float(sub.min()), 4), max=round(float(sub.max()), 4),
                             empirical_P_vs_nestedCV=round(float((1 + (sub >= auc_oof).sum()) / (n_rand + 1)), 4)))
pd.DataFrame(null_summary).to_csv(os.path.join(OUT, "random_gene_control_summary.csv"),
                                  index=False, encoding="utf-8-sig")


# ------------------------------------------------------------------ 负对照：患者内标签置换
# 在每个患者内部随机交换 tumor/normal 标签（保留配对结构、组间平衡与全部预处理），
# 用固定 23 基因重新跑同样的分组 CV。若流程无泄漏，AUC 应接近 0.5。
base_feat = [g for g in sel_genes if g in expr.index]
rng_p = np.random.default_rng(SEED + 7)
perm_aucs = []
for b in range(20):
    yp = y.copy()
    for p in np.unique(patients):
        ix = np.where(patients == p)[0]
        if len(ix) == 2 and rng_p.random() < 0.5:
            yp[ix] = yp[ix][::-1]
    po = np.zeros(len(y))
    for k, (i, j) in enumerate(folds):
        tr_cols = [cols[t] for t in i]
        zz = fit_transform_z(expr.loc[base_feat], tr_cols, cols)
        X = zz.T.to_numpy(dtype=float)
        clf = LogisticRegression(l1_ratio=1.0, solver="liblinear",
                                 C=fold_C[k], max_iter=2000).fit(X[i], yp[i])
        po[j] = clf.predict_proba(X[j])[:, 1]
    perm_aucs.append(roc_auc_score(yp, po))
pd.DataFrame({"permutation": range(1, 21), "auc": perm_aucs}).to_csv(
    os.path.join(OUT, "within_patient_label_permutation.csv"), index=False, encoding="utf-8-sig")
print(f"[perm] within-patient label permutation (n=20): "
      f"median={np.median(perm_aucs):.4f} range=[{min(perm_aucs):.4f},{max(perm_aucs):.4f}]", flush=True)


# ------------------------------------------------------------------ 图
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
ax = axes[0]
for nm, (yy, ss) in curves.items():
    fpr, tpr = roc_points(yy, ss)
    ax.plot(fpr, tpr, lw=1.7, label=f"{nm}  AUC={roc_auc_score(yy, ss):.3f}")
ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=0.8)
ax.set_xlabel("1 - Specificity"); ax.set_ylabel("Sensitivity")
ax.set_title("Tip-derived panel, leakage-free evaluation", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="lower right")
ax = axes[1]
for pname, c in [("DEG7190", "#4a9bc4"), ("candidate221", "#b03a2e"),
                 ("foldDEG", "#e0a32e"), ("genome", "#7b4fa0")]:
    sub = null_df.loc[null_df.pool == pname, "auc"]
    ax.hist(sub, bins=28, alpha=0.55, label=pname, color=c, edgecolor="white", linewidth=0.4)
ax.axvline(auc_oof, color="black", lw=1.4, ls="--", label=f"nested CV = {auc_oof:.3f}")
ax.set_xlabel("Out-of-fold AUC"); ax.set_ylabel("Random gene sets")
ax.set_title("Matched random-gene controls (200 sets each)", fontsize=9)
ax.legend(frameon=False, fontsize=7)
for a in axes:
    a.tick_params(labelsize=7)
    a.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "ROC_panel_revised.png"), dpi=600, bbox_inches="tight")
plt.close()

import sklearn
summary = {
    "design": "nested 5-fold StratifiedGroupKFold grouped by patient; fold-internal DEG + tip intersection; train-fold-only z-score; inner 3-fold grouped CV for C",
    "sklearn_version": sklearn.__version__,
    "old_manuscript_control": "200 random 23-gene sets, median AUC 0.531/0.554/0.579; no generating script or result file exists in the workspace, and the sampling pool was never stated. Withdrawn and replaced by four explicit pools.",
    "n_samples": int(len(y)), "n_patients": int(len(set(patients))),
    "gse53625_nested_cv_oof_auc": round(float(auc_oof), 4),
    "gse53625_ci95_patient_bootstrap": [round(lo, 4), round(hi, 4)],
    "old_resubstitution_auc": 1.000,
    "external": ext_rows,
    "random_control": null_summary,
    "within_patient_label_permutation": {
        "n": 20, "median": round(float(np.median(perm_aucs)), 4),
        "min": round(float(min(perm_aucs)), 4), "max": round(float(max(perm_aucs)), 4)},
    "interpretation": "Out-of-fold AUC stayed near 1 after all leakage was removed, but matched random 23-gene sets also reached ~0.999. Tumor-versus-adjacent discrimination is therefore near-trivial in these array cohorts and the tip-derived panel carries no demonstrable specificity; the earlier claim that random sets gave AUC ~0.53-0.58 is not reproducible and has been withdrawn.",
    "fold_C": fold_C,
    "genes_in_full_descriptive_model": len(sel_genes),
    "genes_selected_in_all_5_folds": stab.loc[stab.folds_selected == 5, "gene"].tolist(),
    "seed": SEED,
    "note": "Table 3 (23 genes) is the full-cohort descriptive model fitted on the discovery cohort; its discrimination in GSE53625 is estimated by nested CV only, not by resubstitution."
}
with open(os.path.join(OUT, "nested_cv_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)
print("\n[done] B1 ->", OUT, flush=True)
