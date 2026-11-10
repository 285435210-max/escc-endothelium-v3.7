# -*- coding: utf-8 -*-
"""
补充两项此前「悬空」的统计证据：
  1) 23 基因评分 Cox 模型的比例风险假设（scaled Schoenfeld 残差），Methods L78 原句无代码/输出
  2) 统一闭合成分口径下，6 个状态与病理分期/T/N 的 Spearman + BH（Figure 6 数字）
"""
import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
DATA = r"D:\龙虾文档\2026-10-21-20-33-54\data"
OUT = os.path.join(ROOT, "08_scRNA", "endothelial", "revised")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- 1. PH 假设
expr = pd.read_csv(os.path.join(DATA, "GSE53625", "expr_gene.csv"), index_col=0)
expr.columns = [c.strip() for c in expr.columns]
expr = expr.groupby(expr.index).mean()
grp = pd.read_csv(os.path.join(DATA, "GSE53625", "group.csv"))
grp["sample"] = grp["sample"].str.strip()
clin = pd.read_csv(os.path.join(ROOT, "data_manifest", "GSE53625_clinical_all.csv"),
                   encoding="utf-8-sig")
clin.columns = [c.strip() for c in clin.columns]
coef = pd.read_csv(os.path.join(ROOT, "09_de", "panel", "panel_coefficients.csv"),
                   encoding="utf-8-sig", index_col=0)["coef"]
sel = coef[coef != 0]
use = [g for g in sel.index if g in expr.index]
z = ((expr.loc[use].T - expr.loc[use].T.mean()) / expr.loc[use].T.std(ddof=1)).T
score = pd.Series(z.T.to_numpy() @ sel[use].to_numpy(), index=expr.columns, name="panel_score")

tc = grp.loc[grp.group == "Tumor", "sample"].tolist()
cl = clin[clin["sample"].astype(str).str.strip().isin(tc)].copy()
cl["sample"] = cl["sample"].astype(str).str.strip()
cl = cl.set_index("sample")
cl = cl.loc[[s for s in tc if s in cl.index]]
cl["T"] = pd.to_numeric(cl["OS.time"], errors="coerce")
cl["E"] = pd.to_numeric(cl["OS"], errors="coerce")
cl["score"] = score.loc[cl.index]
cl = cl.dropna(subset=["T", "E", "score"])
cl = cl[cl["T"] > 0]
cl["score_z"] = (cl["score"] - cl["score"].mean()) / cl["score"].std(ddof=1)

df = cl[["T", "E", "score_z"]].copy()
cph = CoxPHFitter().fit(df, duration_col="T", event_col="E")
hr = float(np.exp(cph.params_["score_z"]))
lo, hi = np.exp(cph.confidence_intervals_.loc["score_z"]).tolist()
p = float(cph.summary.loc["score_z", "p"])
ph = proportional_hazard_test(cph, df, time_transform="rank")
ph_p = float(ph.summary.loc["score_z", "p"])
print(f"[PH] n={len(df)} events={int(cl.E.sum())} HR/SD={hr:.3f} "
      f"({lo:.3f}-{hi:.3f}) P={p:.3f}; Schoenfeld P={ph_p:.4f}", flush=True)
pd.DataFrame({"model": ["panel_score_z"], "n": [len(df)], "events": [int(cl.E.sum())],
              "HR": [hr], "lower95": [lo], "upper95": [hi], "p": [p],
              "schoenfeld_p": [ph_p]}).to_csv(
    os.path.join(OUT, "cox_ph_test.csv"), index=False)

# ---------------------------------------------------------------- 2. 分期关联（闭合口径）
closed = pd.read_csv(os.path.join(OUT, "celltype_proportion_closed.csv"))
wide = closed.pivot(index="sample", columns="celltype", values="frac")
wide = wide[wide.index.str.endswith("T")].copy()
wide.index = wide.index.str.replace(r"T$", "", regex=True)
clx = pd.read_excel(os.path.join(ROOT, "data_manifest", "GSE160269_SupplData1_clinical.xlsx"),
                    sheet_name="S1a_scCohort", skiprows=1)
clx.columns = [str(c).strip() for c in clx.columns]
clx = clx.rename(columns={clx.columns[0]: "patient"})
clx["patient"] = clx["patient"].astype(str).str.strip()
clx["stage"] = pd.to_numeric(clx["Pathologic stage"].map({"I": 1, "II": 2, "III": 3}), errors="coerce")
for k, col in [("T", "T stage"), ("N", "N stage")]:
    clx[k] = pd.to_numeric(clx[col].astype(str).str.replace(r"[^0-9]", "", regex=True), errors="coerce")
clx = clx.set_index("patient")

key = ["Tip_angiogenic", "IFNgamma_EC", "Matrix_EC", "Venous", "Capillary", "Arterial"]
rows = []
for var in ["stage", "T", "N"]:
    rr = []
    for s in key:
        d = pd.DataFrame({"x": wide[s].reindex(clx.index), "v": clx[var]}).dropna()
        r, pv = stats.spearmanr(d["x"], d["v"])
        rr.append(dict(state=s, variable=var, n=len(d), rho=r, p=pv))
    tmp = pd.DataFrame(rr)
    tmp["fdr_bh"] = multipletests(tmp["p"], method="fdr_bh")[1]
    rows.append(tmp)
stage = pd.concat(rows, ignore_index=True)
stage.to_csv(os.path.join(OUT, "stage_association_closed.csv"), index=False)
print("\n[stage] 闭合口径下分期关联\n" +
      stage[["state", "variable", "n", "rho", "p", "fdr_bh"]].round(4).to_string(index=False), flush=True)
print(f"[stage] 最小 FDR = {stage.fdr_bh.min():.4f}; 最小名义 P = {stage.p.min():.4f}", flush=True)

# 癌旁 IFNgamma 具体比例（供文字同步）
norm = closed[(closed.type == "Normal") & (closed.celltype == "IFNgamma_EC")]
print("\n[normal IFNgamma_EC frac]", norm[["sample", "n", "n_total", "frac"]].to_string(index=False), flush=True)

with open(os.path.join(OUT, "ph_and_stage_summary.json"), "w", encoding="utf-8") as f:
    json.dump({"cox": {"n": len(df), "events": int(cl.E.sum()), "HR_per_SD": hr,
                       "ci95": [lo, hi], "p": p, "schoenfeld_p": ph_p},
               "stage_min_fdr": float(stage.fdr_bh.min()),
               "stage_min_p": float(stage.p.min()),
               "nominal_hits": stage.loc[stage.p < 0.05].to_dict("records")},
              f, ensure_ascii=False, indent=2)
print("\n[done] 38", flush=True)
