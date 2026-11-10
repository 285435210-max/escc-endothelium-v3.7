# -*- coding: utf-8 -*-
"""v3.7 supporting analyses: 7-state stage correlations and chemotherapy scores."""
import os, io, json, re
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
OUT = os.path.join(BASE, "13_refcheck", "out_v37")
os.makedirs(OUT, exist_ok=True)
DATA = r"D:/escc_sc/data"

STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed",
          "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC"]
DISPLAY = {"Arterial": "Arterial", "Capillary": "Capillary", "Venous": "Venous",
           "Venous_inflamed": "Inflamed venule", "Tip_angiogenic": "Tip-like",
           "Matrix_EC": "Matrix-remodeling", "IFNgamma_EC": "Interferon-responsive"}

cl = pd.read_excel(os.path.join(DATA, "GSE160269_SupplData1_clinical.xlsx"),
                   sheet_name="S1a_scCohort", skiprows=1)
cl.columns = [str(c).strip() for c in cl.columns]
cl = cl.rename(columns={cl.columns[0]: "patient"})
cl["patient"] = cl["patient"].astype(str).str.strip()
stage_map = {"I": 1, "II": 2, "III": 3}
cl["stage"] = cl["Pathologic stage"].astype(str).str.strip().map(stage_map)
cl["Tstage"] = pd.to_numeric(cl["T stage"].astype(str).str.replace(r"[^0-9]", "", regex=True),
                             errors="coerce")
cl["Nstage"] = pd.to_numeric(cl["N stage"].astype(str).str.replace(r"[^0-9]", "", regex=True),
                             errors="coerce")

prop = pd.read_csv(os.path.join(DATA, "celltype_proportion_closed.csv"))
prop["patient"] = prop["sample"].astype(str).str.replace(r"T$", "", regex=True)
prop = prop[prop["type"] == "Tumor"]

rows = []
for st in STATES:
    sub = prop[prop["celltype"] == st][["patient", "frac"]].copy()
    # complete grid: absent state = 0
    grid = pd.DataFrame({"patient": cl["patient"]})
    grid = grid.merge(sub, on="patient", how="left").fillna({"frac": 0.0})
    m = grid.merge(cl[["patient", "stage", "Tstage", "Nstage"]], on="patient", how="inner")
    for var in ["stage", "Tstage", "Nstage"]:
        d = m[["frac", var]].dropna()
        if len(d) < 10:
            continue
        rho, p = stats.spearmanr(d["frac"], d[var])
        rows.append(dict(state=st, display=DISPLAY[st], variable=var,
                         n_patients=len(d), rho=round(float(rho), 4),
                         p=float(p)))
res = pd.DataFrame(rows)
res["q_BH_within_variable"] = np.nan
for var, ii in res.groupby("variable").groups.items():
    res.loc[list(ii), "q_BH_within_variable"] = multipletests(
        res.loc[list(ii), "p"], method="fdr_bh")[1]
res["significant"] = res.q_BH_within_variable < 0.05
res.to_csv(os.path.join(OUT, "S14_stage_7states.csv"), index=False, encoding="utf-8-sig")

chemo = pd.read_csv(r"D:/escc_sc/gse225178/out/GSE225178_result.csv", encoding="utf-8-sig")
chemo = chemo.rename(columns={"signature": "score", "wilcox_P": "p_mannwhitney",
                              "AUC": "auc"})
chemo["score_label"] = chemo["score"].map({
    "tip487": "487-gene tip signature",
    "tip_core": "7-gene core score",
    "breach": "7-gene matrix-remodeling score"})
chemo["n_responders"] = chemo["n_resp"]
chemo["n_nonresponders"] = chemo["n_nonresp"]
chemo = chemo[["score", "score_label", "n_responders", "n_nonresponders",
               "mean_resp", "mean_nonresp", "p_mannwhitney", "auc"]]
chemo.to_csv(os.path.join(OUT, "S14_chemo_scores.csv"), index=False, encoding="utf-8-sig")

stage_csv = os.path.join(OUT, "stage_lab_7states.csv")
summary = dict(
    n_states_stage=len(STATES),
    n_stage_tests=len(res),
    min_q=float(res.q_BH_within_variable.min()),
    any_significant=bool(res.significant.any()),
    capillary_stage=res[(res.state == "Capillary") & (res.variable == "stage")][
        ["rho", "p", "q_BH_within_variable"]].to_dict(orient="records"),
    inflamed_stage=res[(res.state == "Venous_inflamed") & (res.variable == "stage")][
        ["rho", "p", "q_BH_within_variable"]].to_dict(orient="records"),
    chemo=chemo[["score", "p_mannwhitney", "auc"]].to_dict(orient="records"),
    r_stage_table_exists=os.path.exists(stage_csv),
)
with io.open(os.path.join(OUT, "S14_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

print("[OK] 7-state stage tests:", len(res), "min q =", round(summary["min_q"], 4),
      "any significant =", summary["any_significant"])
print("[OK] capillary stage:", summary["capillary_stage"])
print("[OK] inflamed venule stage:", summary["inflamed_stage"])
print("[OK] chemo scores:", [(r["score"], round(r["auc"], 3),
                              round(r["p_mannwhitney"], 3)) for r in summary["chemo"]])
