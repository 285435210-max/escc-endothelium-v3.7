# -*- coding: utf-8 -*-
"""
Step 2: cross-platform prognostic re-check of the 23-gene panel in TCGA-ESCA.

Pre-specified design:
  * Only primary tumours (-01) with histology "Esophagus Squamous Cell Carcinoma".
  * The 23 coefficients were fitted in GSE53625; NO refitting in TCGA-ESCA.
  * Gene-wise z-score ACROSS SAMPLES within each cohort (same rule as discovery),
    then score = z %*% coef. Primary standardisation uses TUMOUR samples only;
    all-sample standardisation is reported as a sensitivity analysis because the
    earlier B-layer script (38_ph_and_stage.py) used all 358 GSE53625 arrays.
  * Primary endpoint: overall survival. DSS, PFI, all-histology ESCA and an
    age/stage-adjusted model are sensitivity analyses.
  * Negative controls: 200 random 23-gene sets with random +/-1 weights, and
    200 sign-permutations of the original 23 coefficients.
  * lifelines' concordance_index interprets its second argument as a predicted
    EVENT TIME (higher = longer survival). A Cox linear predictor therefore has
    to be negated; an un-negated value is exactly 1 - the correct C-index.

Outputs: 15_tcga_esca/out/TableS10_*.csv, tcga_esca_summary.json
"""
import os
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test, logrank_test
from lifelines.utils import concordance_index

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
TCGA = os.path.join(ROOT, "15_tcga_esca")
RAW = os.path.join(TCGA, "raw")
OUT = os.path.join(TCGA, "out")
G536 = r"D:\龙虾文档\2026-10-21-20-33-54\data\GSE53625"
SEED = 20261107
NPERM = 200
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(SEED)


def wr(name, df):
    df.to_csv(os.path.join(OUT, name), index=False, encoding="utf-8-sig")
    print("  ->", name, df.shape, flush=True)


def make_score(expr, coef, z_on=None):
    """expr: genes x samples. z_on: columns used to estimate gene mean/sd."""
    use = [g for g in coef.index if g in expr.index]
    E = expr.loc[use].T.astype(float)                        # samples x genes
    ref = E.loc[z_on] if z_on is not None else E
    sd = ref.std(ddof=1).replace(0, np.nan)
    z = ((E - ref.mean()) / sd).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return pd.Series(z.to_numpy() @ coef[use].to_numpy(), index=E.index, name="score"), use


def cox_one(score, T, E, label, extra=None):
    d = pd.DataFrame({"score": score, "T": T, "E": E}).dropna()
    d = d[(d["T"] > 0) & d["E"].isin([0, 1])]
    if len(d) < 20 or d["E"].sum() < 10:
        return None, d
    d["score_z"] = (d["score"] - d["score"].mean()) / d["score"].std(ddof=1)
    cols = ["T", "E", "score_z"]
    if extra is not None:
        for k, v in extra.items():
            d[k] = pd.to_numeric(v, errors="coerce").reindex(d.index)
        cols += list(extra.keys())
        d = d.dropna(subset=cols)
    fit = CoxPHFitter().fit(d[cols], duration_col="T", event_col="E")
    beta = float(fit.params_["score_z"])
    ci = np.exp(fit.confidence_intervals_.loc["score_z"]).tolist()
    hr, lo, hi = float(np.exp(beta)), ci[0], ci[1]
    p = float(fit.summary.loc["score_z", "p"])
    # lifelines: higher predicted score = LONGER survival -> negate risk score
    c_cox = float(concordance_index(d["T"], -beta * d["score_z"], d["E"]))
    c_raw = float(concordance_index(d["T"], -d["score_z"], d["E"]))
    try:
        ph = proportional_hazard_test(fit, d[cols], time_transform="rank")
        ph_p = float(ph.summary.loc["score_z", "p"])
    except Exception:
        ph_p = np.nan
    hi_sc = d["score_z"] > d["score_z"].median()
    lr = logrank_test(d.loc[hi_sc, "T"], d.loc[~hi_sc, "T"],
                      d.loc[hi_sc, "E"], d.loc[~hi_sc, "E"])
    res = dict(model=label, n=len(d), events=int(d["E"].sum()),
               HR_per_SD=hr, lower95=lo, upper95=hi, p=p,
               cindex_cox=c_cox, cindex_raw=c_raw, schoenfeld_p=ph_p,
               logrank_p=float(lr.p_value), beta=beta)
    return res, d


def random_null(expr, T, E, mode="genes", z_on=None):
    genes = np.array(expr.index)
    coef0 = COEF.loc[COEF.index.intersection(genes)]
    vals = []
    for _ in range(NPERM):
        if mode == "genes":
            g = rng.choice(genes, size=len(coef0), replace=False)
            w = pd.Series(rng.choice([-1.0, 1.0], size=len(coef0)), index=g)
        else:
            w = coef0 * rng.choice([-1.0, 1.0], size=len(coef0))
        s, _ = make_score(expr, w, z_on=z_on)
        r, _ = cox_one(s, T.reindex(s.index), E.reindex(s.index), "null")
        if r is not None:
            vals.append(r["cindex_cox"])
    return np.array(vals)


# ---------------------------------------------------------------- 0. coefficients
COEF = pd.read_csv(os.path.join(ROOT, "09_de", "panel", "panel_coefficients.csv"),
                   encoding="utf-8-sig", index_col=0)["coef"]
COEF = COEF[COEF != 0].astype(float)
print("[0] panel genes:", len(COEF), flush=True)

# ---------------------------------------------------------------- 1. TCGA-ESCA
print("[1] loading TCGA-ESCA ...", flush=True)
ex = pd.read_csv(os.path.join(RAW, "HiSeqV2.gz"), sep="\t", index_col=0)
ex = ex.groupby(ex.index).mean()
ex.columns = [c.strip() for c in ex.columns]
cl = pd.read_csv(os.path.join(RAW, "ESCA_clinicalMatrix"), sep="\t", low_memory=False)
cl = cl.set_index("sampleID")
sv = pd.read_csv(os.path.join(RAW, "ESCA_survival.txt"), sep="\t").set_index("sample")

tum = [c for c in ex.columns if c.endswith("-01")]
escc = [c for c in tum
        if c in cl.index and "Squamous" in str(cl.loc[c, "histological_type"])]
print("    expression", ex.shape, "| primary tumours", len(tum),
      "| ESCC primary", len(escc), flush=True)

score, matched = make_score(ex, COEF, z_on=escc)
score = score.reindex(escc)
print("    matched genes", len(matched), "/", len(COEF), flush=True)

T_os = pd.to_numeric(sv["OS.time"], errors="coerce")
E_os = pd.to_numeric(sv["OS"], errors="coerce")
res_main, d_main = cox_one(score, T_os.reindex(score.index), E_os.reindex(score.index),
                           "TCGA-ESCA ESCC, OS, tumour-only z-score (PRIMARY)")
assert res_main is not None
print("[1] PRIMARY:", json.dumps({k: round(v, 4) if isinstance(v, float) else v
                                 for k, v in res_main.items()}, ensure_ascii=False), flush=True)

samples = d_main[["score", "score_z", "T", "E"]].copy()
samples.insert(0, "sample", samples.index)
samples["histology"] = cl["histological_type"].reindex(samples.index).values
samples["age"] = pd.to_numeric(cl["age_at_initial_pathologic_diagnosis"], errors="coerce") \
                   .reindex(samples.index).values
samples["clinical_stage"] = cl["clinical_stage"].reindex(samples.index).values
wr("TableS10_TCGA_ESCC_samples.csv", samples.round(5))

# ---------------------------------------------------------------- 2. GSE53625 identical pipeline
print("[2] re-scoring GSE53625 with the identical function ...", flush=True)
gx = pd.read_csv(os.path.join(G536, "expr_gene.csv"), index_col=0)
gx.columns = [c.strip() for c in gx.columns]
gx = gx.groupby(gx.index).mean()
grp = pd.read_csv(os.path.join(G536, "group.csv"))
grp["sample"] = grp["sample"].astype(str).str.strip()
clin = pd.read_csv(os.path.join(ROOT, "data_manifest", "GSE53625_clinical_all.csv"),
                   encoding="utf-8-sig")
clin.columns = [c.strip() for c in clin.columns]
clin["sample"] = clin["sample"].astype(str).str.strip()
clin = clin.set_index("sample")
tc = [s for s in grp.loc[grp.group == "Tumor", "sample"].tolist() if s in gx.columns]
g_T = pd.to_numeric(clin["OS.time"], errors="coerce")
g_E = pd.to_numeric(clin["OS"], errors="coerce")

g_score, g_matched = make_score(gx, COEF, z_on=tc)
g_score = g_score.reindex(tc)
res_g536, _ = cox_one(g_score, g_T.reindex(tc), g_E.reindex(tc),
                      "GSE53625 ESCC, OS, tumour-only z-score (PRIMARY)")
assert res_g536 is not None
print("[2] GSE53625 same-pipeline:", json.dumps(
    {k: round(v, 4) if isinstance(v, float) else v for k, v in res_g536.items()},
    ensure_ascii=False), flush=True)

# all-sample z-score sensitivity (reproduces the earlier B-layer script)
g_all, _ = make_score(gx, COEF, z_on=None)
res_g_all, _ = cox_one(g_all.reindex(tc), g_T.reindex(tc), g_E.reindex(tc),
                       "GSE53625 ESCC, OS, all-array z-score (B-layer reproducibility)")
t_all, _ = make_score(ex, COEF, z_on=None)
res_t_all, _ = cox_one(t_all.reindex(escc), T_os.reindex(escc), E_os.reindex(escc),
                       "TCGA-ESCA ESCC, OS, all-sample z-score")
zsens = pd.DataFrame([x for x in [res_g536, res_g_all, res_main, res_t_all] if x])
wr("TableS10_zscore_sensitivity.csv", zsens.round(5))
print("[2] GSE53625 all-array z-score HR/P/ Schoenfeld:",
      round(res_g_all["HR_per_SD"], 3), round(res_g_all["p"], 3),
      round(res_g_all["schoenfeld_p"], 4), flush=True)

# ---------------------------------------------------------------- 3. negative controls
print("[3] negative controls (200 each) ...", flush=True)
rows = []
for mode, lab in [("genes", "random 23 genes, random +/-1 weights"),
                  ("signs", "same 23 genes, random coefficient signs")]:
    v = random_null(ex[escc], T_os, E_os, mode, z_on=escc)
    q = np.quantile(v, [0.025, 0.25, 0.5, 0.75, 0.975])
    obs = res_main["cindex_cox"]
    rows.append(dict(cohort="TCGA-ESCA ESCC", null=lab, n_perm=len(v),
                     median_C=float(np.median(v)), q025=q[0], q975=q[4],
                     min=float(v.min()), max=float(v.max()), observed_C=obs,
                     empirical_P_one_sided=(1 + np.sum(v >= obs)) / (1 + len(v))))
    print("   ", lab, "median", round(float(np.median(v)), 4),
          "95%", np.round(q[[0, 4]], 4),
          "P_emp", round(rows[-1]["empirical_P_one_sided"], 4), flush=True)
v_g = random_null(gx[tc], g_T, g_E, "genes", z_on=tc)
q_g = np.quantile(v_g, [0.025, 0.5, 0.975])
rows.append(dict(cohort="GSE53625", null="random 23 genes, random +/-1 weights",
                 n_perm=len(v_g), median_C=float(np.median(v_g)),
                 q025=q_g[0], q975=q_g[2], min=float(v_g.min()), max=float(v_g.max()),
                 observed_C=res_g536["cindex_cox"],
                 empirical_P_one_sided=(1 + np.sum(v_g >= res_g536["cindex_cox"])) /
                                        (1 + len(v_g))))
null_tab = pd.DataFrame(rows)
wr("TableS10_random_controls.csv", null_tab.round(5))
print("    GSE53625 null median", round(float(np.median(v_g)), 4),
      "P_emp", round(rows[-1]["empirical_P_one_sided"], 4), flush=True)

# ---------------------------------------------------------------- 4. sensitivity analyses
print("[4] sensitivity analyses ...", flush=True)
sens = []
s_a, _ = make_score(ex, COEF, z_on=tum)
r, _ = cox_one(s_a.reindex(tum), T_os.reindex(tum), E_os.reindex(tum),
               "All ESCA primary tumours (ESCC + adenocarcinoma), OS")
sens.append(r)
for ep in ["DSS", "PFI"]:
    Te = pd.to_numeric(sv[ep + ".time"], errors="coerce").reindex(score.index)
    Ee = pd.to_numeric(sv[ep], errors="coerce").reindex(score.index)
    r, _ = cox_one(score, Te, Ee, "TCGA-ESCA ESCC, " + ep)
    sens.append(r)
stage_map = {"Stage I": 1, "Stage IA": 1, "Stage IB": 1,
             "Stage II": 2, "Stage IIA": 2, "Stage IIB": 2,
             "Stage III": 3, "Stage IIIA": 3, "Stage IIIB": 3, "Stage IIIC": 3,
             "Stage IV": 4}
adj = {"age": cl["age_at_initial_pathologic_diagnosis"],
       "stage_num": cl["clinical_stage"].map(stage_map)}
r, _ = cox_one(score, T_os.reindex(score.index), E_os.reindex(score.index),
               "TCGA-ESCA ESCC, OS, age+stage adjusted", extra=adj)
sens.append(r)
all_models = pd.concat([pd.DataFrame([res_g536, res_g_all, res_main, res_t_all]),
                        pd.DataFrame([x for x in sens if x])], ignore_index=True)
wr("TableS10_cox_all_models.csv", all_models.round(5))
for r in sens:
    if r:
        print("   ", r["model"], "n", r["n"], "events", r["events"],
              "HR/SD", round(r["HR_per_SD"], 3), "P", round(r["p"], 3),
              "C", round(r["cindex_cox"], 3), flush=True)

# ---------------------------------------------------------------- 5. cross-cohort table
cross = pd.DataFrame([
    dict(cohort="GSE53625 (discovery; microarray)", platform="Affymetrix array",
         histology="ESCC", n=res_g536["n"], events=res_g536["events"], endpoint="OS",
         HR_per_SD=res_g536["HR_per_SD"],
         CI95="%.3f-%.3f" % (res_g536["lower95"], res_g536["upper95"]),
         p=res_g536["p"], cindex=res_g536["cindex_cox"],
         schoenfeld_p=res_g536["schoenfeld_p"], logrank_p=res_g536["logrank_p"]),
    dict(cohort="TCGA-ESCA (external; RNA-seq)", platform="Illumina RNA-seq (Xena HiSeqV2)",
         histology="ESCC only", n=res_main["n"], events=res_main["events"], endpoint="OS",
         HR_per_SD=res_main["HR_per_SD"],
         CI95="%.3f-%.3f" % (res_main["lower95"], res_main["upper95"]),
         p=res_main["p"], cindex=res_main["cindex_cox"],
         schoenfeld_p=res_main["schoenfeld_p"], logrank_p=res_main["logrank_p"]),
])
wr("TableS10_cross_cohort_comparison.csv", cross.round(5))

prov = pd.DataFrame([
    ["analysis_seed", SEED], ["n_permutations", NPERM],
    ["expression_file", "tcga.xenahubs.net/download/TCGA.ESCA.sampleMap/HiSeqV2.gz"],
    ["clinical_file", "tcga.xenahubs.net/download/TCGA.ESCA.sampleMap/ESCA_clinicalMatrix"],
    ["survival_file", "tcga.xenahubs.net/download/survival/ESCA_survival.txt"],
    ["genes_matched", "%d/%d" % (len(matched), len(COEF))],
    ["endpoint_primary", "OS (OS.time, OS)"],
    ["histology_filter", "Esophagus Squamous Cell Carcinoma, primary tumour (-01) only"],
    ["standardisation_primary", "gene-wise z-score across TUMOUR samples, within each cohort"],
    ["cindex_convention", "lifelines concordance_index: higher predicted score = longer survival; risk score negated"],
    ["coefficients", "GSE53625 full-cohort descriptive L1 model; no refitting in TCGA-ESCA"],
    ["multiple_testing", "OS is the single primary endpoint; all other models are sensitivity analyses"],
], columns=["item", "value"])
wr("TableS10_provenance.csv", prov)

with open(os.path.join(OUT, "tcga_esca_summary.json"), "w", encoding="utf-8") as f:
    json.dump(dict(seed=SEED, n_perm=NPERM, genes_matched=len(matched),
                   primary=res_main, gse53625_same_pipeline=res_g536,
                   gse53625_allarray=res_g_all, tcga_all_sample=res_t_all,
                   sensitivity=sens, null_controls=rows), f, ensure_ascii=False, indent=2)
print("[done] 48 corrected", flush=True)
