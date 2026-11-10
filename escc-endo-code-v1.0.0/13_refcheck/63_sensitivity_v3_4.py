# -*- coding: utf-8 -*-
"""
v3.4 sensitivity analyses requested by the external reviewer report.

All inputs are local, already-used public data. No new data, no experiments.

Outputs (13_refcheck/out_sensitivity/):
  S11_leave_one_normal_out.csv      LOO over the 4 adjacent-normal libraries
  S11_retain_529_cluster.csv        composition with the 529-cell cluster retained
  S11_effect_sizes_bootstrap.csv    Hodges-Lehmann shift + patient-level bootstrap CI
  S11_four_adjacent_normal_values.csv  the 4 actual normal-library fractions
  S11_exact_vs_asymptotic.csv       exact vs asymptotic Mann-Whitney + BH
  S11_summary.json                  numbers quoted in the manuscript
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SRC = os.path.join(BASE, "08_scRNA", "endothelial", "celltype_proportion_by_sample.csv")
OUT = os.path.join(BASE, "13_refcheck", "out_sensitivity")
os.makedirs(OUT, exist_ok=True)
SEED = 20261108
NBOOT = 2000

STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
          "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"]
KEY = ["Venous", "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC"]
raw = pd.read_csv(SRC)
raw["celltype"] = raw["celltype"].astype(str)
meta = raw.groupby("sample")["type"].first()


def closed(wanted):
    """Closed composition over `wanted` cell types; samples with zero total dropped."""
    sub = raw[raw["celltype"].isin(wanted)]
    cnt = sub.pivot_table(index="sample", columns="celltype", values="n",
                          aggfunc="sum", fill_value=0)
    for c in wanted:
        if c not in cnt.columns:
            cnt[c] = 0
    cnt = cnt[wanted]
    tot = cnt.sum(axis=1)
    keep = tot > 0
    frac = cnt.loc[keep].div(tot.loc[keep], axis=0)
    ty = meta.loc[frac.index]
    return frac, ty


def group_tests(frac, ty, states):
    rows = []
    for s in states:
        t = frac.loc[ty == "Tumor", s].to_numpy()
        n = frac.loc[ty == "Normal", s].to_numpy()
        p_as = stats.mannwhitneyu(t, n, alternative="two-sided", method="asymptotic").pvalue
        p_ex = stats.mannwhitneyu(t, n, alternative="two-sided", method="exact").pvalue
        rows.append(dict(state=s, n_tumor=len(t), n_normal=len(n),
                         med_tumor=float(np.median(t)), med_normal=float(np.median(n)),
                         median_diff=float(np.median(t) - np.median(n)),
                         p_asymptotic=float(p_as), p_exact=float(p_ex)))
    df = pd.DataFrame(rows)
    df["q_asymptotic"] = multipletests(df.p_asymptotic, method="fdr_bh")[1]
    df["q_exact"] = multipletests(df.p_exact, method="fdr_bh")[1]
    df["sig_asymptotic"] = df.q_asymptotic < 0.05
    df["sig_exact"] = df.q_exact < 0.05
    return df


# ---------------------------------------------------------------- 0. baseline
base_frac, base_ty = closed(STATES)
base = group_tests(base_frac, base_ty, STATES)
base.to_csv(os.path.join(OUT, "S11_baseline_group_tests.csv"), index=False)

# ---------------------------------------------------------------- 1. LOO
normals = sorted(base_frac.index[base_ty == "Normal"])
loo_rows = []
for drop in normals:
    f = base_frac.drop(index=drop)
    t = base_ty.drop(index=drop)
    g = group_tests(f, t, STATES)
    for s in KEY:
        r = g[g.state == s].iloc[0]
        b = base[base.state == s].iloc[0]
        loo_rows.append(dict(dropped_normal=drop, state=s,
                             n_normal=int(r.n_normal),
                             med_normal=float(r.med_normal),
                             med_tumor=float(r.med_tumor),
                             median_diff=float(r.median_diff),
                             baseline_diff=float(b.median_diff),
                             direction_preserved=bool(np.sign(r.median_diff) == np.sign(b.median_diff)),
                             p_asymptotic=float(r.p_asymptotic),
                             q_asymptotic=float(r.q_asymptotic),
                             still_significant=bool(r.sig_asymptotic)))
loo = pd.DataFrame(loo_rows)
loo.to_csv(os.path.join(OUT, "S11_leave_one_normal_out.csv"), index=False)

# ---------------------------------------------------------------- 2. retain 529
with529 = STATES + ["Contam_fibroblast"]
f529, t529 = closed(with529)
g529 = group_tests(f529, t529, STATES)
cmp529 = g529[["state", "med_tumor", "med_normal", "median_diff",
               "p_asymptotic", "q_asymptotic", "sig_asymptotic"]].copy()
cmp529.columns = ["state", "med_tumor_529", "med_normal_529", "median_diff_529",
                  "p_529", "q_529", "sig_529"]
cmp529 = cmp529.merge(
    base[["state", "med_tumor", "med_normal", "median_diff", "p_asymptotic",
          "q_asymptotic", "sig_asymptotic"]].rename(
        columns={"med_tumor": "med_tumor_base", "med_normal": "med_normal_base",
                 "median_diff": "median_diff_base", "p_asymptotic": "p_base",
                 "q_asymptotic": "q_base", "sig_asymptotic": "sig_base"}),
    on="state")
cmp529["diff_change_pp"] = (cmp529.median_diff_529 - cmp529.median_diff_base) * 100
cmp529["direction_preserved"] = np.sign(cmp529.median_diff_529) == np.sign(cmp529.median_diff_base)
cmp529["sig_unchanged"] = cmp529.sig_529 == cmp529.sig_base
cmp529.to_csv(os.path.join(OUT, "S11_retain_529_cluster.csv"), index=False)

# ---------------------------------------------------------------- 3. effect sizes
rng = np.random.default_rng(SEED)
es_rows = []
for s in STATES:
    t = base_frac.loc[base_ty == "Tumor", s].to_numpy()
    n = base_frac.loc[base_ty == "Normal", s].to_numpy()
    diffs = (t[:, None] - n[None, :]).ravel()          # all 60 x 4 pairwise differences
    hl = float(np.median(diffs))                        # Hodges-Lehmann two-sample shift
    boot = np.empty(NBOOT)
    for b in range(NBOOT):
        tb = rng.choice(t, size=len(t), replace=True)
        nb = rng.choice(n, size=len(n), replace=True)
        boot[b] = np.median(tb) - np.median(nb)
    es_rows.append(dict(state=s,
                        med_tumor_pct=100 * float(np.median(t)),
                        med_normal_pct=100 * float(np.median(n)),
                        median_diff_pp=100 * float(np.median(t) - np.median(n)),
                        hodges_lehmann_pp=100 * hl,
                        hl_ci95_low_pp=100 * float(np.percentile(boot, 2.5)),
                        hl_ci95_high_pp=100 * float(np.percentile(boot, 97.5)),
                        bootstrap_median_pp=100 * float(np.median(boot)),
                        n_bootstrap=NBOOT))
es = pd.DataFrame(es_rows)
es.to_csv(os.path.join(OUT, "S11_effect_sizes_bootstrap.csv"), index=False)

# ---------------------------------------------------------------- 4. the four normals
four = base_frac.loc[normals, KEY].copy() * 100
four.insert(0, "library", four.index)
four.insert(1, "n_endothelial_cells",
            [int(raw[(raw["sample"] == s) & (raw["celltype"].isin(STATES))]["n"].sum())
             for s in normals])
four = four.round(4)
four.to_csv(os.path.join(OUT, "S11_four_adjacent_normal_values.csv"), index=False)

# ---------------------------------------------------------------- 5. exact vs asymptotic
ex = base[["state", "p_asymptotic", "p_exact", "q_asymptotic", "q_exact",
           "sig_asymptotic", "sig_exact"]].copy()
ex["judgement_flipped"] = ex.sig_asymptotic != ex.sig_exact
ex.to_csv(os.path.join(OUT, "S11_exact_vs_asymptotic.csv"), index=False)

# ---------------------------------------------------------------- summary
def k(df, state, col):
    return float(df.loc[df.state == state, col].iloc[0])

summary = {
    "seed": SEED, "n_bootstrap": NBOOT,
    "n_normal_libraries": len(normals), "normal_libraries": normals,
    "baseline": {s: {"med_normal_pct": 100 * k(base, s, "med_normal"),
                     "med_tumor_pct": 100 * k(base, s, "med_tumor"),
                     "median_diff_pp": 100 * k(base, s, "median_diff"),
                     "p_asymptotic": k(base, s, "p_asymptotic"),
                     "q_asymptotic": k(base, s, "q_asymptotic")} for s in KEY},
    "loo_direction_preserved_all": bool(loo.direction_preserved.all()),
    "loo_still_significant_all": bool(loo.still_significant.all()),
    "loo_n_tests": int(len(loo)),
    "loo_q_range": [float(loo.q_asymptotic.min()), float(loo.q_asymptotic.max())],
    "loo_diff_ranges_pp": {s: [100 * float(loo.loc[loo.state == s, "median_diff"].min()),
                               100 * float(loo.loc[loo.state == s, "median_diff"].max())]
                           for s in KEY},
    "retain_529": {
        "max_abs_diff_change_pp": float(cmp529.diff_change_pp.abs().max()),
        "direction_preserved_all": bool(cmp529.direction_preserved.all()),
        "sig_unchanged_all": bool(cmp529.sig_unchanged.all()),
        "venous_diff_change_pp": float(cmp529.loc[cmp529.state == "Venous", "diff_change_pp"].iloc[0]),
        "tip_diff_change_pp": float(cmp529.loc[cmp529.state == "Tip_angiogenic", "diff_change_pp"].iloc[0]),
    },
    "effect_sizes": {s: {"median_diff_pp": float(es.loc[es.state == s, "median_diff_pp"].iloc[0]),
                         "hodges_lehmann_pp": float(es.loc[es.state == s, "hodges_lehmann_pp"].iloc[0]),
                         "ci95_pp": [float(es.loc[es.state == s, "hl_ci95_low_pp"].iloc[0]),
                                     float(es.loc[es.state == s, "hl_ci95_high_pp"].iloc[0])]}
                     for s in KEY},
    "four_normals": four.to_dict("records"),
    "exact_flipped_states": ex.loc[ex.judgement_flipped, "state"].tolist(),
    "capillary_q_asymptotic": k(base, "Capillary", "q_asymptotic"),
    "capillary_q_exact": k(base, "Capillary", "q_exact"),
}
with open(os.path.join(OUT, "S11_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)

print("=== LOO ===")
print(loo[["dropped_normal", "state", "median_diff", "q_asymptotic",
           "direction_preserved", "still_significant"]].round(4).to_string(index=False))
print("\n=== retain 529 ===")
print(cmp529[["state", "median_diff_base", "median_diff_529", "diff_change_pp",
              "q_base", "q_529", "sig_base", "sig_529"]].round(4).to_string(index=False))
print("\n=== effect sizes ===")
print(es[es.state.isin(KEY)].round(2).to_string(index=False))
print("\n=== four normals (pct) ===")
print(four.to_string(index=False))
print("\n=== exact vs asymptotic ===")
print(ex.round(6).to_string(index=False))
print("\n=== summary ===")
print(json.dumps(summary, ensure_ascii=False, indent=2))
