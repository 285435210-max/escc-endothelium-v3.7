# -*- coding: utf-8 -*-
"""
v3.5 statistical redesign requested by the external reviewer (Major comment 2).

GSE160269: 60 tumor libraries + 4 adjacent-normal libraries from 60 individuals.
Four tumors (P126T, P127T, P128T, P130T) are paired with the four normal libraries
(P126N, P127N, P128N, P130N), so a pooled rank-sum test over 60 vs 4 violates
independence. Primary analysis is therefore 56 UNMATCHED tumors vs 4 normals.
The four paired tumors are presented separately as paired directions only.
The original 60 vs 4 comparison is retained as a sensitivity analysis.

Also addresses:
  - exact Mann-Whitney as the primary small-sample test (asymptotic as sensitivity)
  - explicit BH test families, including the 16 leave-one-normal-out comparisons
  - Hodges-Lehmann shift and patient-level bootstrap CI under the primary design
  - equal-cell-count downsampling recomputed under 56 vs 4

Outputs: 13_refcheck/out_stats_v35/S12_*.csv and S12_summary.json
"""
import os, json, re
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SRC = os.path.join(BASE, "08_scRNA", "endothelial", "celltype_proportion_by_sample.csv")
OUT = os.path.join(BASE, "13_refcheck", "out_stats_v35")
os.makedirs(OUT, exist_ok=True)
SEED = 20261109
NBOOT = 2000
NDOWN = 200

STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
          "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"]
KEY = ["Venous", "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC"]
DISPLAY = {"Arterial": "Arterial", "Capillary": "Capillary", "Venous": "Venous",
           "Venous_inflamed": "Inflamed venule", "Tip_angiogenic": "Tip-like",
           "Matrix_EC": "Matrix-remodeling", "IFNgamma_EC": "Interferon-responsive",
           "Proliferating": "Proliferating", "Lymphatic": "Lymphatic"}

raw = pd.read_csv(SRC)
raw["celltype"] = raw["celltype"].astype(str)
meta = raw.groupby("sample")["type"].first()
rng = np.random.default_rng(SEED)


def closed(wanted=STATES):
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
    return frac, ty, cnt.loc[keep]


def tests(frac, ty, tumor_samples=None, normal_samples=None, states=STATES):
    """Exact and asymptotic Mann-Whitney; BH family = `states` (9 by default)."""
    if tumor_samples is None:
        tumor_samples = frac.index[ty == "Tumor"]
    if normal_samples is None:
        normal_samples = frac.index[ty == "Normal"]
    rows = []
    for s in states:
        t = frac.loc[tumor_samples, s].to_numpy()
        n = frac.loc[normal_samples, s].to_numpy()
        p_as = stats.mannwhitneyu(t, n, alternative="two-sided", method="asymptotic").pvalue
        p_ex = stats.mannwhitneyu(t, n, alternative="two-sided", method="exact").pvalue
        diffs = (t[:, None] - n[None, :]).ravel()
        rows.append(dict(state=s, display=DISPLAY[s],
                         n_tumor=len(t), n_normal=len(n),
                         med_tumor_pct=100 * float(np.median(t)),
                         med_normal_pct=100 * float(np.median(n)),
                         median_diff_pp=100 * float(np.median(t) - np.median(n)),
                         hodges_lehmann_pp=100 * float(np.median(diffs)),
                         p_exact=float(p_ex), p_asymptotic=float(p_as)))
    df = pd.DataFrame(rows)
    df["q_exact"] = multipletests(df.p_exact, method="fdr_bh")[1]
    df["q_asymptotic"] = multipletests(df.p_asymptotic, method="fdr_bh")[1]
    df["sig_exact"] = df.q_exact < 0.05
    df["sig_asymptotic"] = df.q_asymptotic < 0.05
    return df


frac, ty, cnt = closed()
all_tumors = sorted(frac.index[ty == "Tumor"])
all_normals = sorted(frac.index[ty == "Normal"])
normal_patients = sorted({re.sub(r"N$", "", s) for s in all_normals})
paired_tumors = [p + "T" for p in normal_patients]
unmatched_tumors = [s for s in all_tumors if s not in paired_tumors]
assert len(all_tumors) == 60 and len(all_normals) == 4
assert len(paired_tumors) == 4 and all(s in all_tumors for s in paired_tumors)
assert len(unmatched_tumors) == 56

# ---------------------------------------------------------------- 1. pairing map
pair_rows = []
for s in frac.index:
    pat = re.sub(r"[TN]$", "", s)
    t = ty[s]
    pair_rows.append(dict(sample=s, patient=pat, type=t,
                          paired_with_normal=(s in paired_tumors),
                          has_paired_tumor=(t == "Normal" and pat + "T" in all_tumors),
                          n_endothelial_cells=int(cnt.loc[s].sum()),
                          included_in_primary=(s in unmatched_tumors or t == "Normal")))
pairing = pd.DataFrame(pair_rows).sort_values(["type", "sample"])
pairing.to_csv(os.path.join(OUT, "S12_sample_pairing.csv"), index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- 2. primary 56 vs 4
primary = tests(frac, ty, tumor_samples=unmatched_tumors, normal_samples=all_normals)

# bootstrap CI for the median difference under the primary design
boot_rows = []
for s in STATES:
    t = frac.loc[unmatched_tumors, s].to_numpy()
    n = frac.loc[all_normals, s].to_numpy()
    b = np.empty(NBOOT)
    for i in range(NBOOT):
        tb = rng.choice(t, size=len(t), replace=True)
        nb = rng.choice(n, size=len(n), replace=True)
        b[i] = np.median(tb) - np.median(nb)
    boot_rows.append(dict(state=s,
                          ci95_low_pp=100 * float(np.percentile(b, 2.5)),
                          ci95_high_pp=100 * float(np.percentile(b, 97.5)),
                          bootstrap_median_pp=100 * float(np.median(b)),
                          n_bootstrap=NBOOT))
primary = primary.merge(pd.DataFrame(boot_rows), on="state")
primary.to_csv(os.path.join(OUT, "S12_primary_56v4.csv"), index=False)

# ---------------------------------------------------------------- 3. paired directions
pr = []
for pat in normal_patients:
    ts, ns = pat + "T", pat + "N"
    for s in STATES:
        ft, fn = float(frac.loc[ts, s]), float(frac.loc[ns, s])
        pr.append(dict(patient=pat, tumor_sample=ts, normal_sample=ns,
                       state=s, display=DISPLAY[s],
                       tumor_pct=100 * ft, normal_pct=100 * fn,
                       tumor_minus_normal_pp=100 * (ft - fn),
                       n_tumor_cells=int(cnt.loc[ts].sum()),
                       n_normal_cells=int(cnt.loc[ns].sum())))
paired = pd.DataFrame(pr)
paired.to_csv(os.path.join(OUT, "S12_paired_4pairs.csv"), index=False)

# ---------------------------------------------------------------- 4. sensitivity 60 vs 4
sens = tests(frac, ty, tumor_samples=all_tumors, normal_samples=all_normals)
sens.to_csv(os.path.join(OUT, "S12_sensitivity_60v4.csv"), index=False)

# ---------------------------------------------------------------- 5. LOO with explicit BH families
loo_rows = []
for drop in all_normals:
    nd = [s for s in all_normals if s != drop]
    g = tests(frac, ty, tumor_samples=unmatched_tumors, normal_samples=nd)
    for s in KEY:
        r = g[g.state == s].iloc[0]
        b = primary[primary.state == s].iloc[0]
        loo_rows.append(dict(dropped_normal=drop, state=s, display=DISPLAY[s],
                             n_normal=int(r.n_normal),
                             median_diff_pp=float(r.median_diff_pp),
                             baseline_diff_pp=float(b.median_diff_pp),
                             direction_preserved=bool(np.sign(r.median_diff_pp) ==
                                                       np.sign(b.median_diff_pp)),
                             p_exact=float(r.p_exact),
                             q_within_9states=float(r.q_exact),
                             sig_within_9states=bool(r.q_exact < 0.05)))
loo = pd.DataFrame(loo_rows)
# joint BH across the 16 prespecified exploratory comparisons (4 drops x 4 states)
loo["q_joint_16"] = multipletests(loo.p_exact, method="fdr_bh")[1]
loo["sig_joint_16"] = loo.q_joint_16 < 0.05
# BH within each drop over the 4 prespecified key states
loo["q_within_4key"] = np.nan
for d, idx in loo.groupby("dropped_normal").groups.items():
    loo.loc[idx, "q_within_4key"] = multipletests(loo.loc[idx, "p_exact"],
                                                  method="fdr_bh")[1]
loo.to_csv(os.path.join(OUT, "S12_loo_bh_family.csv"), index=False)

# ---------------------------------------------------------------- 6. exact vs asymptotic
ex = primary[["state", "display", "p_exact", "q_exact", "sig_exact",
              "p_asymptotic", "q_asymptotic", "sig_asymptotic"]].copy()
ex["judgement_flipped"] = ex.sig_exact != ex.sig_asymptotic
ex["design"] = "56 unmatched tumors vs 4 normals"
ex.to_csv(os.path.join(OUT, "S12_exact_vs_asymptotic.csv"), index=False)

# ---------------------------------------------------------------- 7. downsampling
def downsample(tumor_samples, normal_samples, N, nrep=NDOWN):
    keep = [s for s in tumor_samples + normal_samples if cnt.loc[s].sum() >= N]
    kt = [s for s in keep if s in tumor_samples]
    kn = [s for s in keep if s in normal_samples]
    if len(kn) < 2 or len(kt) < 2:
        return None
    pv = np.full((nrep, len(STATES)), np.nan)
    ti = [keep.index(s) for s in kt]
    ni = [keep.index(s) for s in kn]
    counts = cnt.loc[keep, STATES].to_numpy(dtype=int)
    for b in range(nrep):
        # one draw per sample per replicate, reused across all 9 states
        v = np.array([rng.multivariate_hypergeometric(counts[i], N)
                      for i in range(len(keep))], dtype=float) / N
        a, c = v[ti, :], v[ni, :]
        for j, st in enumerate(STATES):
            pv[b, j] = stats.mannwhitneyu(a[:, j], c[:, j], alternative="two-sided",
                                          method="asymptotic").pvalue
    q = np.apply_along_axis(lambda x: multipletests(np.nan_to_num(x, nan=1.0),
                                                    method="fdr_bh")[1], 1, pv)
    out = []
    for j, st in enumerate(STATES):
        out.append(dict(state=st, display=DISPLAY[st], cells_per_sample=N,
                        n_tumor=len(kt), n_normal=len(kn), n_repeats=nrep,
                        fdr_pass_rate=float(np.mean(q[:, j] < 0.05))))
    return pd.DataFrame(out)


ds_rows = []
for label, ts, ns in [("primary 56 vs 4", unmatched_tumors, all_normals),
                      ("sensitivity 60 vs 4", all_tumors, all_normals)]:
    for N in [20, 30]:
        d = downsample(ts, ns, N)
        if d is not None:
            d.insert(0, "design", label)
            ds_rows.append(d)
ds = pd.concat(ds_rows, ignore_index=True)
ds.to_csv(os.path.join(OUT, "S12_downsampling.csv"), index=False)

# ---------------------------------------------------------------- summary
def k(df, state, col):
    return float(df.loc[df.state == state, col].iloc[0])


summary = {
    "design": {
        "n_tumor_libraries": len(all_tumors),
        "n_normal_libraries": len(all_normals),
        "normal_libraries": all_normals,
        "paired_patients": normal_patients,
        "paired_tumors": paired_tumors,
        "n_unmatched_tumors": len(unmatched_tumors),
        "primary_test": "exact two-sided Mann-Whitney, 56 unmatched tumors vs 4 normals",
        "bh_family_primary": "9 endothelial states within the primary comparison",
        "paired_handling": "4 paired tumors shown as directions only, no P value",
        "sensitivity": "60 tumors vs 4 normals (original pooling) recomputed with exact test",
        "loo_family": "4 drops x 4 prespecified states = 16 exploratory tests; "
                      "BH applied jointly across 16, and separately within each drop",
    },
    "primary": {s: {
        "med_normal_pct": k(primary, s, "med_normal_pct"),
        "med_tumor_pct": k(primary, s, "med_tumor_pct"),
        "median_diff_pp": k(primary, s, "median_diff_pp"),
        "hodges_lehmann_pp": k(primary, s, "hodges_lehmann_pp"),
        "bootstrap_ci95_pp": [k(primary, s, "ci95_low_pp"), k(primary, s, "ci95_high_pp")],
        "p_exact": k(primary, s, "p_exact"),
        "q_exact": k(primary, s, "q_exact"),
        "q_asymptotic": k(primary, s, "q_asymptotic"),
    } for s in KEY},
    "primary_sig_states": primary.loc[primary.sig_exact, "state"].tolist(),
    "exact_flipped_states": ex.loc[ex.judgement_flipped, "state"].tolist(),
    "loo": {
        "n_tests": int(len(loo)),
        "direction_preserved_all": bool(loo.direction_preserved.all()),
        "sig_joint16_all": bool(loo.sig_joint_16.all()),
        "sig_within9_all": bool(loo.sig_within_9states.all()),
        "q_joint16_range": [float(loo.q_joint_16.min()), float(loo.q_joint_16.max())],
        "q_within9_range": [float(loo.q_within_9states.min()),
                            float(loo.q_within_9states.max())],
    },
    "paired_direction": {
        s: int((paired[paired.state == s].tumor_minus_normal_pp > 0).sum())
        for s in KEY
    },
    "downsampling": {
        "%s_%d" % (r.design, r.cells_per_sample): 1.0
        for r in ds.itertuples()
    },
}
# explicit downsampling pass rates by state
summary["downsampling_pass_rates"] = {
    "%s/%d/%s" % (r.design, r.cells_per_sample, r.state): round(r.fdr_pass_rate, 4)
    for r in ds.itertuples() if r.state in KEY
}
summary["four_normal_venous_pct"] = [
    round(100 * float(frac.loc[s, "Venous"]), 2) for s in all_normals]

with open(os.path.join(OUT, "S12_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)

print("=== pairing ===")
print("paired tumors:", paired_tumors, "unmatched:", len(unmatched_tumors))
print("\n=== PRIMARY 56 vs 4 (exact) ===")
print(primary[["state", "med_normal_pct", "med_tumor_pct", "median_diff_pp",
               "hodges_lehmann_pp", "ci95_low_pp", "ci95_high_pp",
               "p_exact", "q_exact", "sig_exact"]].round(4).to_string(index=False))
print("\n=== SENSITIVITY 60 vs 4 ===")
print(sens[["state", "p_exact", "q_exact", "sig_exact"]].round(6).to_string(index=False))
print("\n=== PAIRED DIRECTION (tumor - normal, pp) ===")
print(paired[paired.state.isin(KEY)].pivot(index="patient", columns="state",
                                           values="tumor_minus_normal_pp").round(2).to_string())
print("\n=== LOO BH families ===")
print(loo[["dropped_normal", "state", "median_diff_pp", "p_exact",
           "q_within_9states", "q_joint_16", "q_within_4key",
           "direction_preserved", "sig_joint_16"]].round(5).to_string(index=False))
print("\n=== EXACT vs ASYMPTOTIC (primary) ===")
print(ex.round(6).to_string(index=False))
print("\n=== DOWNSAMPLING (key states) ===")
print(ds[ds.state.isin(KEY)].round(3).to_string(index=False))
print("\n=== four normal venous % ===", summary["four_normal_venous_pct"])
