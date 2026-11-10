# -*- coding: utf-8 -*-
"""Tie-aware EXACT label-permutation tests for the v3.6 revision.

External review (v3.5, comment 1): SciPy's method='exact' Mann-Whitney does not
correct for tied ranks, and the closed composition data contain many ties (zeros).
We therefore enumerate every label assignment:

  primary  56 tumours vs 4 adjacent libraries -> C(60,4) = 487,635 assignments
  LOO      56 tumours vs 3 adjacent libraries -> C(59,3) =  32,509 assignments

For each state the Mann-Whitney U statistic is computed from the RAW pairwise
comparison matrix (ties contribute 0.5), so no rank-tie assumption is needed.
Two-sided P is defined symmetrically about the null centre n1*n2/2:
  P = #{|U - centre| >= |U_obs - centre|} / #assignments
A second conventional definition (2 x smaller tail, capped at 1) is also
reported. All arithmetic uses integer 2U to avoid floating-point drift.

Outputs: 13_refcheck/out_stats_v36/S13_*.csv and S13_summary.json
"""
import os, io, json, re, itertools, math
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SRC = os.path.join(BASE, "08_scRNA", "endothelial", "celltype_proportion_by_sample.csv")
S12 = os.path.join(BASE, "13_refcheck", "out_stats_v35")
OUT = os.path.join(BASE, "13_refcheck", "out_stats_v36")
os.makedirs(OUT, exist_ok=True)

STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
          "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"]
KEY4 = ["Venous", "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC"]
KEY5 = KEY4 + ["Venous_inflamed"]
DISPLAY = {"Arterial": "Arterial", "Capillary": "Capillary", "Venous": "Venous",
           "Venous_inflamed": "Inflamed venule", "Tip_angiogenic": "Tip-like",
           "Matrix_EC": "Matrix-remodeling", "IFNgamma_EC": "Interferon-responsive",
           "Proliferating": "Proliferating", "Lymphatic": "Lymphatic"}

raw = pd.read_csv(SRC)
raw["celltype"] = raw["celltype"].astype(str)
meta = raw.groupby("sample")["type"].first()
cnt = (raw[raw["celltype"].isin(STATES)]
       .pivot_table(index="sample", columns="celltype", values="n",
                    aggfunc="sum", fill_value=0))
for c in STATES:
    if c not in cnt.columns:
        cnt[c] = 0
cnt = cnt[STATES]
tot = cnt.sum(axis=1)
frac = cnt.loc[tot > 0].div(tot.loc[tot > 0], axis=0)
ty = meta.loc[frac.index]

all_tumors = sorted(frac.index[ty == "Tumor"])
all_normals = sorted(frac.index[ty == "Normal"])
paired_tumors = [re.sub(r"N$", "T", s) for s in all_normals]
unmatched = [s for s in all_tumors if s not in paired_tumors]
assert len(all_tumors) == 60 and len(all_normals) == 4 and len(unmatched) == 56


def combinations_array(n, k):
    count = int(math.comb(n, k))
    a = np.fromiter(itertools.chain.from_iterable(itertools.combinations(range(n), k)),
                    dtype=np.int16, count=count * k)
    return a.reshape(count, k), count


def permutation_state(x, normal_pos, idx):
    """Return exact symmetric two-sided P and 2*min-tail P for one state."""
    n = len(x)
    k = len(normal_pos)
    # 2*U contributions: 2 if a>b, 1 if tied, 0 otherwise
    T = (x[:, None] > x[None, :]).astype(np.int16) * 2 \
        + (x[:, None] == x[None, :]).astype(np.int16)
    S = T.sum(axis=0).astype(np.int64)                       # 2U if group B were {j} only
    all_u2 = S[idx].sum(axis=1)
    corr = np.zeros(len(idx), dtype=np.int64)
    for a in range(k):
        for b in range(k):
            corr += T[idx[:, a], idx[:, b]]
    u2 = all_u2 - corr                                       # integer 2U, A = complement
    centre2 = (n - k) * k                                    # 2 * n1*n2/2
    obs = np.zeros(n, dtype=bool)
    obs[normal_pos] = True
    obs_idx = np.flatnonzero(np.all(obs[idx], axis=1))
    assert len(obs_idx) == 1, "observed label assignment must be in the enumeration"
    u_obs2 = int(u2[obs_idx[0]])
    d_obs = abs(u_obs2 - centre2)
    p_sym = float(np.count_nonzero(np.abs(u2 - centre2) >= d_obs) / len(u2))
    lo = float(np.count_nonzero(u2 <= u_obs2) / len(u2))
    hi = float(np.count_nonzero(u2 >= u_obs2) / len(u2))
    p_2min = float(min(1.0, 2.0 * min(lo, hi)))
    return dict(u2_obs=u_obs2, u_obs=u_obs2 / 2.0, centre=centre2 / 2.0,
                p_perm_sym=p_sym, p_perm_2min=p_2min,
                n_assignments=len(u2), n_ties=int(np.count_nonzero(np.triu(T == 1, 1))))


def scipy_pair(x, normal_pos, n, k):
    obs = np.zeros(n, dtype=bool)
    obs[normal_pos] = True
    t, nn = x[~obs], x[obs]
    return (float(stats.mannwhitneyu(t, nn, alternative="two-sided", method="exact").pvalue),
            float(stats.mannwhitneyu(t, nn, alternative="two-sided", method="asymptotic").pvalue),
            float(np.median(t) - np.median(nn)) * 100)


# ------------------------------------------------------------------ primary 56 vs 4
samples = unmatched + all_normals
idx, ncomb = combinations_array(60, 4)
assert ncomb == 487635
normal_pos = list(range(56, 60))
rows = []
for st in STATES:
    x = frac.loc[samples, st].to_numpy(dtype=float)
    r = permutation_state(x, normal_pos, idx)
    p_ex, p_as, med_diff = scipy_pair(x, normal_pos, 60, 4)
    r.update(state=st, display=DISPLAY[st], n_tumor=56, n_normal=4,
             median_diff_pp=med_diff, p_exact_scipy=p_ex, p_asymptotic_tiecorrected=p_as)
    rows.append(r)
primary = pd.DataFrame(rows)
primary["q_perm"] = multipletests(primary.p_perm_sym, method="fdr_bh")[1]
primary["q_exact"] = multipletests(primary.p_exact_scipy, method="fdr_bh")[1]
primary["q_asymptotic"] = multipletests(primary.p_asymptotic_tiecorrected, method="fdr_bh")[1]
primary["sig_perm"] = primary.q_perm < 0.05
primary["sig_exact"] = primary.q_exact < 0.05
primary["sig_asymptotic"] = primary.q_asymptotic < 0.05
primary["concordant_all_three"] = (primary.sig_perm == primary.sig_exact) & \
                                  (primary.sig_perm == primary.sig_asymptotic)
primary.to_csv(os.path.join(OUT, "S13_permutation_primary.csv"), index=False,
               encoding="utf-8-sig")

# ------------------------------------------------------------------ pooled 60 vs 4 sensitivity
samples64 = all_tumors + all_normals
idx64, n64 = combinations_array(64, 4)
assert n64 == 635376
pos64 = [60, 61, 62, 63]
rows64 = []
for st in STATES:
    x = frac.loc[samples64, st].to_numpy(dtype=float)
    r = permutation_state(x, pos64, idx64)
    p_ex, p_as, med_diff = scipy_pair(x, pos64, 64, 4)
    r.update(state=st, display=DISPLAY[st], n_tumor=60, n_normal=4,
             median_diff_pp=med_diff, p_exact_scipy=p_ex,
             p_asymptotic_tiecorrected=p_as)
    rows64.append(r)
pooled = pd.DataFrame(rows64)
pooled["q_perm"] = multipletests(pooled.p_perm_sym, method="fdr_bh")[1]
pooled["q_exact"] = multipletests(pooled.p_exact_scipy, method="fdr_bh")[1]
pooled["sig_perm"] = pooled.q_perm < 0.05
pooled["sig_exact"] = pooled.q_exact < 0.05
pooled["concordant"] = pooled.sig_perm == pooled.sig_exact
pooled.to_csv(os.path.join(OUT, "S13_permutation_60v4.csv"), index=False,
              encoding="utf-8-sig")

# ------------------------------------------------------------------ LOO 56 vs 3 (5 states)
loo_rows = []
for dropped in all_normals:
    normals3 = [s for s in all_normals if s != dropped]
    smp = unmatched + normals3
    ix, nc = combinations_array(59, 3)
    assert nc == 32509
    pos = [56, 57, 58]
    for st in KEY5:
        x = frac.loc[smp, st].to_numpy(dtype=float)
        r = permutation_state(x, pos, ix)
        p_ex, p_as, med_diff = scipy_pair(x, pos, 59, 3)
        r.update(dropped_normal=dropped, state=st, display=DISPLAY[st],
                 n_tumor=56, n_normal=3, median_diff_pp=med_diff,
                 baseline_diff_pp=None, p_exact_scipy=p_ex,
                 p_asymptotic_tiecorrected=p_as)
        loo_rows.append(r)
loo = pd.DataFrame(loo_rows)
base = primary.set_index("state").median_diff_pp.to_dict()
loo["baseline_diff_pp"] = loo.state.map(base)
loo["direction_preserved"] = np.sign(loo.median_diff_pp) == np.sign(loo.baseline_diff_pp)
loo16 = loo[loo.state.isin(KEY4)].copy().reset_index(drop=True)
loo16["q_joint_16"] = multipletests(loo16.p_perm_sym, method="fdr_bh")[1]
loo16["q_within_4key"] = np.nan
for d, ii in loo16.groupby("dropped_normal").groups.items():
    loo16.loc[list(ii), "q_within_4key"] = multipletests(
        loo16.loc[list(ii), "p_perm_sym"], method="fdr_bh")[1]
loo16["sig_joint_16"] = loo16.q_joint_16 < 0.05
loo16.to_csv(os.path.join(OUT, "S13_permutation_loo_16.csv"), index=False,
             encoding="utf-8-sig")
inf = loo[loo.state == "Venous_inflamed"].copy().reset_index(drop=True)
inf["q_within_4"] = multipletests(inf.p_perm_sym, method="fdr_bh")[1]
inf["sig_within_4"] = inf.q_within_4 < 0.05
inf.to_csv(os.path.join(OUT, "S13_inflamed_venule_loo.csv"), index=False,
           encoding="utf-8-sig")

# ------------------------------------------------------------------ four paired states
pairs = pd.read_csv(os.path.join(S12, "S12_paired_4pairs.csv"))
pair_rows = []
for st in KEY4:
    sub = pairs[pairs.state == st]
    vals = sub.tumor_minus_normal_pp.to_numpy(dtype=float)
    exp = -1 if st == "Venous" else 1
    pair_rows.append(dict(state=st, display=DISPLAY[st], total_pairs=4,
                          pairs_expected_direction=int(np.sum(vals * exp > 0)),
                          pairs_opposite_direction=int(np.sum(vals * exp < 0)),
                          pairs_tied_at_zero=int(np.sum(vals == 0)),
                          values_pp="; ".join("%.2f" % v for v in vals)))
pd.DataFrame(pair_rows).to_csv(os.path.join(OUT, "S13_paired_four_states.csv"),
                               index=False, encoding="utf-8-sig")

# ------------------------------------------------------------------ down-sampling n
down = pd.read_csv(os.path.join(S12, "S12_downsampling.csv"))
ns = (down.groupby(["design", "cells_per_sample"])
      .agg(n_tumor=("n_tumor", "first"), n_normal=("n_normal", "first"),
           n_repeats=("n_repeats", "first")).reset_index())
ns["samples_excluded_from_60_tumors"] = ns.n_tumor.map(lambda v: 60 - int(v))
ns["same_sample_set_across_replicates"] = True
ns.to_csv(os.path.join(OUT, "S13_downsampling_sample_counts.csv"), index=False,
          encoding="utf-8-sig")

# ------------------------------------------------------------------ methods / README
methods = [
    ("Primary design", "56 unmatched tumours vs 4 adjacent non-malignant libraries (60 libraries)."),
    ("Primary enumeration", "All C(60,4) = 487,635 label assignments; 4 of 60 libraries labelled adjacent."),
    ("LOO design", "56 unmatched tumours vs 3 adjacent libraries (59 libraries), for each of 4 deletions."),
    ("LOO enumeration", "All C(59,3) = 32,509 label assignments per deletion."),
    ("Statistic", "Mann-Whitney U computed on raw pairwise comparisons: 1 for a>b, 0.5 for a tie, 0 otherwise. Integer 2U used internally."),
    ("Tie handling", "Ties are retained at 0.5; no untied-rank distribution is assumed. Tied observations (mainly zeros) are therefore handled exactly."),
    ("Two-sided P (primary definition)",
     "Symmetric about the null centre n1*n2/2: P = #{|U - centre| >= |U_obs - centre|} / #assignments. The observed assignment is included in the enumeration."),
    ("Two-sided P (secondary definition)",
     "2 x min(P(U <= U_obs), P(U >= U_obs)), capped at 1; reported for comparison."),
    ("BH families",
     "Primary: 9 states. LOO: jointly across the 16 prespecified comparisons (4 states x 4 deletions) and separately within each deletion across 4 states; the additional inflamed-venule LOO uses its own 4-test family."),
    ("Comparators",
     "SciPy 1.18.1 mannwhitneyu(method='exact'), which does not correct ties, and method='asymptotic', which does; both are reported with BH q values."),
    ("Down-sampling", "Samples with fewer than the target cell count were excluded; the same eligible sample set was used in every replicate."),
    ("Software", "Python 3.13.14, numpy 2.5.4, scipy 1.18.1, statsmodels 0.15.0; script 13_refcheck/86_perm_v3_6.py; seed-free (complete enumeration)."),
]
pd.DataFrame(methods, columns=["Item", "Description"]).to_csv(
    os.path.join(OUT, "S13_methods.csv"), index=False, encoding="utf-8-sig")

summary = dict(
    primary_assignments=int(ncomb), pooled_assignments=int(n64),
    loo_assignments_per_deletion=32509,
    pooled_capillary_q_perm=float(pooled.loc[pooled.state == "Capillary", "q_perm"].iloc[0]),
    pooled_capillary_q_exact=float(pooled.loc[pooled.state == "Capillary", "q_exact"].iloc[0]),
    primary_states=int(len(primary)),
    all_three_concordant=bool(primary.concordant_all_three.all()),
    n_sig_perm=int(primary.sig_perm.sum()),
    n_sig_exact=int(primary.sig_exact.sum()),
    n_sig_asymptotic=int(primary.sig_asymptotic.sum()),
    min_p_perm=float(primary.p_perm_sym.min()),
    max_q_perm=float(primary.q_perm.max()),
    loo16_directions_preserved=int(loo16.direction_preserved.sum()),
    loo16_joint_q_min=float(loo16.q_joint_16.min()),
    loo16_joint_q_max=float(loo16.q_joint_16.max()),
    inflamed_loo_q_max=float(inf.q_within_4.max()),
    downsampling=ns.to_dict(orient="records"),
)
with io.open(os.path.join(OUT, "S13_summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

print("[OK] primary 56v4 : %d assignments, %d states" % (ncomb, len(primary)))
print("[OK] LOO 56v3     : 32509 assignments x 4 deletions x 5 states")
print("[OK] concordant all three methods:", summary["all_three_concordant"])
print("[OK] significant: perm=%d exact=%d asymptotic=%d"
      % (summary["n_sig_perm"], summary["n_sig_exact"], summary["n_sig_asymptotic"]))
print("[OK] LOO16 directions preserved: %d/16, joint q %.5f-%.5f"
      % (summary["loo16_directions_preserved"],
         summary["loo16_joint_q_min"], summary["loo16_joint_q_max"]))
print("[OK] inflamed-venule LOO max q: %.6f" % summary["inflamed_loo_q_max"])
print("--- primary permutation results ---")
print(primary[["state", "p_perm_sym", "q_perm", "p_exact_scipy", "p_asymptotic_tiecorrected",
               "sig_perm", "concordant_all_three"]].to_string(index=False))
