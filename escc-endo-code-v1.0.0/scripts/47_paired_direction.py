# -*- coding: utf-8 -*-
"""Patient-level paired direction test across the two independent ESCC cohorts."""
import os
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

BASE = r"D:/龙虾文档/2026-11-02-09-16-41/ESCC_endo"
OUT = os.path.join(BASE, "14_external_validation", "out")
d = pd.read_csv(os.path.join(OUT, "TableS9_state_proportions_by_sample.csv"),
                encoding="utf-8-sig")
d = d[d.method == "centroid_rank"].copy()
STATES = [("Venous", "N>T"), ("Tip_angiogenic", "N<T"),
          ("IFNgamma_EC", "N<T")]

rows = []
for c in ["GSE196756", "GSE203067"]:
    x = d[d.cohort == c].copy()
    x["pair"] = (x["sample"].str.extract(r"_(\d)[TN]")[0] if c == "GSE196756"
                 else x["sample"].str.extract(r"[NT](\d)")[0])
    for k, g in x.groupby("pair"):
        gn = g[g.group == "Normal"]
        gt = g[g.group == "Tumor"]
        if len(gn) != 1 or len(gt) != 1:
            continue
        for s, expect in STATES:
            n, t = float(gn[s].iloc[0]), float(gt[s].iloc[0])
            ok = (n > t) if expect == "N>T" else (n < t)
            rows.append(dict(cohort=c, pair=k, state=s, expected=expect,
                             frac_normal=round(n, 4), frac_tumor=round(t, 4),
                             absolute_diff=round(t - n, 4),
                             n_cells_normal=int(gn.n_cells.iloc[0]),
                             n_cells_tumor=int(gt.n_cells.iloc[0]),
                             direction_ok=bool(ok)))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "TableS9_paired_direction.csv"),
                          index=False, encoding="utf-8-sig")

r = pd.DataFrame(rows)
out = []
for s, expect in STATES:
    g = r[r.state == s]
    k, n = int(g.direction_ok.sum()), len(g)
    p = binomtest(k, n, 0.5, alternative="greater").pvalue
    out.append(dict(state=s, expected=expect, n_pairs=n, pairs_concordant=k,
                    p_one_sided=p, mean_abs_diff=round(float(g.absolute_diff.mean()), 4)))
q = multipletests([x["p_one_sided"] for x in out], method="fdr_bh")[1]
for x, qq in zip(out, q):
    x["fdr_bh"] = float(qq)
pd.DataFrame(out).to_csv(os.path.join(OUT, "TableS9_paired_sign_test.csv"),
                         index=False, encoding="utf-8-sig")
print(pd.DataFrame(out).to_string(index=False))
