# -*- coding: utf-8 -*-
"""
B2 + B3 修订（审稿意见 M2 / M3 / M8）。

问题：
  1) 旧 frac 的分母 n_total 含后来剔除的 Contam 细胞，9 个状态之和平均仅 0.93，
     与稿件 "9,531 endothelial cells" / "fraction of endothelial cells" 不符。
  2) 成分闭合可能造成静脉-tip 负相关伪影；从未做 CLR 或零模型检验。
  3) 癌旁 4 例 / 856 细胞 vs 肿瘤 60 例 / 8,675 细胞，检出限偏倚；从未做等细胞数下采样。
  4) 6-9 次 Wilcoxon 未做 BH 校正；4 对配对样本弃用配对设计；新口径 P 值未落地脚本。
  5) IFNgamma-EC 与 ex_of_cd8 的相关只报名义 P，FDR>=0.11 引用的是另一套 9x5 结局。

本脚本统一口径：frac = n / sum(9 个保留状态的 n)，每个样本和为 1；未检出记 0。
输出 08_scRNA/endothelial/revised/ 与 11_correlation/revised/，并生成数字变更对照表。
"""
import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SC = r"D:\escc_sc"
OUT = os.path.join(ROOT, "08_scRNA", "endothelial", "revised")
COUT = os.path.join(ROOT, "11_correlation", "revised")
os.makedirs(OUT, exist_ok=True)
os.makedirs(COUT, exist_ok=True)
SEED = 20261105
STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
          "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"]

# ---------------------------------------------------------------- 1. 统一口径：闭合成分
raw = pd.read_csv(os.path.join(ROOT, "08_scRNA", "endothelial",
                               "celltype_proportion_by_sample.csv"))
raw = raw[~raw["celltype"].astype(str).str.startswith("Contam")].copy()
raw["celltype"] = raw["celltype"].astype(str)
meta = raw.groupby("sample")["type"].first()
tot9 = raw.groupby("sample")["n"].sum().rename("n_total9")
grid = pd.MultiIndex.from_product([sorted(raw["sample"].unique()), STATES],
                                  names=["sample", "celltype"]).to_frame(index=False)
grid["type"] = grid["sample"].map(meta)
g = grid.merge(raw[["sample", "type", "celltype", "n", "n_total", "frac"]],
               on=["sample", "type", "celltype"], how="left")
g["n"] = g["n"].fillna(0).astype(int)
g["n_total"] = g["sample"].map(tot9).astype(int)
g["frac"] = np.where(g["n_total"] > 0, g["n"] / g["n_total"], 0.0)
closed = g[["sample", "type", "celltype", "n", "n_total", "frac"]].copy()
assert np.allclose(closed.groupby("sample")["frac"].sum(), 1.0), "闭合失败"
closed.to_csv(os.path.join(OUT, "celltype_proportion_closed.csv"), index=False)
# 下游 R/Python 直接读这个文件（列名与原文件兼容）
closed.to_csv(os.path.join(SC, "data", "celltype_proportion_closed.csv"), index=False)
print(f"[closed] samples={closed['sample'].nunique()} sum(frac)=1.0 "
      f"tumor={sum(meta=='Tumor')} normal={sum(meta=='Normal')}", flush=True)

wide = closed.pivot(index="sample", columns="celltype", values="frac")[STATES]
types = closed.groupby("sample")["type"].first().loc[wide.index]
n9 = closed.groupby("sample")["n_total"].first().loc[wide.index]
cnt = closed.pivot(index="sample", columns="celltype", values="n")[STATES]

# ---------------------------------------------------------------- 2. 组间比较 + BH + 配对
rows = []
for s in STATES:
    t = wide.loc[types == "Tumor", s].to_numpy()
    nn = wide.loc[types == "Normal", s].to_numpy()
    u, p = stats.mannwhitneyu(t, nn, alternative="two-sided")
    # 4 对配对样本（P126/P127/P128/P130）
    pairs, diffs = [], []
    for pat in ["P126", "P127", "P128", "P130"]:
        st, sn = f"{pat}T", f"{pat}N"
        if st in wide.index and sn in wide.index:
            pairs.append((pat, float(wide.loc[st, s]), float(wide.loc[sn, s])))
            diffs.append(float(wide.loc[st, s]) - float(wide.loc[sn, s]))
    sign_p = stats.binomtest(sum(d > 0 for d in diffs), len(diffs), 0.5).pvalue if diffs else np.nan
    try:
        wr_p = stats.wilcoxon([d for d in diffs if d != 0]).pvalue if sum(d != 0 for d in diffs) >= 2 else np.nan
    except Exception:
        wr_p = np.nan
    rows.append(dict(state=s, med_tumor=float(np.median(t)), med_normal=float(np.median(nn)),
                     n_tumor=len(t), n_normal=len(nn), p_raw=p,
                     n_pairs=len(diffs), n_concordant=sum(d > 0 for d in diffs),
                     p_sign=sign_p, p_signedrank=wr_p))
comp = pd.DataFrame(rows)
comp["p_bh"] = multipletests(comp["p_raw"], method="fdr_bh")[1]
comp["sig_bh"] = comp["p_bh"] < 0.05
comp.to_csv(os.path.join(OUT, "composition_group_tests.csv"), index=False)
print("\n[group tests, closed composition + BH]\n" +
      comp[["state", "med_tumor", "med_normal", "p_raw", "p_bh", "sig_bh"]]
      .round(4).to_string(index=False), flush=True)

from math import comb
min_p_unpaired = 2 / comb(64, 4)
min_p_paired = 2 * 0.5 ** 4
print(f"[power] 4 vs 60 rank-sum 理论最小双侧 P = {min_p_unpaired:.3g}; "
      f"n=4 符号检验最小双侧 P = {min_p_paired}", flush=True)

# ---------------------------------------------------------------- 3. CLR + 成分零模型
T = wide.loc[types == "Tumor"]
rho_raw = stats.spearmanr(T["Venous"], T["Tip_angiogenic"])
print(f"[corr] raw closed  Venous~Tip rho={rho_raw.statistic:.4f} P={rho_raw.pvalue:.3g}", flush=True)

# 零值用 0.5 计数替换后再闭合（Bayesian-multiplicative replacement）
cT = cnt.loc[T.index].to_numpy(dtype=float)
pseudo = np.where(cT == 0, 0.5, 0.0)
cb = cT + pseudo
pb = cb / cb.sum(axis=1, keepdims=True)
clr = np.log(pb) - np.log(pb).mean(axis=1, keepdims=True)
rho_clr = stats.spearmanr(clr[:, STATES.index("Venous")], clr[:, STATES.index("Tip_angiogenic")])
print(f"[corr] CLR          Venous~Tip rho={rho_clr.statistic:.4f} P={rho_clr.pvalue:.3g}", flush=True)

# Dirichlet-multinomial 零模型：保留样本真实细胞数 + 成分闭合 + 抽样噪声
rng = np.random.default_rng(SEED)
nT = n9.loc[T.index].to_numpy().astype(int)
p_glob = cT.sum(axis=0) / cT.sum()
obs_var = T["Venous"].var(ddof=1)


def null_rhos_fast(concentration, nrep=4000):
    a = p_glob * concentration
    vi, ti = STATES.index("Venous"), STATES.index("Tip_angiogenic")
    out = np.empty(nrep)
    for b in range(nrep):
        prob = rng.dirichlet(a, size=len(nT))
        v = np.empty(len(nT)); t = np.empty(len(nT))
        for i, ni in enumerate(nT):
            k = rng.multinomial(int(ni), prob[i])
            v[i], t[i] = k[vi] / ni, k[ti] / ni
        out[b] = stats.spearmanr(v, t).statistic
    return out


# 用矩估计选 concentration：使零模型下静脉比例的样本间方差接近观测值
def var_at(c, nrep=200):
    vi = STATES.index("Venous")
    a = p_glob * c
    vals = []
    for _ in range(nrep):
        prob = rng.dirichlet(a, size=len(nT))
        v = np.array([rng.multinomial(int(ni), prob[i])[vi] / ni for i, ni in enumerate(nT)])
        vals.append(v.var(ddof=1))
    return float(np.mean(vals))


grid_c = [5, 20, 50, 100, 200, 400]
var_grid = {c: var_at(c, 60) for c in grid_c}
c_fit = min(var_grid, key=lambda c: abs(var_grid[c] - obs_var))
print(f"[null] observed Var(venous)={obs_var:.4f}; concentration grid {var_grid}; chosen c={c_fit}", flush=True)

null_rows = []
for label, c in [("pure sampling noise", None), (f"Dirichlet-multinomial c={c_fit}", c_fit)]:
    if c is None:
        v, t = [], []
        for _ in range(4000):
            k = np.array([rng.multinomial(int(ni), p_glob) for ni in nT])
            v.append(k[:, STATES.index("Venous")] / nT)
            t.append(k[:, STATES.index("Tip_angiogenic")] / nT)
        rr = np.array([stats.spearmanr(v[b], t[b]).statistic for b in range(4000)])
    else:
        rr = null_rhos_fast(c, 4000)
    emp = (1 + (rr <= rho_raw.statistic).sum()) / (len(rr) + 1)
    null_rows.append(dict(null_model=label, n_rep=len(rr),
                          rho_median=float(np.median(rr)),
                          rho_q025=float(np.percentile(rr, 2.5)),
                          rho_q975=float(np.percentile(rr, 97.5)),
                          empirical_P_one_sided_negative=float(emp)))
    print(f"[null] {label}: median rho={np.median(rr):.3f} "
          f"95%=[{np.percentile(rr,2.5):.3f},{np.percentile(rr,97.5):.3f}] "
          f"P(rho<=observed)={emp:.4f}", flush=True)
null_df = pd.DataFrame(null_rows)
null_df.to_csv(os.path.join(OUT, "compositional_null_models.csv"), index=False)
pd.DataFrame({"model": ["raw closed", "CLR"],
              "rho": [rho_raw.statistic, rho_clr.statistic],
              "p": [rho_raw.pvalue, rho_clr.pvalue]}).to_csv(
    os.path.join(OUT, "venous_tip_correlation.csv"), index=False)

# ---------------------------------------------------------------- 4. 等细胞数下采样敏感性
SUBSETS = ["CD8 T", "CD8 exhausted", "CD4 T", "Treg", "NK-like"]
down_rows = []
for N in [20, 30]:
    keep = [s for s in wide.index if n9[s] >= N and types[s] in ("Tumor", "Normal")]
    kt = [s for s in keep if types[s] == "Tumor"]
    kn = [s for s in keep if types[s] == "Normal"]
    if len(kn) < 2:
        print(f"[downsample] N={N}: 癌旁样本不足，跳过")
        continue
    R = 200
    pv = np.zeros((R, len(STATES)))
    for b in range(R):
        for j in range(len(STATES)):
            vals = []
            for s in keep:
                counts = cnt.loc[s, STATES].to_numpy(dtype=int)
                if counts.sum() < N:
                    vals.append(np.nan); continue
                draw = rng.multivariate_hypergeometric(counts, N) if hasattr(rng, "multivariate_hypergeometric") \
                    else rng.multinomial(N, counts / counts.sum())
                vals.append(draw[j] / N)
            v = np.array(vals, dtype=float)
            a = v[[keep.index(s) for s in kt]]
            c = v[[keep.index(s) for s in kn]]
            pv[b, j] = stats.mannwhitneyu(a, c, alternative="two-sided").pvalue if len(c) >= 2 else np.nan
    fdr = np.apply_along_axis(lambda x: multipletests(np.nan_to_num(x, nan=1.0), method="fdr_bh")[1], 1, pv)
    for j, s in enumerate(STATES):
        down_rows.append(dict(target_cells=N, n_tumor=len(kt), n_normal=len(kn), state=s,
                              median_p=float(np.nanmedian(pv[:, j])),
                              frac_bh_significant=float(np.mean(fdr[:, j] < 0.05)),
                              median_bh=float(np.nanmedian(fdr[:, j]))))
down = pd.DataFrame(down_rows)
down.to_csv(os.path.join(OUT, "downsampling_sensitivity.csv"), index=False)
print("\n[downsample] BH 显著比例（200 次重抽样）\n" +
      down.pivot(index="state", columns="target_cells", values="frac_bh_significant")
      .round(2).to_string(), flush=True)

# ---------------------------------------------------------------- 5. B3：相关性 + BH
tcache = pd.read_csv(os.path.join(ROOT, "11_correlation", "tcell_subset_by_sample.csv"), index_col=0)
tcache.index = [str(i)[:-1] if str(i).endswith("T") else str(i) for i in tcache.index]
tfrac = tcache[SUBSETS].div(tcache["total"], axis=0)
cd8_all_n = tcache["CD8 T"] + tcache["CD8 exhausted"]
tfrac["ex_of_cd8"] = tcache["CD8 exhausted"] / cd8_all_n.replace(0, np.nan)
tfrac["t_n"] = tcache["total"]

ec = wide.loc[types == "Tumor"].copy()
ec.index = [str(i)[:-1] if str(i).endswith("T") else str(i) for i in ec.index]
ec["ec_n"] = [n9[s] for s in wide.loc[types == "Tumor"].index]
m = ec.join(tfrac, how="inner").dropna(subset=["ex_of_cd8"])
m.to_csv(os.path.join(COUT, "merged_proportions_closed.csv"), encoding="utf-8-sig")
# 同步给绘图脚本使用的主文件（旧版另存备份）
main_merged = os.path.join(ROOT, "11_correlation", "merged_proportions.csv")
if os.path.exists(main_merged) and not os.path.exists(main_merged + ".legacy"):
    os.rename(main_merged, main_merged + ".legacy")
m.to_csv(main_merged, encoding="utf-8-sig")

sub = m[(m.ec_n >= 30) & (m.t_n >= 50)]
print(f"\n[B3] all tumors n={len(m)}, adequately sampled n={len(sub)} (post hoc)", flush=True)


def bh_table(d, outcome, label):
    rr = []
    for s in STATES:
        r, p = stats.spearmanr(d[s], d[outcome])
        rr.append(dict(state=s, outcome=outcome, n=len(d), rho=r, p=p))
    o = pd.DataFrame(rr)
    o["fdr_bh"] = multipletests(o["p"], method="fdr_bh")[1]
    o["analysis"] = label
    o.to_csv(os.path.join(COUT, f"{label}.csv"), index=False)
    return o


all60 = bh_table(m, "ex_of_cd8", "ExCD8_all60_BH")
sub47 = bh_table(sub, "ex_of_cd8", "ExCD8_sub47_posthoc_BH")
print("\n[B3] ex_of_cd8, n=60, BH across 9 states\n" +
      all60[["state", "rho", "p", "fdr_bh"]].round(4).to_string(index=False), flush=True)
print("\n[B3] ex_of_cd8, n=47 post hoc\n" +
      sub47[["state", "rho", "p", "fdr_bh"]].round(4).to_string(index=False), flush=True)

# 9x5 原始 T 亚群结局（替代旧 EC_vs_T_correlation.csv 的口径）
all45 = []
for s in STATES:
    for t in SUBSETS:
        r, p = stats.spearmanr(m[s], m[t])
        all45.append(dict(state=s, t_subset=t, n=len(m), rho=r, p=p))
a45 = pd.DataFrame(all45)
a45["fdr_bh"] = multipletests(a45["p"], method="fdr_bh")[1]
a45.to_csv(os.path.join(COUT, "EC_vs_T_9x5_BH.csv"), index=False)
print(f"\n[B3] 9x5 最小 FDR = {a45.fdr_bh.min():.4f}; "
      f"最小名义 P = {a45.p.min():.4f}", flush=True)

# ---------------------------------------------------------------- 6. 关键数字与变更表
key = {
    "composition_definition": "frac = n / sum of the nine retained states per sample (closed, sums to 1)",
    "n_samples": {"tumor": int((types == "Tumor").sum()), "normal": int((types == "Normal").sum())},
    "group_tests": comp.set_index("state")[["med_tumor", "med_normal", "p_raw", "p_bh"]].to_dict("index"),
    "min_p_unpaired_4v60": min_p_unpaired,
    "min_p_paired_n4": min_p_paired,
    "venous_tip": {"raw_rho": float(rho_raw.statistic), "raw_p": float(rho_raw.pvalue),
                   "clr_rho": float(rho_clr.statistic), "clr_p": float(rho_clr.pvalue),
                   "null_models": null_rows},
    "downsampling": down.to_dict("records"),
    "ex_of_cd8": {"n60_rho_IFNg": float(all60.set_index("state").loc["IFNgamma_EC", "rho"]),
                  "n60_p_IFNg": float(all60.set_index("state").loc["IFNgamma_EC", "p"]),
                  "n60_fdr_IFNg": float(all60.set_index("state").loc["IFNgamma_EC", "fdr_bh"]),
                  "n60_min_fdr": float(all60.fdr_bh.min()),
                  "n47_rho_IFNg": float(sub47.set_index("state").loc["IFNgamma_EC", "rho"]),
                  "n47_p_IFNg": float(sub47.set_index("state").loc["IFNgamma_EC", "p"]),
                  "n47_fdr_IFNg": float(sub47.set_index("state").loc["IFNgamma_EC", "fdr_bh"]),
                  "n47_min_fdr": float(sub47.fdr_bh.min()),
                  "n": len(m), "n_sub": len(sub)},
    "min_fdr_9x5": float(a45.fdr_bh.min()),
    "seed": SEED,
}
with open(os.path.join(OUT, "key_numbers.json"), "w", encoding="utf-8") as f:
    json.dump(key, f, ensure_ascii=False, indent=2)

lines = ["# B2/B3 重算：稿件数字变更对照表", "",
         "统一口径：frac = 该样本 9 个保留内皮状态的细胞数占比（和为 1）；未检出记 0。",
         f"样本：肿瘤 {int((types=='Tumor').sum())} 例，癌旁 {int((types=='Normal').sum())} 例。", "",
         "## 组间比较（中位数，Wilcoxon rank-sum，BH 校正 9 次）", "",
         "| 状态 | 肿瘤 | 癌旁 | 名义 P | BH-FDR | 校正后显著 |", "|---|---|---|---|---|---|"]
for _, r in comp.iterrows():
    lines.append(f"| {r.state} | {r.med_tumor*100:.1f}% | {r.med_normal*100:.1f}% | "
                 f"{r.p_raw:.4g} | {r.p_bh:.4g} | {'是' if r.sig_bh else '否'} |")
lines += ["", f"配对：4 对（P126/P127/P128/P130）；n=4 时符号检验最小双侧 P = {min_p_paired}（**配对设计无功效**）。",
          f"4 vs 60 非配对 rank-sum 理论最小双侧 P = {min_p_unpaired:.3g}。", "",
          "## 静脉–tip 相关与成分闭合", "",
          f"- 闭合比例原始 Spearman rho = {rho_raw.statistic:.3f}（P = {rho_raw.pvalue:.3g}）",
          f"- CLR 变换后 rho = {rho_clr.statistic:.3f}（P = {rho_clr.pvalue:.3g}）"]
for r in null_rows:
    lines.append(f"- 零模型「{r['null_model']}」：rho 中位 {r['rho_median']:.3f}，"
                 f"95% 区间 [{r['rho_q025']:.3f}, {r['rho_q975']:.3f}]，"
                 f"P(rho<=观测) = {r['empirical_P_one_sided_negative']:.4f}")
lines += ["", "## ex_of_cd8 相关（BH 校正 9 个状态）", "",
          f"- 全 60 例：IFNgamma-EC rho = {key['ex_of_cd8']['n60_rho_IFNg']:.3f}，"
          f"P = {key['ex_of_cd8']['n60_p_IFNg']:.4g}，FDR = {key['ex_of_cd8']['n60_fdr_IFNg']:.3f}；"
          f"9 状态最小 FDR = {key['ex_of_cd8']['n60_min_fdr']:.3f}",
          f"- 47 例（**事后**充分采样亚组）：rho = {key['ex_of_cd8']['n47_rho_IFNg']:.3f}，"
          f"P = {key['ex_of_cd8']['n47_p_IFNg']:.4g}，FDR = {key['ex_of_cd8']['n47_fdr_IFNg']:.3f}；"
          f"最小 FDR = {key['ex_of_cd8']['n47_min_fdr']:.3f}",
          f"- 9x5 原始 T 亚群结局最小 FDR = {a45.fdr_bh.min():.3f}（稿件原引 FDR>=0.11 属此套分析，且已变更）", "",
          "## 等细胞数下采样敏感性", "",
          "见 `downsampling_sensitivity.csv`（每状态 200 次多元超几何重抽样后 BH 显著比例）。", ""]
with open(os.path.join(OUT, "数字变更对照表.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("\n[done] B2+B3 ->", OUT, flush=True)
