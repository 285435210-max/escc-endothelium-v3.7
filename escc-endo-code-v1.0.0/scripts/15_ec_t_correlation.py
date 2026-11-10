# -*- coding: utf-8 -*-
"""T 细胞亚群比例 vs 内皮状态比例：样本水平关联分析（60 例肿瘤）"""
import gzip, os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

SRC = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo\data_manifest\GSE160269_UMI_matrix_Tcell.txt.gz"
OUT = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo\11_correlation"
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "tcell_subset_by_sample.csv")

PANEL = {"CD8A", "CD8B", "CD4", "FOXP3", "IL2RA", "PDCD1", "HAVCR2", "LAG3",
         "TIGIT", "TOX", "GNLY", "NKG7", "KLRD1", "CD3D", "CD3E"}
SUBSETS = ["CD8 T", "CD8 exhausted", "CD4 T", "Treg", "NK-like"]
STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
          "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"]


def scan_tcells():
    with gzip.open(SRC, "rt", encoding="utf-8", errors="ignore") as f:
        cells = [c.strip() for c in f.readline().rstrip("\n").split(" ") if c.strip()]
        n = len(cells)
        print(f"T cells: {n}", flush=True)
        panel = {}
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            i = line.find(" ")
            g = line[:i]
            if g not in PANEL:
                continue
            v = np.fromstring(line[i + 1:], sep=" ", dtype=np.float32)
            if v.size != n:
                v = np.pad(v, (0, max(0, n - v.size)))[:n]
            panel[g] = v > 0

    def has(g):
        return panel.get(g, np.zeros(n, dtype=bool))

    cd8 = has("CD8A") | has("CD8B")
    cd4 = has("CD4") & ~cd8
    treg = cd4 & has("FOXP3") & has("IL2RA")
    ex = sum(has(g).astype(int) for g in ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"]) >= 2
    nk = (~has("CD3D") & ~has("CD3E")) & (has("GNLY") | has("NKG7") | has("KLRD1"))

    calls = {"CD8 T": cd8 & ~ex, "CD8 exhausted": cd8 & ex,
             "CD4 T": cd4 & ~treg & ~ex, "Treg": treg, "NK-like": nk}
    sample = np.array([c.split("-")[0] for c in cells])
    df = pd.DataFrame({"sample": sample})
    for k, v in calls.items():
        df[k] = v.astype(int)
    agg = df.groupby("sample").sum()
    agg["total"] = agg[SUBSETS].sum(axis=1)
    agg.to_csv(CACHE, encoding="utf-8-sig")
    print("samples:", agg.shape[0], flush=True)
    return agg


if os.path.exists(CACHE):
    agg = pd.read_csv(CACHE, index_col=0)
    print("cached T-subset counts:", agg.shape, flush=True)
else:
    agg = scan_tcells()

# ---------- 内皮状态比例（未检出记 0） ----------
ec = pd.read_csv(r"D:\escc_sc\data\celltype_proportion_by_sample.csv")
ec = ec[~ec["celltype"].str.startswith("Contam") & (ec["type"] == "Tumor")].copy()
ec["patient"] = ec["sample"].str.replace(r"T$", "", regex=True)
grid = pd.MultiIndex.from_product([sorted(ec["patient"].unique()), STATES],
                                  names=["patient", "state"]).to_frame(index=False)
ecf = grid.merge(ec, left_on=["patient", "state"], right_on=["patient", "celltype"], how="left")
ecf["frac"] = ecf["frac"].fillna(0)
ecn = ecf.groupby("patient")["n"].sum().rename("ec_n")
wide = ecf.pivot(index="patient", columns="state", values="frac")

# ---------- 合并 ----------
agg.index = [s[:-1] if s.endswith("T") else s for s in agg.index]
t_n = agg["total"].rename("t_n")
tfrac = agg[SUBSETS].div(agg["total"], axis=0)
tfrac.index.name = "patient"
m = wide.join(tfrac, how="inner").join(ecn, how="inner").join(t_n, how="inner")
print("matched patients:", len(m), flush=True)
m.to_csv(os.path.join(OUT, "merged_proportions.csv"), encoding="utf-8-sig")


def partial_spearman(x, y, z):
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    rxy = np.corrcoef(rx, ry)[0, 1]
    rxz = np.corrcoef(rx, rz)[0, 1]
    ryz = np.corrcoef(ry, rz)[0, 1]
    r = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    k = len(x)
    t = r * np.sqrt((k - 3) / (1 - r ** 2))
    return r, 2 * (1 - stats.t.cdf(abs(t), k - 3))


def corr_table(m, min_ec=0, min_t=0):
    d = m[(m["ec_n"] >= min_ec) & (m["t_n"] >= min_t)]
    rows = []
    for s in STATES:
        for t in SUBSETS:
            if d[s].std() == 0 or d[t].std() == 0:
                continue
            r, p = stats.spearmanr(d[s], d[t])
            rp, pp = partial_spearman(d[s].values, d[t].values, np.log10(d["ec_n"].values + 1))
            rows.append(dict(ec_state=s, t_subset=t, n=len(d), rho=r, p=p,
                             rho_partial=rp, p_partial=pp))
    out = pd.DataFrame(rows)
    out["fdr"] = multipletests(out["p"], method="fdr_bh")[1]
    out["fdr_partial"] = multipletests(out["p_partial"], method="fdr_bh")[1]
    return out


res = corr_table(m)
res.to_csv(os.path.join(OUT, "EC_vs_T_correlation.csv"), index=False, encoding="utf-8-sig")
sens = corr_table(m, min_ec=30, min_t=50)
sens.to_csv(os.path.join(OUT, "EC_vs_T_correlation_sensitivity.csv"), index=False, encoding="utf-8-sig")

print("\n=== raw (n=%d) top 12 ===" % len(m))
print(res.sort_values("p").head(12).round(3).to_string(index=False))
print("\n=== sensitivity (EC>=30, T>=50, n=%d) top 10 ===" % int(sens["n"].iloc[0]))
print(sens.sort_values("p").head(10).round(3).to_string(index=False))

# ---------- 热图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mat = res.pivot(index="ec_state", columns="t_subset", values="rho").loc[STATES, SUBSETS]
fdr = res.pivot(index="ec_state", columns="t_subset", values="fdr").loc[STATES, SUBSETS]
matp = res.pivot(index="ec_state", columns="t_subset", values="rho_partial").loc[STATES, SUBSETS]
fdrp = res.pivot(index="ec_state", columns="t_subset", values="fdr_partial").loc[STATES, SUBSETS]

fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))
for ax, M, F, ttl in [(axes[0], mat, fdr, "Spearman rho"),
                      (axes[1], matp, fdrp, "Partial rho (controlling EC cell number)")]:
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(M.columns, rotation=35, ha="right", fontsize=7)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(M.index, fontsize=7)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            star = "*" if F.values[i, j] < 0.05 else ""
            ax.text(j, i, f"{M.values[i, j]:.2f}{star}", ha="center", va="center", fontsize=5.5)
    ax.set_title(ttl, fontsize=8)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
fig.colorbar(im, ax=list(axes), shrink=0.8, label="correlation")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "FigS2_EC_T_correlation.png"), dpi=600, bbox_inches="tight")
print("\nsaved:", OUT)
