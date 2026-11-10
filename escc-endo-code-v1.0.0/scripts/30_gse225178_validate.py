# -*- coding: utf-8 -*-
"""
GSE225178 外部验证（Python 手动解析版）：
tip 签名 vs ESCC 新辅助化疗组织学反应。
不依赖 GEOquery，直接解析 series matrix 与 GPL SOFT。
"""
import gzip, re, os
import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\escc_sc\gse225178"
SM = os.path.join(BASE, "GSE225178_series_matrix.txt.gz")
GPL = os.path.join(BASE, "GPL21185.soft.gz")
OUT = os.path.join(BASE, "out"); os.makedirs(OUT, exist_ok=True)

# ---------- 1. 解析 series matrix 表头（临床）与表达矩阵 ----------
def read_series_matrix(path):
    header = []
    table_rows = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_rows.append(line)
            else:
                header.append(line)
    return header, table_rows

header, rows = read_series_matrix(SM)

# 表头字段
def hdr_vals(key):
    for l in header:
        if l.startswith("!" + key):
            return [x.strip('"') for x in l.split("\t")[1:]]
    return []

gsm = hdr_vals("Sample_geo_accession")
source = hdr_vals("Sample_source_name_ch1")
# 临床特征
chars = {}
for l in header:
    if l.startswith("!Sample_characteristics_ch1"):
        vals = [x.strip('"') for x in l.split("\t")[1:]]
        if vals:
            k = vals[0].split(":")[0].strip()
            chars[k] = [v.split(":", 1)[-1].strip() for v in vals]

print("样本数:", len(gsm))
print("特征字段:", list(chars.keys()))

# 找反应分级字段
resp_key = None
for k in chars:
    if re.search(r"response|grade", k, re.I):
        resp_key = k
print("反应字段:", resp_key)
print("反应取值分布:", pd.Series(chars[resp_key]).value_counts().to_dict() if resp_key else None)

# 找 tumor/non-tumor 字段
tissue_key = None
for k in chars:
    if re.search(r"tumor|tissue|non-tumor", k, re.I):
        tissue_key = k
print("组织字段:", tissue_key, "| 分布:", pd.Series(chars[tissue_key]).value_counts().to_dict() if tissue_key else None)

# ---------- 2. 表达矩阵 ----------
# 第一行是 ID_REF + 样本列
colnames = [x.strip('"') for x in rows[0].split("\t")]
probes = []
mat = []
for r in rows[1:]:
    parts = r.split("\t")
    if len(parts) != len(colnames):
        continue
    probes.append(parts[0].strip('"'))
    mat.append([float(x.strip('"')) for x in parts[1:]])
expr = pd.DataFrame(mat, index=probes, columns=colnames[1:])
print("表达矩阵:", expr.shape, "| 样本列数:", expr.shape[1])

# 样本与临床对齐
meta = pd.DataFrame({
    "gsm": gsm,
    "source": source,
})
for k in chars:
    meta[k] = chars[k]
meta.index = meta["gsm"]
meta = meta.loc[expr.columns]
print("meta:", meta.shape)

# ---------- 3. GPL 注释：探针 → GENE_SYMBOL ----------
def read_gpl_soft(path):
    mapping = {}
    in_table = False
    header_cols = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                parts = line.split("\t")
                if not header_cols:
                    header_cols = parts
                    continue
                if len(parts) >= len(header_cols):
                    d = dict(zip(header_cols, parts))
                    mapping[d.get("ID", "")] = d
    return mapping

gpl = read_gpl_soft(GPL)
print("\nGPL 探针数:", len(gpl))
# 示例
sample_id = list(gpl.keys())[0]
print("样例:", sample_id, "->", gpl[sample_id].get("GENE_SYMBOL", "") if sample_id in gpl else "")

probe2gene = {}
for pid, d in gpl.items():
    gs = d.get("GENE_SYMBOL", "").strip()
    if gs and not re.match(r"^(DarkCorner|BrightCorner|ETG|E1A|r60|_|Corner|Control)", gs):
        probe2gene[pid] = gs
print("有效探针-基因映射:", len(probe2gene))

# ---------- 4. tip 签名基因 ----------
tip = pd.read_csv(r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo\09_de\Tip_signature_shared.csv")["gene"].unique().tolist()
core = ["ESM1", "APLN", "NID2", "PLVAP", "COL4A1", "COL4A2", "LAMC1"]
breach = ["COL4A1", "COL4A2", "LAMC1", "COL5A1", "COL12A1", "MMP2", "MMP14"]

gene2probes = {}
for pid, gs in probe2gene.items():
    gene2probes.setdefault(gs, []).append(pid)

def score_signature(genes):
    """基因集分数：每个基因取平均表达最高的探针，样本内 z-score 后取均值。"""
    sel = []
    for g in genes:
        if g in gene2probes:
            ps = [p for p in gene2probes[g] if p in expr.index]
            if ps:
                # 取表达最高（平均）的探针
                m = expr.loc[ps].mean(axis=1)
                sel.append((g, m.idxmax()))
    if not sel:
        return None
    sub = expr.loc[[p for _, p in sel]]
    sub.index = [g for g, _ in sel]
    # 样本内 z-score
    z = sub.apply(lambda x: (x - x.mean()) / (x.std(ddof=1) + 1e-12), axis=1)
    s = z.mean(axis=0)
    return pd.Series(s, index=expr.columns), len(sel)

tip_score, tip_n = score_signature(tip)
core_score, core_n = score_signature(core)
breach_score, breach_n = score_signature(breach)
print(f"\ntip487 匹配 {tip_n} 基因；核心 {core_n}；breach {breach_n}")

# ---------- 5. 关联反应 ----------
def parse_grade(v):
    m = re.search(r"[:)]\s*(\d+)", str(v))
    if m:
        return int(m.group(1))
    m2 = re.search(r"(\d+)$", str(v).strip())
    return int(m2.group(1)) if m2 else None

meta["grade"] = meta[resp_key].apply(parse_grade) if resp_key else None
# 只保留肿瘤样本
if tissue_key:
    is_tumor = meta[tissue_key].astype(str).str.contains("tumor", case=False, na=False)
else:
    is_tumor = pd.Series(True, index=meta.index)
meta = meta[is_tumor].copy()
print("\n肿瘤样本数:", len(meta), "| grade 分布:", meta["grade"].value_counts().sort_index().to_dict())

meta["resp"] = meta["grade"].apply(lambda g: "responder" if g >= 2 else ("non-responder" if g in (0, 1) else None))

results = []
for label, sc in [("tip487", tip_score), ("tip_core", core_score), ("breach", breach_score)]:
    if sc is None:
        continue
    d = meta[["resp"]].join(sc.rename("score"))
    d = d[d["resp"].notna()]
    r = d[d.resp == "responder"].score
    n = d[d.resp == "non-responder"].score
    if len(r) < 3 or len(n) < 3:
        print(f"{label}: 样本不足 resp={len(r)} nonresp={len(n)}")
        continue
    U, p = stats.mannwhitneyu(r, n, alternative="two-sided")
    # AUC
    from sklearn.metrics import roc_auc_score
    y = (d.resp == "responder").astype(int)
    auc = roc_auc_score(y, d.score)
    results.append(dict(signature=label, n_resp=len(r), n_nonresp=len(n),
                        mean_resp=r.mean(), mean_nonresp=n.mean(),
                        wilcox_P=p, AUC=auc))
    print(f"\n[{label}] n={len(d)} (resp {len(r)}, nonresp {len(n)})")
    print(f"  responder {r.mean():.3f} vs non-responder {n.mean():.3f}, Mann-Whitney P={p:.4g}, AUC={auc:.3f}")

res = pd.DataFrame(results)
res.to_csv(os.path.join(OUT, "GSE225178_result.csv"), index=False, encoding="utf-8-sig")
print("\n=== 结果表 ===\n", res.to_string(index=False))

# ---------- 6. 图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
for label, sc in [("tip487", tip_score), ("tip_core", core_score), ("breach", breach_score)]:
    if sc is None:
        continue
    d = meta[["resp"]].join(sc.rename("score")).dropna(subset=["resp"])
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for g_, c_ in [("responder", "#b03a2e"), ("non-responder", "#4a9bc4")]:
        ax.boxplot(d[d.resp == g_].score, positions=[0 if g_ == "responder" else 1],
                   widths=0.4, patch_artist=True,
                   boxprops=dict(facecolor=c_, alpha=0.6))
    for i, g_ in enumerate(["responder", "non-responder"]):
        y = d[d.resp == g_].score
        ax.scatter(np.random.normal(i, 0.05, len(y)), y, s=18, alpha=0.7, color="#333")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Responder\n(Grade 2-3)", "Non-responder\n(Grade 0-1)"])
    ax.set_ylabel("Tip signature score (z)")
    from matplotlib.ticker import FormatStrFormatter
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    row = res[res.signature == label]
    if len(row):
        p = row.iloc[0].wilcox_P; auc = row.iloc[0].AUC
        ax.set_title(f"GSE225178: {label}\n(Mann-Whitney P = {p:.2f}, AUC = {auc:.2f})",
                       fontsize=9, pad=8)
    else:
        ax.set_title(f"GSE225178: {label}", fontsize=9, pad=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"GSE225178_{label}.png"), dpi=300, bbox_inches="tight")
    plt.close()

print("\nDONE")
