# -*- coding: utf-8 -*-
"""v3.5: split Tables S1-S12 into standalone workbooks and verify every number
quoted in the manuscript against the underlying outputs.

Produces:
  12_submission/supplementary_tables_individual/Supplementary_Table_S{k}.xlsx
  12_submission/supplementary_tables_individual/INDEX.md
  13_refcheck/out_stats_v35/附表数值一致性核验_v3_5.md
"""
import os, io, re, csv, json, glob, collections, hashlib
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE, "12_submission")
OUT = os.path.join(SUB, "supplementary_tables_individual")
S12 = os.path.join(BASE, "13_refcheck", "out_stats_v35")
CC = r"D:\escc_sc\cellchat_out"
WB = os.path.join(SUB, "Supplementary_Tables_S1-S12_v3_5.xlsx")
EN = os.path.join(BASE, "11_manuscript", "manuscript_FI.md")
ZH = os.path.join(BASE, "11_manuscript", "manuscript_zh.md")
REPORT = os.path.join(S12, "附表数值一致性核验_v3_5.md")

os.makedirs(OUT, exist_ok=True)
log = []
ok = 0
fail = 0


def chk(cond, label, detail=""):
    global ok, fail
    if cond:
        ok += 1
        log.append("| PASS | %s | %s |" % (label, detail))
    else:
        fail += 1
        log.append("| **FAIL** | %s | %s |" % (label, detail))
    return cond


def norm(v):
    if v is None or v == "":
        return ""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float, np.generic)):
        try:
            return float(v)
        except Exception:
            return str(v)
    if isinstance(v, str):
        try:
            return float(v) if re.fullmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", v.strip()) else v
        except Exception:
            return v
    return str(v)


def same(a, b, tol=1e-9):
    a, b = norm(a), norm(b)
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) <= tol * max(1.0, abs(a), abs(b))
    return a == b


# ------------------------------------------------------------------ load
src = openpyxl.load_workbook(WB, read_only=True, data_only=True)
names = list(src.sheetnames)
groups = collections.OrderedDict((str(k), []) for k in range(1, 13))
for n in names:
    m = re.match(r"^S(\d+)", n)
    if m and 1 <= int(m.group(1)) <= 12:
        groups[str(int(m.group(1)))].append(n)
data = {n: [list(r) for r in src[n].iter_rows(values_only=True)] for n in names}
src.close()

chk(len(names) == 69, "工作簿 sheet 总数", "%d（预期 69）" % len(names))
chk(all(groups[k] for k in groups), "S1–S12 每组至少一张表",
    "; ".join("S%s=%d" % (k, len(v)) for k, v in groups.items()))

# ------------------------------------------------------------------ split
DESC = {
    "1": "Quality-control metrics per library",
    "2": "Marker genes of the nine endothelial states",
    "3": "Differential expression: tip-like versus venous",
    "4": "Differential expression: tip-like versus capillary",
    "5": "Differential expression: interferon-responsive versus venous",
    "6": "Differential expression: matrix-remodeling versus venous",
    "7": "Retained ligand-receptor pairs from the v3.5 CellChat re-run",
    "8": "Reproducibility, robustness and revised statistics",
    "9": "External single-cell directional validation",
    "10": "Cross-platform prognostic re-check in TCGA-ESCA",
    "11": "Legacy pooled-design sensitivity analyses",
    "12": "Independence-preserving statistical redesign",
}
index_lines = [
    "# 独立附表文件索引（v3.5）",
    "",
    "主工作簿：`../Supplementary_Tables_S1-S12_v3_5.xlsx`（68 张工作表）。",
    "下列 12 个文件与其中对应前缀的工作表内容逐格一致，供 Frontiers 投稿系统分别上传。",
    "",
    "| 文件 | 内容 | 工作表数 |",
    "| --- | --- | --- |",
]
for k, sheets in groups.items():
    path = os.path.join(OUT, "Supplementary_Table_S%s.xlsx" % k)
    w = openpyxl.Workbook()
    idx = w.active
    idx.title = "INDEX"
    idx.append(["Item", "Description"])
    idx.append(["Table", "Supplementary Table S%s" % k])
    idx.append(["Content", DESC[k]])
    idx.append(["Source workbook", os.path.basename(WB)])
    idx.append(["Sheets", len(sheets)])
    idx.append(["Note", "All values are identical to the corresponding sheets in the "
                        "master workbook. Machine-readable state labels are retained "
                        "for reproducibility; see S12_README for display-name mapping."])
    for s in sheets:
        idx.append([s, "%d rows (including header)" % len(data[s])])
    for c in idx[1]:
        c.font = Font(bold=True)
    idx.column_dimensions["A"].width = 34
    idx.column_dimensions["B"].width = 78
    for s in sheets:
        ws = w.create_sheet(s[:31])
        for r in data[s]:
            ws.append(["" if v is None else v for v in r])
        for c in ws[1]:
            c.font = Font(bold=True)
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "A2"
    w.save(path)
    index_lines.append("| `Supplementary_Table_S%s.xlsx` | %s | %d |" % (k, DESC[k], len(sheets) + 1))
    print("[OK] wrote Supplementary_Table_S%s.xlsx (%d sheets)" % (k, len(sheets) + 1))

index_lines += [
    "",
    "生成脚本：`13_refcheck/78_split_and_verify.py`。",
    "版本：v3.5（2026-11-09）。",
]
io.open(os.path.join(OUT, "INDEX.md"), "w", encoding="utf-8").write("\n".join(index_lines) + "\n")
chk(len(glob.glob(os.path.join(OUT, "Supplementary_Table_S*.xlsx"))) == 12,
    "独立附表文件数", "12 个（S1–S12）")

# verify split files cell-by-cell against master
split_ok = True
for k, sheets in groups.items():
    w = openpyxl.load_workbook(os.path.join(OUT, "Supplementary_Table_S%s.xlsx" % k),
                               read_only=True, data_only=True)
    for s in sheets:
        got = [list(r) for r in w[s[:31]].iter_rows(values_only=True)]
        if len(got) != len(data[s]) or any(
                not same(x, y) for a, b in zip(got, data[s]) for x, y in zip(a, b)):
            split_ok = False
            print("[FAIL] mismatch", k, s)
    w.close()
chk(split_ok, "独立文件与主工作簿逐格一致", "S1–S12 全部工作表")

# ------------------------------------------------------------------ S12 CSV <-> workbook
csv_ok = True
for fn in sorted(glob.glob(os.path.join(S12, "S12_*.csv"))):
    sheet = os.path.splitext(os.path.basename(fn))[0]
    rows = list(csv.reader(io.open(fn, encoding="utf-8-sig")))
    ref = data[sheet]
    if len(rows) != len(ref):
        csv_ok = False
        print("[FAIL] row count", sheet, len(rows), len(ref))
        continue
    for i, (a, b) in enumerate(zip(rows, ref)):
        for x, y in zip(a, b):
            if not same(x, y):
                # CSV numeric strings vs workbook floats handled by same()
                csv_ok = False
                print("[FAIL] cell", sheet, i + 1, repr(x), repr(y))
                break
chk(csv_ok, "S12 CSV 与工作簿逐格一致", "7 张数据表")

# ------------------------------------------------------------------ manuscript numbers
en = io.open(EN, encoding="utf-8").read()
zh = io.open(ZH, encoding="utf-8").read()
summ = json.load(io.open(os.path.join(S12, "S12_summary.json"), encoding="utf-8"))
pr = {r["state"]: r for r in csv.DictReader(
    io.open(os.path.join(S12, "S12_primary_56v4.csv"), encoding="utf-8-sig"))}

pairs = [
    ("Venous 癌旁 74.8%", "74.8%", float(pr["Venous"]["med_normal_pct"]), 74.8),
    ("Venous 肿瘤 28.5%", "28.5%", float(pr["Venous"]["med_tumor_pct"]), 28.5),
    ("Tip-like 14.8%", "14.8%", float(pr["Tip_angiogenic"]["med_tumor_pct"]), 14.8),
    ("Matrix 13.8%", "13.8%", float(pr["Matrix_EC"]["med_tumor_pct"]), 13.8),
    ("IFN 4.6%", "4.6%", float(pr["IFNgamma_EC"]["med_tumor_pct"]), 4.6),
]
for label, token, actual, expected in pairs:
    chk(abs(actual - expected) < 0.05 and token in en and token in zh,
        "正文/附表一致：" + label, "S12 = %.4f%%，中英文均含 %s" % (actual, token))

hl = float(pr["Venous"]["hodges_lehmann_pp"])
lo, hi = float(pr["Venous"]["ci95_low_pp"]), float(pr["Venous"]["ci95_high_pp"])
chk(abs(hl + 44.7881) < 1e-3 and "-44.8" in en.replace("−", "-") and "−44.8" in zh,
    "Hodges–Lehmann 位移 −44.8 pp", "S12 = %.4f" % hl)
chk(abs(lo + 61.2128) < 1e-3 and abs(hi + 28.6822) < 1e-3,
    "bootstrap 95% CI −61.2 至 −28.7", "S12 = %.4f 至 %.4f" % (lo, hi))
for st, q in [("Venous", 0.0005), ("Tip_angiogenic", 0.0007),
              ("Matrix_EC", 0.0005), ("IFNgamma_EC", 0.0007)]:
    qa = float(pr[st]["q_exact"])
    chk(abs(qa - q) < 0.00006, "BH q（%s）" % st, "S12 = %.6f，正文约 %.4f" % (qa, q))
chk(abs(float(pr["Capillary"]["q_exact"]) - 0.0592) < 0.001,
    "Capillary 主分析不显著", "q = %.6f" % float(pr["Capillary"]["q_exact"]))

loo = list(csv.DictReader(io.open(os.path.join(S12, "S12_loo_bh_family.csv"), encoding="utf-8-sig")))
chk(len(loo) == 16 and all(r["direction_preserved"] in ("True", "true", "1") for r in loo),
    "逐一剔除癌旁 16/16 方向保留", "n = %d" % len(loo))
q16 = [float(r["q_joint_16"]) for r in loo]
chk(abs(min(q16) - 0.000984) < 1e-5 and abs(max(q16) - 0.002522) < 1e-5,
    "联合 BH q 0.0010–0.0025", "实测 %.6f–%.6f" % (min(q16), max(q16)))

ds = {(r["design"], r["cells_per_sample"], r["state"]): float(r["fdr_pass_rate"])
      for r in csv.DictReader(io.open(os.path.join(S12, "S12_downsampling.csv"), encoding="utf-8-sig"))}
chk(abs(ds[("primary 56 vs 4", "20", "IFNgamma_EC")] - 0.160) < 1e-9 and
    abs(ds[("primary 56 vs 4", "30", "IFNgamma_EC")] - 0.855) < 1e-9,
    "IFN 下采样 16% / 85.5%", "20 cells = %.3f；30 cells = %.3f" %
    (ds[("primary 56 vs 4", "20", "IFNgamma_EC")], ds[("primary 56 vs 4", "30", "IFNgamma_EC")]))

sp = list(csv.DictReader(io.open(os.path.join(S12, "S12_sample_pairing.csv"), encoding="utf-8-sig")))
chk(len([r for r in sp if r["included_in_primary"] in ("True", "true", "1")]) == 60,
    "主分析纳入 56 肿瘤 + 4 癌旁", "共 %d 个纳入" %
    len([r for r in sp if r["included_in_primary"] in ("True", "true", "1")]))
chk(sorted(summ["design"]["paired_tumors"]) ==
    ["P126T", "P127T", "P128T", "P130T"],
    "4 对配对样本 ID", ", ".join(summ["design"]["paired_tumors"]))

# ------------------------------------------------------------------ CellChat
net = list(csv.DictReader(io.open(os.path.join(CC, "cellchat_net.csv"), encoding="utf-8")))
EC = {"Arterial", "Capillary", "Venous", "Venous_inflamed", "Tip_angiogenic",
      "Matrix_EC", "IFNgamma_EC", "Proliferating", "Lymphatic"}
T5 = {"CD4 T", "CD8 T", "CD8 exhausted", "NK-like", "Treg"}
e2t = [r for r in net if r["source"] in EC and r["target"] in T5]
chk(len(net) == 6132 and len(data["S7_CellChat_interactions"]) == 6133,
    "保留配体–受体对 6,132", "CSV = %d；S7 表 = %d 行（含表头）" %
    (len(net), len(data["S7_CellChat_interactions"])))
chk(len(e2t) == 1556 and "1,556" in en and "1,556" in zh,
    "内皮→T 配对 1,556", "CSV = %d" % len(e2t))

mhc = [r for r in e2t if r["pathway_name"] == "MHC-I"]
mif = [r for r in e2t if r["pathway_name"] == "MIF"]
chk(len(mhc) == 180 and abs(sum(float(r["prob"]) for r in mhc) - 32.2933) < 0.001,
    "MHC-I 180 对 / 32.3 任意单位", "pairs = %d；score = %.4f" %
    (len(mhc), sum(float(r["prob"]) for r in mhc)))
chk(len(mif) == 80 and abs(sum(float(r["prob"]) for r in mif) - 13.7080) < 0.001,
    "MIF 80 对 / 13.7 任意单位", "pairs = %d；score = %.4f" %
    (len(mif), sum(float(r["prob"]) for r in mif)))

tcc = {r["lab"]: int(r["Freq"]) for r in csv.DictReader(
    io.open(os.path.join(CC, "tcell_group_counts_v3_5.csv"), encoding="utf-8"))}
chk(tcc.get("Other/unclassified T") == 2109 and "2,109" in en and "2,109" in zh,
    "Other/unclassified T = 2,109", "计数表 = %d；合计 %d" %
    (tcc.get("Other/unclassified T", 0), sum(tcc.values())))

rob = []
for s in (1, 2):
    rob.append(list(csv.DictReader(io.open(
        os.path.join(CC, "robust_v35", "summary_seed%d.csv" % s), encoding="utf-8")))[0])
chk(all(r["top1"] == "MHC-I" and r["top2"] == "MIF" and r["top3"] == "APP" for r in rob),
    "200 次置换双 seed 前三通路不变", "; ".join(
        "seed%s: %s/%s/%s" % (r["seed"], r["top1"], r["top2"], r["top3"]) for r in rob))

# ------------------------------------------------------------------ S11 <-> S12 legacy design
s11 = list(csv.DictReader(io.open(
    os.path.join(BASE, "13_refcheck", "out_sensitivity", "S11_baseline_group_tests.csv"),
    encoding="utf-8-sig")))
s12s = {r["state"]: r for r in csv.DictReader(
    io.open(os.path.join(S12, "S12_sensitivity_60v4.csv"), encoding="utf-8-sig"))}
legacy_ok = all(abs(float(r["p_asymptotic"]) -
                    float(s12s[r["state"]]["p_asymptotic"])) < 1e-12 for r in s11)
chk(legacy_ok, "S11 旧池化设计与 S12 敏感性 60v4 的渐近 P 值一致",
    "9 个状态逐项比对")

# ------------------------------------------------------------------ stale terms
stale_en = ["None of these studies compared", "venous depletion", "venous-to-tip exchange",
            "All 60 tumors and 4 adjacent samples entered the test", "5,758", "5,692",
            "significant interactions", "Total communication probability",
            "default CD4 T label", "IFN-gamma N<T", "IFNgamma_EC fraction"]
stale_zh = ["这些研究都未在肿瘤与配对癌旁黏膜之间比较", "全部 60 例肿瘤与 4 例癌旁均纳入检验",
            "5,758", "5,692", "默认的 CD4 T 标签", "静脉型向尖端细胞的交换",
            "聚合通讯概率", "IFN-gamma N<T"]
# Exclude the reference list: article titles such as "IFN-gamma-related mRNA
# profile ..." (ref. 29) must stay verbatim and are not stale state labels.
en_body = re.split(r"(?m)^##\s*References\b", en)[0]
zh_body = re.split(r"(?m)^##\s*参考文献", zh)[0]
for label, text, terms in [("英文", en_body, stale_en), ("中文", zh_body, stale_zh)]:
    bad = [t for t in terms if t in text]
    chk(not bad, "%s稿无遗留旧表述" % label, "、".join(bad) if bad else "0 处")

# ------------------------------------------------------------------ Frontiers checks
ab = en.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0].strip()
nwords = len(re.findall(r"[A-Za-z0-9'’-]+", ab))
chk(nwords <= 350, "英文摘要 ≤350 词", "%d 词" % nwords)
chk(len(re.findall(r"(?m)^\d+\. ", en)) == 46 and
    len(re.findall(r"(?m)^\d+\. ", zh)) == 46, "中英文参考文献均 46 条", "46 / 46")
chk([c[:11] for c in re.findall(r"(?m)^\*\*Table \d+\.", en)] ==
    ["**Table 1.", "**Table 2.", "**Table 3."], "英文表题物理顺序 1→2→3", "1, 2, 3")
chk(not re.search(r"TODO|TBD|XXX|\{\{|\[insert", en + zh), "无占位符", "0 处")
rt = re.search(r"(?m)^\*\*Running title:\*\* (.+)$", en)
chk(rt and len(rt.group(1).split()) <= 5, "Running title ≤5 词",
    rt.group(1) if rt else "未找到")
chk("Table S12" in en and "附表 S12" in zh, "补充材料声明含 S12", "中英文均有")
chk("Each table is also supplied as a standalone file" in en and "各附表另以独立文件提供" in zh,
    "独立附表交付声明与实物一致", "12 个文件已生成")

# ------------------------------------------------------------------ figures
figdir = os.path.join(SUB, "figures")
tiffdir = os.path.join(SUB, "figures_tiff")
pngs = sorted(glob.glob(os.path.join(figdir, "Figure*.png")))
tifs = sorted(glob.glob(os.path.join(tiffdir, "*.tif")))
chk(len(pngs) >= 11, "投稿图件 PNG 数量", "%d 张" % len(pngs))
chk(len(tifs) == 11, "TIFF 图件数量", "%d 张" % len(tifs))
from PIL import Image
dpi_ok = True
for t in tifs:
    im = Image.open(t)
    d = im.info.get("dpi", (0, 0))
    if im.mode != "RGB" or min(d) < 300:
        dpi_ok = False
        print("[FAIL] TIFF", os.path.basename(t), im.mode, d)
chk(dpi_ok, "TIFF 均为 RGB 且 ≥300 dpi", "11 张")

# ------------------------------------------------------------------ report
header = [
    "# 附表 S1–S12 与正文数值一致性核验（v3.5）",
    "",
    "执行日期：2026-11-09。脚本：`13_refcheck/78_split_and_verify.py`。",
    "",
    "核验对象：",
    "- 主工作簿 `12_submission/Supplementary_Tables_S1-S12_v3_5.xlsx`（69 张工作表）",
    "- 独立附表 `12_submission/supplementary_tables_individual/Supplementary_Table_S1–S12.xlsx`",
    "- 英文源稿 `11_manuscript/manuscript_FI.md`、中文源稿 `11_manuscript/manuscript_zh.md`",
    "- 统计输出 `13_refcheck/out_stats_v35/`、CellChat 输出 `D:/escc_sc/cellchat_out/`",
    "",
    "**结论：%d 项通过，%d 项未通过。**" % (ok, fail),
    "",
    "| 结果 | 核验项 | 实测 |",
    "| --- | --- | --- |",
]
io.open(REPORT, "w", encoding="utf-8").write("\n".join(header + log + [
    "",
    "## 说明",
    "",
    "- S11 使用旧池化 60 对 4 设计与渐近检验，保留为稳健性检查；v3.5 主分析为 S12 的 56 例未配对肿瘤对 4 例癌旁的精确检验。两者渐近 P 值逐项一致，设计差异已在 S11_README 与 S12_README 中标注。",
    "- S7 与 S8 的 CellChat 工作表已按新增 `Other/unclassified T` 的重跑结果替换；旧版保存在 `Supplementary_Tables_S1-S10_pre_S12_20261109.xlsx`。",
    "- 正文 1,556 对内皮→T 配对来自 nboot = 50 的主分析；S8 的 200 次置换双 seed 重复分别为 1,558 与 1,551 对，属不同运行，不构成矛盾。",
    "- 附表中保留机器标签（如 `IFNgamma_EC`）以便复现；稿件与图件使用显示名，映射见 S12_README。",
]) + "\n")
print("[REPORT]", REPORT)
print("[SUMMARY] %d passed, %d failed" % (ok, fail))
if fail:
    raise SystemExit("[FAIL] consistency check has %d failed items" % fail)
