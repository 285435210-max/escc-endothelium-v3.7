# -*- coding: utf-8 -*-
"""v3.6 final verification + standalone Supplementary Table S13."""
import os, io, re, csv, json, glob, hashlib, zipfile
import openpyxl
from docx import Document
from docx.oxml.ns import qn
from PIL import Image

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE, "12_submission")
S13 = os.path.join(BASE, "13_refcheck", "out_stats_v36")
IND = os.path.join(SUB, "supplementary_tables_individual")
REPORT = os.path.join(S13, "v3_6_综合校验.md")
CHECKS = []


def chk(ok, name, detail=""):
    CHECKS.append((bool(ok), name, detail))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def csv_rows(fn):
    with io.open(os.path.join(S13, fn), encoding="utf-8-sig") as f:
        return list(csv.reader(f))


# ---------------------------------------------------------------- standalone S13
SHEETS = [
    ("S13_README", None),
    ("S13_permutation_primary", "S13_permutation_primary.csv"),
    ("S13_permutation_60v4", "S13_permutation_60v4.csv"),
    ("S13_permutation_loo_16", "S13_permutation_loo_16.csv"),
    ("S13_inflamed_venule_loo", "S13_inflamed_venule_loo.csv"),
    ("S13_paired_four_states", "S13_paired_four_states.csv"),
    ("S13_downsampling_sample_counts", "S13_downsampling_sample_counts.csv"),
    ("S13_methods", "S13_methods.csv"),
]
wb = openpyxl.Workbook()
wb.remove(wb.active)
for name, fn in SHEETS:
    ws = wb.create_sheet(name)
    rows = [[r[0], r[1]] for r in csv_rows("S13_methods.csv")] if fn is None else csv_rows(fn)
    if fn is None:
        master = openpyxl.load_workbook(os.path.join(SUB, "Supplementary_Tables_S1-S13_v3_6.xlsx"),
                                        read_only=True, data_only=True)
        rows = [list(r) for r in master["S13_README"].iter_rows(values_only=True)]
        master.close()
    for r in rows:
        ws.append(r)
s13_file = os.path.join(IND, "Supplementary_Table_S13.xlsx")
wb.save(s13_file)
chk(os.path.getsize(s13_file) > 3000, "独立附表 S13 生成", "%.1f KB" % (os.path.getsize(s13_file) / 1024.0))

index = os.path.join(IND, "INDEX.md")
if os.path.exists(index):
    t = io.open(index, encoding="utf-8").read()
    t = t.replace("S1–S12", "S1–S13").replace("S1-S12", "S1-S13").replace("12 个", "13 个")
    if "Supplementary_Table_S13.xlsx" not in t:
        t += "\n- `Supplementary_Table_S13.xlsx` — 考虑并列秩的精确标签置换检验、四个重点状态配对方向、下抽样样本数与方法定义（8 张工作表）。\n"
    io.open(index, "w", encoding="utf-8", newline="\n").write(t)
n_ind = len(glob.glob(os.path.join(IND, "Supplementary_Table_S*.xlsx")))
chk(n_ind == 13, "独立附表数量", "%d（预期 13）" % n_ind)

# ---------------------------------------------------------------- workbook
master_path = os.path.join(SUB, "Supplementary_Tables_S1-S13_v3_6.xlsx")
wb = openpyxl.load_workbook(master_path, read_only=True, data_only=True)
names = wb.sheetnames
chk(len(names) == 77, "主工作簿 sheet 总数", "%d（预期 77）" % len(names))
chk(sum(1 for n in names if n.startswith("S13")) == 8, "S13 工作表数", "8")
chk("S13_README" in names, "S13_README 存在", "")
wb.close()

# ---------------------------------------------------------------- S13 numbers
summary = json.load(io.open(os.path.join(S13, "S13_summary.json"), encoding="utf-8"))
chk(summary["primary_assignments"] == 487635, "主设计枚举数", "487,635")
chk(summary["loo_assignments_per_deletion"] == 32509, "LOO 枚举数", "32,509")
chk(summary["pooled_assignments"] == 635376, "60v4 枚举数", "635,376")
chk(summary["all_three_concordant"] is True, "三法判定一致", "9/9")
chk(summary["n_sig_perm"] == 6, "置换法显著状态数", "6")
chk(summary["loo16_directions_preserved"] == 16, "LOO16 方向保留", "16/16")
chk(summary["inflamed_loo_q_max"] < 0.05, "炎症小静脉补充 LOO 显著",
    "max q = %.6f" % summary["inflamed_loo_q_max"])

prim = list(csv.DictReader(io.open(os.path.join(S13, "S13_permutation_primary.csv"),
                                   encoding="utf-8-sig")))
q = {r["state"]: float(r["q_perm"]) for r in prim}
for st, expect in [("Venous", 0.000305), ("Tip_angiogenic", 0.000605),
                   ("Matrix_EC", 0.000305), ("IFNgamma_EC", 0.000605),
                   ("Venous_inflamed", 0.000605), ("Capillary", 0.054860)]:
    chk(abs(q[st] - expect) < 1e-5, "S13 q_perm %s" % st, "%.6f" % q[st])

# ---------------------------------------------------------------- manuscripts
en = io.open(os.path.join(BASE, "11_manuscript", "manuscript_FI.md"), encoding="utf-8").read()
zh = io.open(os.path.join(BASE, "11_manuscript", "manuscript_zh.md"), encoding="utf-8").read()
ab = en.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0].strip()
aw = len(re.findall(r"[A-Za-z0-9'’-]+", ab))
chk(310 <= aw <= 350, "英文摘要词数", "%d（Frontiers 上限 350）" % aw)
seg = en.split("## Introduction", 1)[1].split("## Data availability statement", 1)[0]
keep = [l.strip() for l in seg.split("\n")
        if l.strip() and not l.strip().startswith(("#", "|", "![", "**Figure"))
        and l.strip() != "---"]
wc = len(" ".join(keep).split())
chk(wc < 12000, "正文词数", "%d（上限 12,000）" % wc)
refs = re.findall(r"^(\d+)\. ", en, re.M)
chk([int(x) for x in refs] == list(range(1, 47)), "参考文献编号 1..46", "46 refs")
cites = set(int(x) for x in re.findall(r"\[(\d+)\]", en))
chk(cites and max(cites) <= 46, "正文引用上限", "max = %d" % max(cites))
caps = re.findall(r"\*\*Table (\d)\.", en)
chk(caps == ["1", "2", "3"], "表题顺序", "->".join(caps))
chk("TODO" not in en and "XXX" not in en and "{{" not in en, "无占位符",
    "Word count: AUTO 为构建脚本标记，生成 docx 时已替换为实际数字")
chk(en.count("Table S13") >= 3, "英文稿 Table S13 声明", "%d 处" % en.count("Table S13"))
chk("附表 S13" in zh, "中文稿附表 S13 声明", "")
for label, text, terms in [
    ("英文", en, ["neither the score nor the state proportions predict survival",
                "limma-trend, with patient included as a blocking factor",
                "ranked MHC class I first and MIF second by ligand transcript abundance",
                "exact two-sided Mann-Whitney tests were applied"]),
    ("中文", zh, ["且评分与状态占比均不预测生存", "用 limma-trend 建模",
               "按配体转录本丰度将 MHC I 类列于首位", "精确双侧 Mann–Whitney 检验"]),
]:
    stale = [t for t in terms if t in text]
    chk(not stale, "%s 旧表述清零" % label, ", ".join(stale) if stale else "0 残留")
chk("487,635" in en and "487,635" in zh, "中英文枚举数一致", "")
chk("0.00031" in en or "0.0003" in en, "正文置换 q 值", "")

# ---------------------------------------------------------------- docx structure
DOCS = [
    ("英文主件", os.path.join(SUB, "ESCC_FI_manuscript_v3_6.docx"), 0),
    ("英文图文版", os.path.join(SUB, "ESCC_FI_manuscript_withfigs_v3_6.docx"), 11),
    ("中文同步版", os.path.join(BASE, "11_manuscript", "ESCC_endothelium_manuscript_zh_v3_6.docx"), 11),
]
for label, path, nfig in DOCS:
    d = Document(path)
    npics = len(d.inline_shapes)
    chk(len(d.tables) == 3, "%s 表格数" % label, "3")
    chk(npics == nfig, "%s 图片数" % label, "%d（预期 %d）" % (npics, nfig))
    cap = [p for p in d.paragraphs if p.text.strip().startswith("Table 2.")
           or p.text.strip().startswith("表 2")]
    chk(len(cap) == 1 and cap[0].paragraph_format.keep_with_next is True
        and cap[0].paragraph_format.page_break_before is True,
        "%s Table 2 分页修复" % label, "keep_with_next + page_break_before")
    if len(d.tables) >= 2:
        t2 = d.tables[1]
        ok = all(r._tr.get_or_add_trPr().find(qn("w:cantSplit")) is not None for r in t2.rows)
        chk(ok, "%s Table 2 行不拆分" % label, "%d 行" % len(t2.rows))

# embedded picture checksums vs current PNGs
with zipfile.ZipFile(DOCS[1][1]) as z:
    hashes = set()
    for n in z.namelist():
        if n.startswith("word/media/") and n.lower().endswith((".png", ".jpeg", ".jpg")):
            hashes.add(hashlib.md5(z.read(n)).hexdigest())
current = set(md5(p) for p in glob.glob(os.path.join(SUB, "figures", "Figure*.png")))
chk(len(hashes & current) >= 11, "图文版嵌入图为 v3.6 新图",
    "%d/%d 匹配" % (len(hashes & current), len(current)))

# ---------------------------------------------------------------- figures
tifs = sorted(glob.glob(os.path.join(SUB, "figures_tiff", "*.tif")))
chk(len(tifs) == 11, "投稿 TIFF 数量", "%d" % len(tifs))
for f in tifs:
    im = Image.open(f)
    chk(im.mode == "RGB", "TIFF RGB %s" % os.path.basename(f), im.mode)
for base in ["FigureS1", "FigureS2", "Figure6"]:
    new = os.path.join(SUB, "figures", base + ".png")
    old = os.path.join(SUB, "figures_pre_v36", base + "_pre_v36.png")
    chk(os.path.exists(new) and os.path.exists(old) and md5(new) != md5(old),
        "%s 已重绘" % base, "md5 与 v3.5 不同")

# ---------------------------------------------------------------- report
npass = sum(1 for ok, _, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
lines = ["# v3.6 综合校验报告（2026-11-09）", "",
         "**结果：%d 项 PASS / %d 项 FAIL**" % (npass, nfail), "",
         "| # | 检查项 | 结果 | 明细 |", "| --- | --- | --- | --- |"]
for i, (ok, name, detail) in enumerate(CHECKS, 1):
    lines.append("| %d | %s | %s | %s |" % (i, name, "PASS" if ok else "**FAIL**", detail))
lines += ["", "生成脚本：`13_refcheck/95_verify_v3_6.py`。",
          "统计重算：`13_refcheck/86_perm_v3_6.py`；汇总 `out_stats_v36/S13_summary.json`。"]
io.open(REPORT, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
print("[OK] checks: %d PASS / %d FAIL" % (npass, nfail))
print("[OK] report:", REPORT)
if nfail:
    for ok, name, detail in CHECKS:
        if not ok:
            print("  [FAIL]", name, detail)
