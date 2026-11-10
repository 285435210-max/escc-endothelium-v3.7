# -*- coding: utf-8 -*-
"""v3.7 comprehensive verification."""
import os, io, re, json, csv, glob, hashlib, zipfile
import openpyxl
from docx import Document
from docx.oxml.ns import qn
from PIL import Image

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE, "12_submission")
OUT = os.path.join(BASE, "13_refcheck", "out_v37")
IND = os.path.join(SUB, "supplementary_tables_individual")
REPORT = os.path.join(OUT, "v3_7_综合校验.md")
C = []


def chk(ok, name, detail=""):
    C.append((bool(ok), name, detail))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(65536), b""):
            h.update(c)
    return h.hexdigest()


def rows(p):
    with io.open(p, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


# workbook
p = os.path.join(SUB, "Supplementary_Tables_S1-S14_v3_7.xlsx")
wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
n = len(wb.sheetnames)
s14 = [x for x in wb.sheetnames if x.startswith("S14")]
wb.close()
chk(n == 85, "主工作簿 sheet 总数", "%d（预期 85）" % n)
chk(len(s14) == 8, "S14 工作表数", "8")
chk(len(glob.glob(os.path.join(IND, "Supplementary_Table_S*.xlsx"))) == 14,
    "独立附表数量", "14")

# identity check
co = rows(os.path.join(OUT, "inflamed_coexpression.csv"))
allrow = [r for r in co if r["subset"] == "all Venous_inflamed"][0]
chk(float(allrow["pct_with_ge2_panEC"]) >= 99.0, "炎症小静脉 ≥2 泛内皮", allrow["pct_with_ge2_panEC"] + "%")
chk(float(allrow["pct_with_ge2_activation"]) >= 99.0, "≥2 活化标志", allrow["pct_with_ge2_activation"] + "%")
chk(float(allrow["pct_with_ge1_leukocyte"]) < 3.0, "白细胞标志污染", allrow["pct_with_ge1_leukocyte"] + "%")
chk(float(allrow["pct_with_ge1_erythroid"]) == 0.0, "红细胞标志污染", "0%")
mk = rows(os.path.join(OUT, "inflamed_marker_by_state.csv"))
inf = {r["gene"]: float(r["pct_detected"]) for r in mk if r["state"] == "Venous_inflamed"}
for g, v in [("PECAM1", 99.6), ("VWF", 95.9), ("SELE", 91.1), ("ICAM1", 97.0)]:
    chk(abs(inf[g] - v) < 0.2, "标志检出率 %s" % g, "%.1f%%" % inf[g])

# stage + chemo
stage = rows(os.path.join(OUT, "S14_stage_7states.csv"))
chk(len(stage) == 21, "七状态分期检验数", "21")
chk(min(float(r["q_BH_within_variable"]) for r in stage) > 0.05, "分期无显著关联",
    "min q = %.3f" % min(float(r["q_BH_within_variable"]) for r in stage))
chemo = rows(os.path.join(OUT, "S14_chemo_scores.csv"))
chk(len(chemo) == 3, "化疗三评分", "3")
chk(all(0.54 < float(r["auc"]) < 0.58 for r in chemo), "化疗 AUC 区间",
    ", ".join("%.3f" % float(r["auc"]) for r in chemo))

# manuscripts
en = io.open(os.path.join(BASE, "11_manuscript", "manuscript_FI.md"), encoding="utf-8").read()
zh = io.open(os.path.join(BASE, "11_manuscript", "manuscript_zh.md"), encoding="utf-8").read()
body = en.split("## References")[0]
ab = en.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0].strip()
aw = len(re.findall(r"[A-Za-z0-9'’-]+", ab))
chk(300 <= aw <= 350, "英文摘要词数", "%d" % aw)
seg = en.split("## Introduction", 1)[1].split("## Data availability statement", 1)[0]
wc = len(" ".join(l.strip() for l in seg.split("\n")
                  if l.strip() and not l.strip().startswith(("#", "|", "![", "**Figure"))
                  and l.strip() != "---").split())
chk(8500 < wc < 12000, "正文词数", "%d" % wc)
refs = [int(x) for x in re.findall(r"^(\d+)\. ", en, re.M)]
chk(refs == list(range(1, 47)), "参考文献 1..46", "46")
chk(re.findall(r"\*\*Table (\d)\.", en) == ["1", "2", "3"], "表题顺序", "1->2->3")
chk(en.count("Table S14") >= 5, "英文稿 Table S14", "%d 处" % en.count("Table S14"))
chk("附表 S14" in zh, "中文稿附表 S14", "")
for lab, txt, terms in [
    ("EN", body, ["primary paired design", "ranks pairs by ligand abundance",
                  "no independent paired ESCC adjacent-normal cohort was analysed",
                  "The figures, supplementary tables and analysis scripts retain",
                  "Tables S1-S12, Figures", "inflamed venule (SELL, ICAM4),",
                  "tumour", "analysed", "binarised", "oesophago-gastric"]),
    ("ZH", zh, ["与主要配对设计分开", "图件、附表与分析脚本沿用",
               "炎症小静脉（SELL、ICAM4）", "附表 S1–S13、附图"])]:
    stale = [t for t in terms if t in txt]
    chk(not stale, "%s 旧表述/英式拼写清零" % lab, ", ".join(stale) if stale else "0 残留")
chk("99.9%" in en and "99.9%" in zh, "身份核查数字中英文一致", "")
chk("487,635" in en, "置换枚举数保留", "")

# docx
DOCS = [("英文主件", os.path.join(SUB, "ESCC_FI_manuscript_v3_7.docx"), 0),
        ("英文图文版", os.path.join(SUB, "ESCC_FI_manuscript_withfigs_v3_7.docx"), 11),
        ("中文同步版", os.path.join(BASE, "11_manuscript", "ESCC_endothelium_manuscript_zh_v3_7.docx"), 11)]
for lab, path, nfig in DOCS:
    d = Document(path)
    chk(len(d.tables) == 3, "%s 表格数" % lab, "3")
    chk(len(d.inline_shapes) == nfig, "%s 图片数" % lab, "%d" % len(d.inline_shapes))
    cap = [x for x in d.paragraphs if x.text.strip().startswith("Table 2.")
           or x.text.strip().startswith("表 2")]
    chk(len(cap) == 1 and cap[0].paragraph_format.keep_with_next is True
        and cap[0].paragraph_format.page_break_before is True,
        "%s Table 2 分页" % lab, "")
    t2 = d.tables[1]
    chk(all(r._tr.get_or_add_trPr().find(qn("w:cantSplit")) is not None for r in t2.rows),
        "%s Table 2 行不拆分" % lab, "")
    chk(round(t2.rows[0].cells[4].width.cm, 1) == 4.6, "%s Table 2 AUC 列宽" % lab, "4.6 cm")
    concl = [x for x in d.paragraphs if x.text.strip().startswith(
        "Based on four adjacent non-malignant libraries,")
        or x.text.strip().startswith("基于 4 例癌旁非恶性文库")]
    chk(len(concl) == 1 and concl[0].paragraph_format.keep_together is True,
        "%s Conclusion 不跨页断词" % lab, "")

with zipfile.ZipFile(DOCS[1][1]) as z:
    hashes = set(hashlib.md5(z.read(x)).hexdigest() for x in z.namelist()
                 if x.startswith("word/media/") and x.lower().endswith(".png"))
cur = set(md5(x) for x in glob.glob(os.path.join(SUB, "figures", "Figure*.png")))
chk(len(hashes & cur) >= 11, "图文版嵌入图为 v3.7 新图",
    "%d/%d" % (len(hashes & cur), len(cur)))

tifs = sorted(glob.glob(os.path.join(SUB, "figures_tiff", "*.tif")))
chk(len(tifs) == 11, "投稿 TIFF 数量", "11")
chk(all(Image.open(x).mode == "RGB" for x in tifs), "TIFF 全部 RGB", "")
for b in ["Figure1", "Figure2", "Figure6"]:
    src = {"Figure1": "Fig1_atlas.png", "Figure2": "Fig2_composition.png",
           "Figure6": "Fig5_clinical.png"}[b]
    chk(md5(os.path.join(SUB, "figures", b + ".png")) == md5(os.path.join(r"D:/escc_sc/figs", src)),
        "%s 与重绘源一致" % b, "")

np_ = sum(1 for ok, _, _ in C if ok)
nf = len(C) - np_
lines = ["# v3.7 综合校验报告（2026-11-09）", "",
         "**结果：%d 项 PASS / %d 项 FAIL**" % (np_, nf), "",
         "| # | 检查项 | 结果 | 明细 |", "| --- | --- | --- | --- |"]
for i, (ok, name, detail) in enumerate(C, 1):
    lines.append("| %d | %s | %s | %s |" % (i, name, "PASS" if ok else "**FAIL**", detail))
lines += ["", "生成脚本：`13_refcheck/107_verify_v3_7.py`。"]
io.open(REPORT, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
print("[OK] %d PASS / %d FAIL" % (np_, nf))
for ok, name, detail in C:
    if not ok:
        print("  [FAIL]", name, detail)
