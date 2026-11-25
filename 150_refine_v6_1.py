# -*- coding: utf-8 -*-
"""v6 -> v6_1：第四轮复核意见的表述精修（不改任何科学内容/数值/图件）。
  T1 Methods stage family: add '21 tests in total' explicit
  T2 GO 措辞分层: six leading terms concordant (not 'all identical'); adjusted P explicit;
                 11,919 SYMBOL -> 10,836 unique Entrez (90.9%) defined
  T3 signature 数字关系: 487 / 241 / 235 / Jaccard 235/493 = 0.477 explicit
  T4 doublet 表述: pre-extracted matrix, complete unfiltered CD45-negative population unavailable
  T5 结果层次定位: composition = core; ECM theme = robust; 487 = exploratory gene set
  T6 整改记录分类: 4 项称「已确认的问题」而非都叫「统计学硬伤」
  T7 附表首页 README 加「当前有效表清单」；S8/S15 nestedCV 表标注当前/历史
  T8 校验报告措辞: 24 PASS 标为内部 QA，不进正式 Supplementary
"""
import io, os, re, shutil, zipfile
from docx import Document
import openpyxl

ROOT = r"D:\龙虾文档\2026-11-02-09-16-41"
BASE = os.path.join(ROOT, "ESCC_endo")
SUB = os.path.join(BASE, "12_submission")
V6 = os.path.join(SUB, "v6_bmcgenomics")
V61 = os.path.join(SUB, "v6_1_bmcgenomics")
os.makedirs(V61, exist_ok=True)

# ---------- copy v6 ----------
for f in os.listdir(V6):
    p = os.path.join(V6, f)
    if os.path.isfile(p): shutil.copy2(p, os.path.join(V61, f))
for sub in ("figures_tiff", "figures_png", "figures_pdf"):
    src, dst = os.path.join(V6, sub), os.path.join(V61, sub)
    if os.path.exists(src) and not os.path.exists(dst): shutil.copytree(src, dst)
print("[OK] copied v6 -> v6_1")

# ---------- source md refinements ----------
s = io.open(os.path.join(BASE, "11_manuscript", "manuscript_BMC_v6.md"), encoding="utf-8").read()
log = []
def rep(old, new, tag, expect=1):
    global s
    n = s.count(old)
    if n != expect: raise SystemExit("[FAIL] %s expect %d found %d :: %r" % (tag, expect, n, old[:120]))
    s = s.replace(old, new); log.append("%-56s x%d" % (tag, n))

# T1 stage family explicit total
rep("the seven state proportions against each of three stage variables (seven tests per variable), the retained CellChat pairs",
    "the seven state proportions against each of three stage variables, with Benjamini-Hochberg correction performed separately within each variable (seven tests per variable; 21 tests in total), the retained CellChat pairs",
    "T1 Methods 21 tests total")

# T2a Methods GO: precise universe definition + Entrez mapping + six-leading-terms wording
rep("As a sensitivity analysis, enrichment was repeated with the universe restricted to the 11,995 and 12,561 genes actually tested in the two tip-cell contrasts (intersection 11,919 genes; Table S16). The leading terms were identical under both backgrounds, so the functional interpretation does not depend on this choice.",
    "As a sensitivity analysis, enrichment was repeated with a custom universe restricted to the genes that passed the expression filter and actually entered differential testing in both tip-cell contrasts (12,561 genes in tip versus venous, 11,995 in tip versus capillary; 11,919 genes in the common tested universe; Table S16). After Symbol-to-Entrez conversion with org.Hs.eg.db and removal of duplicate identifiers, 10,836 unique Entrez identifiers (90.9% of the 11,919 input symbols) formed the usable background. The six leading GO biological-process terms were concordant between the two signature definitions and between the two background choices, so the functional interpretation does not depend on this choice.",
    "T2a Methods GO universe exact")

# T2b Results GO: explicit adjusted P + six leading terms + no 'all identical'
rep("Repeating enrichment with the universe restricted to the genes actually tested in the two contrasts gave the same first six terms for both the 487-gene and 241-gene signatures (extracellular matrix organization adjusted P = 1.9 x 10^-12 and 1.8 x 10^-11 respectively; Table S16); only the absolute P values changed.",
    "Repeating enrichment with the tested-gene universe gave concordant six leading biological-process terms for both the 487-gene and 241-gene signatures; for extracellular matrix organization the Benjamini-Hochberg-adjusted P values were 1.9 x 10^-12 and 1.8 x 10^-11 respectively (nominal P values were smaller; Table S16). The six leading terms were concordant, but this does not imply that every enrichment term was identical between backgrounds or signatures.",
    "T2b Results GO six leading terms")

# T3 signature numbers relationship (append to Discussion 487/241 sentence)
rep("re-running the pseudobulk contrasts with TMM normalization and voom precision weights shortened the shared signature from 487 to 241 genes while leaving the four leading Gene Ontology terms unchanged (Table S15).",
    "re-running the pseudobulk contrasts with TMM normalization and voom precision weights shortened the shared signature from 487 to 241 genes while leaving the four leading Gene Ontology terms unchanged (Table S15). The two gene sets share 235 genes, giving a Jaccard index of 235/(487+241-235) = 235/493 = 0.477: gene-level overlap was moderate, whereas the leading extracellular-matrix-related functional themes were consistent across the two differential-expression models.",
    "T3 signature 487/241/235/Jaccard")

# T4 doublet precise limitation statement
rep("No independent doublet-detection step (DoubletFinder, scDblFinder or equivalent) was run; potential doublets were addressed only indirectly through the gene-count, mitochondrial and contamination-cluster filters.",
    "The analysis was performed on a pre-extracted CD45-negative endothelial-cell matrix. Because the complete, unfiltered CD45-negative cell population was unavailable, a comprehensive sample-level doublet-detection procedure could not be performed; potential contamination was assessed indirectly through quality-control metrics, canonical lineage markers and exclusion of non-endothelial clusters.",
    "T4 doublet precise limitation")

# T5 result hierarchy in Discussion (prepend to tip gene list sentence)
rep("The tip gene list should therefore be read as a representative of an extracellular-matrix program, not a fixed set of markers:",
    "We therefore treat the venous, tip-like and matrix-remodeling composition differences as the core findings, the extracellular-matrix-related transcriptional theme as a comparatively robust functional result, and the 487-gene list as an exploratory gene set rather than a stable signature. The tip gene list should therefore be read as a representative of an extracellular-matrix program, not a fixed set of markers:",
    "T5 result hierarchy")

io.open(os.path.join(BASE, "11_manuscript", "manuscript_BMC_v6_1.md"), "w",
        encoding="utf-8", newline="\n").write(s)
for x in log: print("  [OK]", x)
print("[OK] source saved %d chars" % len(s))

# ---------- rebuild docx: copy 103+122 chain from 148 but point at v6_1 ----------
import subprocess, sys
# reuse 148's approach: derive from 145 with source/cover sections cut
src145 = io.open(os.path.join(BASE, "13_refcheck", "145_build_v5_bmcgenomics.py"), encoding="utf-8").read()
t = src145
t = t.replace('V43 = os.path.join(SUB, "v4_3_do")', 'V43 = os.path.join(SUB, "v6_bmcgenomics")')
t = t.replace('V5 = os.path.join(SUB, "v5_bmcgenomics")', 'V5 = os.path.join(SUB, "v6_1_bmcgenomics")')
t = t.replace('"manuscript_BMC_v5.md"', '"manuscript_BMC_v6_1.md"')
t = t.replace('"v5_bmcgenomics", "ESCC_BMC_Genomics_manuscript_v5.docx"',
              '"v6_1_bmcgenomics", "ESCC_BMC_Genomics_manuscript_v6_1.docx"')
t = t.replace('DOCX = os.path.join(V5, "ESCC_BMC_Genomics_manuscript_v5.docx")',
              'DOCX = os.path.join(V5, "ESCC_BMC_Genomics_manuscript_v6_1.docx")')
t = t.replace('"Figure_legends_BMC_v5.txt"', '"Figure_legends_BMC_v6_1.txt"')
t = t.replace('Figure_legends_BMC_v5.txt', 'Figure_legends_BMC_v6_1.txt')
t = t.replace('Cover_letter_BMC_v5.txt', 'Cover_letter_BMC_v6_1.txt')
t = t.replace('_tmp_103_v5.py', '_tmp_103_v61.py').replace('_tmp_122_v5.py', '_tmp_122_v61.py')
a = t.find('# ---------------------------------------------------------------- 1 源稿 5 项')
b = t.find('# ---------------------------------------------------------------- 2 主稿 docx')
assert a != -1 and b > a
t = t[:a] + '# source md prepared by 150\n' + t[b:]
c = t.find('# ---------------------------------------------------------------- 4 cover letter')
if c != -1: t = t[:c] + '\n# cover letter retained via copy\n'
sp = os.path.join(BASE, "13_refcheck", "_tmp_145_v61.py")
io.open(sp, "w", encoding="utf-8").write(t)
exec(compile(t, sp, "exec"), {"__name__": "__main__"})
os.remove(sp)
# restore source in case touched
io.open(os.path.join(BASE, "11_manuscript", "manuscript_BMC_v6_1.md"), "w",
        encoding="utf-8", newline="\n").write(s)

# ---------- table layout: repeating header + cantSplit + widths ----------
from docx.shared import Cm
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
DOCX = os.path.join(V61, "ESCC_BMC_Genomics_manuscript_v6_1.docx")
d = Document(DOCX)
WIDTHS = [Cm(2.7), Cm(2.5), Cm(2.9), Cm(2.1), Cm(4.6)]
for p in d.paragraphs:
    if p.text.strip().startswith("Table 2."): p.paragraph_format.keep_with_next = True
t2 = d.tables[1]; t2.autofit = False
for row in t2.rows:
    for j, c in enumerate(row.cells):
        if j < len(WIDTHS): c.width = WIDTHS[j]
for tb in d.tables:
    for ri, row in enumerate(tb.rows):
        trPr = row._tr.get_or_add_trPr()
        if trPr.find(qn("w:cantSplit")) is None: trPr.append(OxmlElement("w:cantSplit"))
        if ri == 0 and trPr.find(qn("w:tblHeader")) is None:
            h = OxmlElement("w:tblHeader"); h.set(qn("w:val"), "true"); trPr.append(h)
d.save(DOCX)

# ---------- supplementary: README current-sheet list + mark nestedCV current/historical ----------
XLSX_SRC = os.path.join(V61, "Supplementary_Tables_S1-S14_v6.xlsx")
XLSX_DST = os.path.join(V61, "Supplementary_Tables_S1-S14_v6_1.xlsx")
wb = openpyxl.load_workbook(XLSX_SRC)
ws = wb["README"]
# append current-effective table list after existing rows
start = ws.max_row + 2
lines = [
 ["CURRENT EFFECTIVE TABLES (v6_1, 2026-11-24)"],
 ["S1-S7", "Primary data, marker definitions, DE contrasts, CellChat retained pairs"],
 ["S8", "Reproducibility/robustness. NOTE: S8_nestedCV_folds is the ORIGINAL outer-fold workflow; the complete inner-loop re-selection variant is the CURRENT sensitivity workflow in S15."],
 ["S9-S14", "External validation, QC/identity checks, permutation tests, stage associations, response scores"],
 ["S15", "CURRENT model-sensitivity results: TMM+voom signature (487->241), GO leading themes, nested CV inner-loop re-selection"],
 ["S16", "CURRENT GO tested-gene-universe sensitivity: 11,919 input symbols -> 10,836 unique Entrez (90.9%); six leading terms concordant"],
 ["Archive", "Supplementary_Tables_archive_v3_5.xlsx is HISTORICAL and must not be uploaded."],
 ["QA note", "24-item automated check (校验报告) is internal QA only; not part of the Supplementary Material."],
]
for i, row in enumerate(lines):
    for j, v in enumerate(row): ws.cell(start+i, j+1, v)
# mark S8 historical fold sheet explicitly in its own sheet
if "S8_nestedCV_folds" in wb.sheetnames:
    w8 = wb["S8_nestedCV_folds"]; w8.cell(1, 10, "STATUS: ORIGINAL workflow (see S15 for current inner-loop re-selection)")
wb.save(XLSX_DST)
os.remove(XLSX_SRC)
print("[OK] workbook with current/historical labels: %d sheets" % len(wb.sheetnames))

# ---------- cleanup stale names, cover letter, points ----------
for fn in ["Cover_letter_DO_v4_3.docx","Cover_letter_DO_v4_3.txt","ESCC_DO_manuscript_v4_3.docx",
           "Figure_legends_DO_v4_3.txt","投稿操作要点_DO_v4_3.md","Cover_letter_BMC_v5.txt",
           "ESCC_BMC_Genomics_manuscript_v5.docx","Figure_legends_BMC_v5.txt","投稿要点_BMC_v5.md"]:
    p = os.path.join(V61, fn)
    if os.path.exists(p): os.remove(p)
# rename generated files from v5 to v6_1
for old,new in [("Cover_letter_BMC_v6_1.txt","Cover_letter_BMC_v6_1.txt")]: pass
cov = os.path.join(V61, "Cover_letter_BMC_v6_1.txt")
if os.path.exists(cov):
    x=io.open(cov,encoding="utf-8").read().replace("24 November 2026","24 November 2026").replace("v5","v6")
    io.open(cov,"w",encoding="utf-8",newline="\n").write(x)
# points file
pts = os.path.join(V61, "投稿要点_BMC_v5.md")
if os.path.exists(pts):
    x=io.open(pts,encoding="utf-8").read().replace("v5","v6_1").replace("BMC_v5","BMC_v6_1")
    io.open(os.path.join(V61,"投稿要点_BMC_v6_1.md"),"w",encoding="utf-8",newline="\n").write(x)
    os.remove(pts)
print("[OK] cleanup done")
print("[OK] FINAL:", os.path.basename(DOCX))
