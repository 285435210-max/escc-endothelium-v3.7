# -*- coding: utf-8 -*-
"""v7 -> v8_bmcgenomics: complete AI disclosure (add second AI reviewer per external audit)
   + prepare v1.1.0 code bundle including GO universe script (147) and reproducibility chain.
   NO science changes.
"""
import io, os, re, shutil, hashlib, csv, zipfile
BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE,"12_submission")
V7 = os.path.join(SUB,"v7_bmcgenomics")
V8 = os.path.join(SUB,"v8_bmcgenomics")
os.makedirs(V8,exist_ok=True)
for f in os.listdir(V7):
    p=os.path.join(V7,f)
    if os.path.isfile(p): shutil.copy2(p,os.path.join(V8,f))
for sub in ("figures_tiff","figures_png","figures_pdf"):
    src,dst=os.path.join(V7,sub),os.path.join(V8,sub)
    if os.path.exists(src) and not os.path.exists(dst): shutil.copytree(src,dst)
print("[OK] copied v7 -> v8")

s=io.open(os.path.join(BASE,"11_manuscript","manuscript_BMC_v7.md"),encoding="utf-8").read()
log=[]
def rep(old,new,tag):
    global s
    assert s.count(old)==1,(tag,s.count(old))
    s=s.replace(old,new); log.append(tag)

AI_M = ("A generative AI assistant (Tencent Hunyuan, model identifier Hy, developed by Tencent and accessed through the WorkBuddy "
 "desktop application; the interface did not expose a numeric model version) was used for language editing, drafting support and "
 "formatting of tables and figure legends. A second generative AI assistant (ChatGPT, OpenAI; the exact model version was not "
 "identified in the responses) was used for critical reading of the manuscript, methodological commentary and suggestions for "
 "refinement of selected text. Neither tool was used to generate data, statistical results, images or references. All analyses "
 "were performed by the authors with the software specified above, and every AI-assisted passage was verified against the "
 "underlying analysis outputs and cited literature. The authors take full responsibility for the content; no generative AI tool "
 "is listed as an author.")

AI_A = ("The authors used a generative artificial intelligence assistant (Tencent Hunyuan, model identifier Hy, developed by Tencent "
 "and accessed through the WorkBuddy desktop application, workbuddy.cn; the interface did not expose a numeric model version; "
 "accessed between 2 and 9 November 2026) for language editing, drafting of the manuscript text and formatting of tables and "
 "figure legends. A second generative AI assistant (ChatGPT, OpenAI; exact model version not identified in the responses; "
 "accessed during the November 2026 revision) provided critical review, methodological commentary and suggestions for refinement "
 "of selected sections. Neither tool was used to generate data, statistical results, images or references; all analyses were "
 "performed by the authors with the software specified in Methods; every AI-assisted passage was verified against the underlying "
 "analysis outputs and the cited literature; and the authors take full responsibility for the accuracy, originality and integrity "
 "of the content. No generative AI tool is listed as an author.")

rep("A generative AI assistant (Tencent Hunyuan, model identifier Hy, developed by Tencent and accessed through the WorkBuddy desktop application; the interface did not expose a numeric model version) was used for language editing, drafting support and formatting of tables and figure legends. The tool was not used to generate data, statistical results, images or references. All analyses were performed by the authors with the software specified above, and every AI-assisted passage was verified against the underlying analysis outputs and cited literature. No generative AI tool is listed as an author.", AI_M, "M-AI Methods")
rep("The authors used a generative artificial intelligence assistant (Tencent Hunyuan, model identifier Hy, developed by Tencent and accessed through the WorkBuddy desktop application, workbuddy.cn; the interface did not expose a numeric model version; accessed between 2 and 9 November 2026) to assist with language editing, drafting of the manuscript text and formatting of tables and figure legends; the tool was not used to generate data, statistical results, images or references, all analyses were performed by the authors with the software specified in Methods, every AI-assisted passage was verified against the underlying analysis outputs and the cited literature, and the authors take full responsibility for the accuracy, originality and integrity of the content.", AI_A, "M-AI Acknowledgements")
io.open(os.path.join(BASE,"11_manuscript","manuscript_BMC_v8.md"),"w",encoding="utf-8",newline="\n").write(s)
for x in log: print("  [OK]",x)

# ---------- build docx using 145 chain ----------
src145=io.open(os.path.join(BASE,"13_refcheck","145_build_v5_bmcgenomics.py"),encoding="utf-8").read()
t=src145
t=t.replace('V43 = os.path.join(SUB, "v4_3_do")','V43 = os.path.join(SUB, "v7_bmcgenomics")')
t=t.replace('V5 = os.path.join(SUB, "v5_bmcgenomics")','V5 = os.path.join(SUB, "v8_bmcgenomics")')
t=t.replace('"manuscript_BMC_v5.md"','"manuscript_BMC_v8.md"')
t=t.replace('"v5_bmcgenomics", "ESCC_BMC_Genomics_manuscript_v5.docx"','"v8_bmcgenomics", "ESCC_BMC_Genomics_manuscript_v8.docx"')
t=t.replace('DOCX = os.path.join(V5, "ESCC_BMC_Genomics_manuscript_v5.docx")','DOCX = os.path.join(V5, "ESCC_BMC_Genomics_manuscript_v8.docx")')
t=t.replace('Figure_legends_BMC_v5.txt','Figure_legends_BMC_v8.txt')
t=t.replace('Cover_letter_BMC_v5.txt','Cover_letter_BMC_v8.txt')
t=t.replace('_tmp_103_v5.py','_tmp_103_v8.py').replace('_tmp_122_v5.py','_tmp_122_v8.py')
marker="    ('FONT = \"Times New Roman\"', 'FONT = \"Times New Roman\"'),"
t=t.replace(marker,"    ('\"## Introduction\", 1', '\"## Background\", 1'),\n    ('\"## Data availability statement\", 1', '\"## Declarations\", 1'),\n"+marker)
a=t.find('# ---------------------------------------------------------------- 1 源稿 5 项'); b=t.find('# ---------------------------------------------------------------- 2 主稿 docx')
t=t[:a]+'# source md prepared by 155\n'+t[b:]
c=t.find('# ---------------------------------------------------------------- 4 cover letter')
if c!=-1: t=t[:c]+'\n# cover letter handled by 155\n'
sp=os.path.join(BASE,"13_refcheck","_tmp_145_v8.py"); io.open(sp,"w",encoding="utf-8").write(t)
exec(compile(t,sp,"exec"),{"__name__":"__main__"}); os.remove(sp)
io.open(os.path.join(BASE,"11_manuscript","manuscript_BMC_v8.md"),"w",encoding="utf-8",newline="\n").write(s)
# table layout
from docx import Document
from docx.shared import Cm
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
DOCX=os.path.join(V8,"ESCC_BMC_Genomics_manuscript_v8.docx")
d=Document(DOCX); W=[Cm(2.7),Cm(2.5),Cm(2.9),Cm(2.1),Cm(4.6)]
for p in d.paragraphs:
    if p.text.strip().startswith("Table 2."): p.paragraph_format.keep_with_next=True
t2=d.tables[1]; t2.autofit=False
for row in t2.rows:
    for j,cc in enumerate(row.cells):
        if j<len(W): cc.width=W[j]
for tb in d.tables:
    for ri,row in enumerate(tb.rows):
        trPr=row._tr.get_or_add_trPr()
        if trPr.find(qn("w:cantSplit")) is None: trPr.append(OxmlElement("w:cantSplit"))
        if ri==0 and trPr.find(qn("w:tblHeader")) is None:
            h=OxmlElement("w:tblHeader"); h.set(qn("w:val"),"true"); trPr.append(h)
d.save(DOCX)
print("[OK] docx built:",os.path.basename(DOCX))

# ---------- v1.1.0 code bundle ----------
OLD=r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo\12_submission\code_release\escc-endo-code-v1.0.0"
NEWDIR=os.path.join(SUB,"code_release","escc-endo-code-v1.1.0")
if os.path.exists(NEWDIR): shutil.rmtree(NEWDIR)
shutil.copytree(OLD,NEWDIR)
# add reproducibility scripts for v4-v7 analyses
add=[
 "140_voom_edgeR_sensitivity.R","141_nested_cv_inner_reselect.py",
 "147_go_tested_universe.R","150_refine_v6_1.py","152_build_v7_final.py","155_v8_ai_and_codebundle.py",
]
for f in add:
    src=os.path.join(BASE,"13_refcheck",f)
    assert os.path.exists(src),f
    shutil.copy2(src,os.path.join(NEWDIR,"13_refcheck",f))
# metadata updates
rd=io.open(os.path.join(NEWDIR,"README.md"),encoding="utf-8").read()
rd=rd.replace("v1.0.0","v1.1.0")
rd += "\n\n## v1.1.0 additions (2026-11-25)\n\n"
rd += "- `13_refcheck/140_voom_edgeR_sensitivity.R`: TMM+voom pseudobulk sensitivity (tip signature 487 -> 241 genes).\n"
rd += "- `13_refcheck/141_nested_cv_inner_reselect.py`: nested CV with feature re-selection in every inner training split.\n"
rd += "- `13_refcheck/147_go_tested_universe.R`: GO enrichment with tested-gene universe (11,919 SYMBOL -> 10,836 Entrez -> 9,889 GO-BP background); produces Supplementary Table S16. Run with R 4.4.3 (clusterProfiler 4.14.6, org.Hs.eg.db) on `D:/escc_sc/out/endo_annotated.rds` (public GSE160269-derived object).\n"
rd += "- `13_refcheck/150_refine_v6_1.py`, `152_build_v7_final.py`, `155_v8_ai_and_codebundle.py`: manuscript/table building and reporting-chain scripts.\n"
rd += "- This release corresponds to manuscript v8 / BMC Genomics and Supplementary Tables S1-S16.\n"
io.open(os.path.join(NEWDIR,"README.md"),"w",encoding="utf-8",newline="\n").write(rd)
cff=io.open(os.path.join(NEWDIR,"CITATION.cff"),encoding="utf-8").read().replace("v1.0.0","v1.1.0").replace("24096140","PENDING_NEW_VERSION_DOI")
io.open(os.path.join(NEWDIR,"CITATION.cff"),"w",encoding="utf-8",newline="\n").write(cff)
# regenerate manifest
man=os.path.join(NEWDIR,"MANIFEST.tsv")
with io.open(man,"w",encoding="utf-8",newline="\n") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["path","sha256"])
    for r,_,fs in os.walk(NEWDIR):
        for f in sorted(fs):
            p=os.path.join(r,f); rel=os.path.relpath(p,NEWDIR).replace(os.sep,"/")
            if rel=="MANIFEST.tsv": continue
            w.writerow([rel,hashlib.sha256(open(p,"rb").read()).hexdigest()])
zp=os.path.join(SUB,"code_release","escc-endo-code-v1.1.0.zip")
with zipfile.ZipFile(zp,"w",zipfile.ZIP_DEFLATED) as z:
    for r,_,fs in os.walk(NEWDIR):
        for f in sorted(fs):
            p=os.path.join(r,f); z.write(p,os.path.join("escc-endo-code-v1.1.0",os.path.relpath(p,NEWDIR)))
print("[OK] code bundle v1.1.0:",len([1 for r,_,fs in os.walk(NEWDIR) for f in fs]),"files")
# cleanup stale names in V8
for f in ["Cover_letter_DO_v4_3.docx","Cover_letter_DO_v4_3.txt","ESCC_DO_manuscript_v4_3.docx",
          "Figure_legends_DO_v4_3.txt","投稿操作要点_DO_v4_3.md","Cover_letter_BMC_v7.txt",
          "ESCC_BMC_Genomics_manuscript_v7.docx","Figure_legends_BMC_v7.txt","投稿要点_BMC_v7.md",
          "校验报告_v7.md","整改记录_v7.md"]:
    p=os.path.join(V8,f)
    if os.path.exists(p): os.remove(p)
V7P=os.path.join(SUB,"v7_bmcgenomics")
cov=io.open(os.path.join(V7P,"Cover_letter_BMC_v7.txt"),encoding="utf-8").read().replace("v7","v8")
io.open(os.path.join(V8,"Cover_letter_BMC_v8.txt"),"w",encoding="utf-8",newline="\n").write(cov)
pts=io.open(os.path.join(V7P,"投稿要点_BMC_v7.md"),encoding="utf-8").read().replace("v7","v8")
io.open(os.path.join(V8,"投稿要点_BMC_v8.md"),"w",encoding="utf-8",newline="\n").write(pts)
print("[OK] FINAL docx + code bundle ready")
