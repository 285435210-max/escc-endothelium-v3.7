# -*- coding: utf-8 -*-
"""v6_1 -> v7_bmcgenomics FINAL submission version.
External-review targeted revision: 6 must-fix + all minors.
NO new analyses. All changes are text/structure/metadata.
"""
import io, os, re, shutil
BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE, "12_submission")
V61 = os.path.join(SUB, "v6_1_bmcgenomics")
V7 = os.path.join(SUB, "v7_bmcgenomics")
os.makedirs(V7, exist_ok=True)
for f in os.listdir(V61):
    p = os.path.join(V61, f)
    if os.path.isfile(p): shutil.copy2(p, os.path.join(V7, f))
for sub in ("figures_tiff","figures_png","figures_pdf"):
    src,dst = os.path.join(V61,sub), os.path.join(V7,sub)
    if os.path.exists(src) and not os.path.exists(dst): shutil.copytree(src,dst)
print("[OK] copied v6_1 -> v7")

s = io.open(os.path.join(BASE,"11_manuscript","manuscript_BMC_v6_1.md"),encoding="utf-8").read()
log=[]
def rep(old,new,tag,expect=1):
    global s
    n=s.count(old)
    if n!=expect: raise SystemExit("[FAIL] %s expect %d found %d :: %r"%(tag,expect,n,old[:150]))
    s=s.replace(old,new); log.append("%-58s x%d"%(tag,n))

# ============ MUST-FIX 1: GO cross-background wording ============
rep("The six leading GO biological-process terms were concordant between the two signature definitions and between the two background choices, so the functional interpretation does not depend on this choice.",
    "Within the tested-gene universe, the six leading GO biological-process terms were concordant between the 487-gene and 241-gene signatures, with a rank exchange between angiogenesis and blood-vessel development. Across the genome-wide and tested-gene backgrounds, only the three leading extracellular-matrix and extracellular-structure organization terms were retained; subsequent terms differed. The_input_SYMBOL pool comprised 11,919 genes; 10,836 unique Entrez identifiers were obtained after Symbol-to-Entrez conversion and deduplication (90.9% unique-Entrez yield relative to input symbols), of which 9,889 carried a GO biological-process annotation and formed the effective background; 422 and 214 input genes entered enrichment for the 487-gene and 241-gene signatures, respectively (Table S16).",
    "M1a Methods GO cross-background precise")
rep("The six leading terms were concordant, but this does not imply that every enrichment term was identical between backgrounds or signatures.",
    "The six leading terms were concordant between the two signatures within the tested-gene universe, but only the first three extracellular-matrix and extracellular-structure terms were shared with the genome-wide background; secondary terms and their rankings differed.",
    "M1b Results GO cross-background precise")
rep("a tested-gene-universe sensitivity analysis gave identical leading terms but is itself restricted to genes passing the expression filter (Table S16).",
    "a tested-gene-universe sensitivity analysis retained the leading extracellular-matrix and extracellular-structure organization terms, although secondary terms and their rankings differed from those obtained with the genome-wide background (Table S16).",
    "M1c Limitations GO wording")
# Figure 3B legend: specify genome-wide background
rep("**(B)** Gene Ontology biological process enrichment of the 487-gene tip signature; dot size indicates gene count and color indicates -log10 adjusted P (brick red is more significant)",
    "**(B)** Gene Ontology biological process enrichment of the 487-gene tip signature using the default genome-wide org.Hs.eg.db background; dot size indicates gene count and color indicates -log10 adjusted P (brick red is more significant). Table S16 reports the corresponding analysis with a tested-gene universe.",
    "M1d Fig3B legend background")
# ECM terms as one theme
rep("The first three terms reflect one overlapping extracellular-matrix organization theme rather than three independent biological mechanisms.",
    "The first three terms reflect one overlapping extracellular-matrix organization theme rather than three independent biological mechanisms.") if "The first three terms reflect one overlapping" in s else None
rep("This is not simply a proliferation signature.",
    "These three terms share largely overlapping gene sets and are interpreted as one extracellular-matrix organization theme rather than three independent mechanisms; this is not simply a proliferation signature.",
    "M1e ECM theme single interpretation")

# ============ MUST-FIX 2: Data availability S1-S15 -> S1-S16 ============
rep("Tables S1-S15, Figures S1-S4", "Tables S1-S16, Figures S1-S4", "M2a DA table range 1")
rep("Supplementary Tables S1-S15 and Supplementary Figures S1-S4", "Supplementary Tables S1-S16 and Supplementary Figures S1-S4", "M2b DA table range 2")
rep("Table S15 documents the pseudobulk and nested cross-validation sensitivity analyses.",
    "Table S15 documents the pseudobulk and nested cross-validation sensitivity analyses, and Table S16 the tested-gene-universe GO enrichment sensitivity analysis.",
    "M2c DA add S16")

# ============ MUST-FIX 4: doublet accurate GEO statement ============
rep("The analysis was performed on a pre-extracted CD45-negative endothelial-cell matrix. Because the complete, unfiltered CD45-negative cell population was unavailable, a comprehensive sample-level doublet-detection procedure could not be performed; potential contamination was assessed indirectly through quality-control metrics, canonical lineage markers and exclusion of non-endothelial clusters.",
    "The analysis used the pre-extracted endothelial-cell UMI matrix supplied with GSE160269. A broader processed CD45-negative UMI matrix and cell-level metadata are publicly available from GEO, and raw reads are deposited in the GSA repository under HRA000195. We did not reprocess the complete CD45-negative compartment or perform a new comprehensive sample-level doublet-detection analysis; the analysis therefore relied on the supplied endothelial-cell subset, quality-control metrics, canonical lineage-marker assessment and exclusion of non-endothelial clusters.",
    "M4 doublet accurate GEO statement")

# ============ MUST-FIX 5: arterial P + inflamed CI ============
rep("Arterial endothelium was unchanged (12.2% versus 13.4%, P = 0.67).",
    "Arterial endothelium was unchanged in the primary comparison (median 12.1% in adjacent mucosa versus 13.4% in tumors; tie-aware permutation P = 0.616, BH q = 0.616).",
    "M5a arterial P 0.616")
rep("+4.6 for the interferon-responsive state (3.2 to 6.4).",
    "+4.6 for the interferon-responsive state (3.2 to 6.4) and +9.9 for the inflamed venule state (8.3 to 11.5).",
    "M5b inflamed venule CI")
# Limitations downsampling add inflamed venule
rep("venous, tip-like and matrix-remodeling states passed correction in 100% of down-sampling replicates at both cell thresholds",
    "venous, matrix-remodeling and inflamed-venule states passed correction in 100% of down-sampling replicates at both cell thresholds, and tip-like endothelium in 99.5% at 20 cells and 100% at 30 cells",
    "M5c Limitations downsampling add inflamed") if "venous, tip-like and matrix-remodeling states passed correction in 100% of down-sampling replicates at both cell thresholds" in s else None

# ============ MINOR writing ============
rep("Tip cells carry an extracellular-matrix program",
    "The tip-like endothelial state carries an extracellular-matrix program", "m1 tip-like state")
rep("under a comparison that respects the paired structure of the libraries.",
    "under a primary comparison that excluded the four tumor libraries sharing patients with adjacent samples, thereby avoiding treatment of overlapping paired libraries as independent observations.",
    "m2 paired design wording")
rep("Tumor endothelium in esophageal squamous cell carcinoma (ESCC) is usually treated as a homogeneous population, and its endothelial-state composition has not been systematically re-evaluated at finer transcriptional resolution with explicit assessment of sample imbalance and compositional constraints.",
    "Although previous esophageal squamous cell carcinoma (ESCC) single-cell atlases have described endothelial heterogeneity, the robustness of finer endothelial-state composition differences has not been systematically evaluated with explicit assessment of sample imbalance and compositional constraints.",
    "m3 abstract opening")
# Introduction duplicated purpose
rep(" We therefore re-analysed the CD45-negative endothelial compartment with the venous, tip-like, matrix-remodeling and interferon-responsive programs separated, and assessed how the four adjacent libraries, the paired design and the closure constraint affect the estimated differences. We therefore re-examined composition with a refined transcriptional-state framework, quantified robustness to sample imbalance and compositional bias, and tested the specificity and clinical limitations of a tip-derived expression score.",
    " We therefore re-examined composition with a refined transcriptional-state framework separating venous, tip-like, matrix-remodeling and interferon-responsive programs, quantified robustness to the four adjacent libraries, the paired design and compositional closure, and tested the specificity and clinical limitations of a tip-derived expression score.",
    "m4 intro duplicate purpose")
# inflamed venule spelling unify
s = s.replace("inflamed-venule", "inflamed venule")
log.append("m5 inflamed venule hyphen unified x%d" % s.count("inflamed venule"))

# ============ MUST-FIX 6: BMC structure ============
# Abstract: split into Background / Results / Conclusions
i = s.find("## Abstract\n"); j = s.find("## Introduction", i)
assert i>0 and j>i
ab = s[i+len("## Abstract\n"):j].strip()
# find sentence boundaries
m1 = re.search(r"(?<=[.])\s+(?=We addressed this gap)", ab)
assert m1, "abstract split1 not found"
background = ab[:m1.start()].strip()
rest = ab[m1.end():].strip()
m2 = re.search(r"(?<=[.])\s+(?=As an analysis of existing data)", rest)
assert m2, "abstract split2 not found"
results = rest[:m2.start()].strip()
conclusions = rest[m2.end():].strip()
new_ab = "## Abstract\n\n**Background:** " + background + "\n\n**Results:** " + results + "\n\n**Conclusions:** " + conclusions + "\n\n"
s = s[:i] + new_ab + s[j:]
log.append("M6a Abstract split Background/Results/Conclusions")

# Section renames
s = s.replace("## Introduction\n", "## Background\n", 1)
s = s.replace("## Materials and methods\n", "## Methods\n", 1)
log.append("M6b Introduction->Background, Materials and methods->Methods")

# Author contributions initials
rep("Bin Liu conceived the study and performed the analyses; Bin Liu, Zhen Ni and Runbin Liang interpreted the data; Bin Liu drafted the manuscript; Zhen Ni and Runbin Liang revised it critically. All authors read and approved the final version.",
    "BL conceived the study and performed the analyses; BL, ZN and RL interpreted the data; BL drafted the manuscript; ZN and RL critically revised the manuscript. All authors read and approved the final version.",
    "M6c Author contributions initials")

# Methods AI disclosure: append a short subsection in Methods (before Results)
i = s.find("## Results")
ai_methods = ("### Generative AI use\n\n"
              "A generative AI assistant (Tencent Hunyuan, model identifier Hy, developed by Tencent and accessed through the WorkBuddy desktop application; the interface did not expose a numeric model version) was used for language editing, drafting support and formatting of tables and figure legends. The tool was not used to generate data, statistical results, images or references. All analyses were performed by the authors with the software specified above, and every AI-assisted passage was verified against the underlying analysis outputs and cited literature. No generative AI tool is listed as an author.\n\n")
s = s[:i] + ai_methods + s[i:]
log.append("M6d Methods generative-AI subsection")

# Declarations: rename sections and add umbrella heading before Ethics
rep("## Ethics statement\n\nThis study used publicly available, de-identified data; no ethical approval or informed consent was required.",
    "## Declarations\n\n### Ethics approval and consent to participate\n\nThis study used publicly available, de-identified data; no ethical approval or informed consent was required.",
    "M6e Declarations umbrella + ethics heading")
rep("## Data availability statement", "### Availability of data and materials", "M6f DA heading")
rep("## Consent for publication", "### Consent for publication", "M6g consent heading")
rep("## Competing interests", "### Competing interests", "M6h competing heading")
rep("## Funding\n\nNone.", "### Funding\n\nNone.", "M6i funding heading")
rep("## Author contributions", "### Authors' contributions", "M6j authors heading")
rep("## Acknowledgements", "### Acknowledgements", "M6k ack heading")
# Abbreviations stay outside declarations as a normal section heading before Declarations
rep("## Abbreviations", "## Abbreviations", "M6l abbreviations kept", expect=1) if False else None

# Software/OS note in Methods
rep("scikit-learn 1.9.1, statsmodels 0.15.0, lifelines 0.30.3, matplotlib 3.11.2).",
    "scikit-learn 1.9.1, statsmodels 0.15.0, lifelines 0.30.3, matplotlib 3.11.2), run on Windows 11.",
    "m12 OS note")

# GEO access date + full Zenodo link
n=s.count("Zenodo under the version DOI 10.5281/zenodo.24096140")
s=s.replace("Zenodo under the version DOI 10.5281/zenodo.24096140","Zenodo under the version DOI https://doi.org/10.5281/zenodo.24096140")
log.append("m14 Zenodo full link x%d"%n)

io.open(os.path.join(BASE,"11_manuscript","manuscript_BMC_v7.md"),"w",encoding="utf-8",newline="\n").write(s)
for x in log: print("  [OK]",x)
print("[OK] source saved %d chars"%len(s))
