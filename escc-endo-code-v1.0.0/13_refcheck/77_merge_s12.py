# -*- coding: utf-8 -*-
"""v3.5 supplementary tables: add S12, refresh S7/S8 CellChat outputs, fix README.

Follows the second external review:
  * primary comparison changed to 56 unmatched tumours vs 4 adjacent libraries;
  * CellChat re-run with an explicit Other/unclassified T group;
  * ligand-receptor results downgraded to a pooled transcript-level ranking.

The historical workbook filename is retained for script compatibility; a clean
v3.5 delivery copy is written under a new name.
"""
import os, shutil
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment

BASE = r"D:\龙虾文档\2026-11-02-09-16-41\ESCC_endo"
SUB = os.path.join(BASE, "12_submission")
WB = os.path.join(SUB, "Supplementary_Tables_S1-S10.xlsx")
BAK = os.path.join(SUB, "Supplementary_Tables_S1-S10_pre_S12_20261109.xlsx")
DELIVERY = os.path.join(SUB, "Supplementary_Tables_S1-S12_v3_5.xlsx")
S12 = os.path.join(BASE, "13_refcheck", "out_stats_v35")
CC = r"D:\escc_sc\cellchat_out"
ROB = os.path.join(CC, "robust_v35")

if not os.path.exists(BAK):
    shutil.copy2(WB, BAK)
    print("[OK] backup:", os.path.basename(BAK))

before = openpyxl.load_workbook(BAK, read_only=True, data_only=True)
old_names = list(before.sheetnames)
assert len(old_names) == 57, "unexpected pre-S12 sheet count %d" % len(old_names)
wb = openpyxl.load_workbook(WB)

CHANGED = set()


def native(v):
    if v is None:
        return ""
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, float) and pd.isna(v):
        return ""
    return v


def csv_rows(path):
    df = pd.read_csv(path)
    rows = [[native(c) for c in df.columns]]
    arr = df.astype(object).where(pd.notna(df), None).to_numpy().tolist()
    for rec in arr:
        rows.append([native(v) for v in rec])
    return rows


def write_sheet(name, rows, index=None):
    """Replace a sheet in place (keeping its position) or append a new one."""
    if name in wb.sheetnames:
        idx = wb.sheetnames.index(name)
        del wb[name]
    else:
        idx = index
    ws = wb.create_sheet(name, idx)
    for r in rows:
        ws.append(r)
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"
    ws.column_dimensions["A"].width = 34
    ws.column_dimensions["B"].width = 70
    CHANGED.add(name)
    print("[OK] sheet %-38s rows=%d" % (name, len(rows)))
    return ws


# ------------------------------------------------------------------ S12
S12_README = [
    ("Item", "Description"),
    ("Scope", "Supplementary Table S12 documents the v3.5 re-analysis of sample pairing "
              "and the exact-test family requested in the second external review "
              "(2026-11-09). No new data were generated; all inputs are the public "
              "GSE160269 composition table already used in S1-S11."),
    ("Design", "GSE160269 has 60 tumour and 4 adjacent non-malignant libraries. The four "
               "adjacent libraries (P126N, P127N, P128N, P130N) come from four patients "
               "who also contributed tumours P126T, P127T, P128T and P130T. Pooling all "
               "60 tumours with the 4 adjacent libraries therefore violates independence. "
               "The primary comparison uses the 56 unmatched tumours and the 4 adjacent "
               "libraries; the 4 pairs are reported as directions only (no P value), and "
               "the pooled 60-versus-4 comparison is retained as a sensitivity analysis."),
    ("Primary test", "Exact two-sided Mann-Whitney test on closed state fractions "
                     "(frac = n / sum of nine retained states). Benjamini-Hochberg "
                     "correction is applied within one pre-defined family of nine "
                     "endothelial states. Effect sizes are median differences, "
                     "Hodges-Lehmann two-sample shifts and 2,000-replicate patient-level "
                     "bootstrap 95% confidence intervals, all in percentage points."),
    ("Leave-one-out family", "Four prespecified states x four deletions of one adjacent "
                             "library = 16 exploratory comparisons. BH is reported both "
                             "jointly across all 16 tests and separately within each "
                             "deletion (four states). All 16 directions are preserved; "
                             "joint BH q = 0.00098-0.00252."),
    ("Down-sampling", "Equal cell counts of 20 and 30 endothelial cells per sample, 200 "
                      "replicates each. The candidate interferon-responsive state survives "
                      "BH correction in 16.0% and 85.5% of replicates at 20 and 30 cells, "
                      "respectively, whereas venous, tip-like and matrix-remodeling "
                      "differences are retained in 99.5-100% of replicates."),
    ("Display names", "Machine labels in the data columns map to the manuscript terms: "
                      "Venous -> Venous; Venous_inflamed -> Inflamed venule; "
                      "Tip_angiogenic -> Tip-like; Matrix_EC -> Matrix-remodeling; "
                      "IFNgamma_EC -> Interferon-responsive; Arterial, Capillary, "
                      "Proliferating and Lymphatic are unchanged. Machine labels are kept "
                      "in data tables for reproducibility; figures and text use display names."),
    ("S12_sample_pairing", "All 64 libraries with patient ID, tissue type, paired-tumour "
                           "flag, endothelial cell count and primary-analysis inclusion."),
    ("S12_primary_56v4", "Primary comparison: 56 unmatched tumours vs 4 adjacent libraries, "
                         "exact and asymptotic P values, BH q values, effect sizes and "
                         "bootstrap CIs for all nine states."),
    ("S12_paired_4pairs", "The four matched pairs, nine states each. Differences are "
                          "directions only; no inferential test is applied to four pairs "
                          "(minimum attainable two-sided sign-test P = 0.125)."),
    ("S12_sensitivity_60v4", "Original pooled 60-versus-4 design recomputed with the exact "
                             "test. Capillary endothelium crosses the 0.05 boundary "
                             "(q = 0.049); this state is not a stable result."),
    ("S12_loo_bh_family", "Leave-one-adjacent-normal-out results with both BH families."),
    ("S12_exact_vs_asymptotic", "Exact versus asymptotic Mann-Whitney P values and BH q "
                                "values for the primary design."),
    ("S12_downsampling", "FDR-pass rates under equal cell counts for both designs."),
    ("S12_analysis_code_inventory", "Paths, roles and outputs of the R and Python scripts "
                                    "used for the analyses reported here, so that a reviewer "
                                    "can map every number to a script. Scripts are available "
                                    "from the corresponding author on reasonable request; "
                                    "deposition in a public repository with a persistent "
                                    "identifier is an author action still pending."),
    ("Relation to S11", "S11 reports legacy sensitivity analyses computed under the pooled "
                        "60-versus-4 design with asymptotic tests. It remains valid as a "
                        "robustness check but is not the primary analysis; the v3.5 primary "
                        "comparison is S12_primary_56v4. S11 and S12 share the same closed "
                        "composition denominator."),
    ("Relation to S7/S8", "S7 and the CellChat sheets in S8 were regenerated from the v3.5 "
                          "re-run, in which marker-negative T cells are labelled "
                          "Other/unclassified T instead of being assigned to CD4 T."),
    ("Script", "13_refcheck/71_stats_v3_5.py (seed 20261109); scripts_r/cellchat_run_v3_5.R "
               "and scripts_r/cellchat_robust_v3_5.R. Machine-readable summary: "
               "13_refcheck/out_stats_v35/S12_summary.json."),
    ("Interpretation", "The venous, tip-like and matrix-remodeling differences survive the "
                       "independence-preserving design, exact testing, deletion of any "
                       "adjacent library and equal-count down-sampling. The interferon-"
                       "responsive state is a candidate/exploratory finding with "
                       "inconsistent paired directions and low down-sampling stability."),
]

write_sheet("S12_README", S12_README)

S12_FILES = [
    ("S12_sample_pairing", "S12_sample_pairing.csv"),
    ("S12_primary_56v4", "S12_primary_56v4.csv"),
    ("S12_paired_4pairs", "S12_paired_4pairs.csv"),
    ("S12_sensitivity_60v4", "S12_sensitivity_60v4.csv"),
    ("S12_loo_bh_family", "S12_loo_bh_family.csv"),
    ("S12_exact_vs_asymptotic", "S12_exact_vs_asymptotic.csv"),
    ("S12_downsampling", "S12_downsampling.csv"),
]
for sheet, fn in S12_FILES:
    write_sheet(sheet, csv_rows(os.path.join(S12, fn)))

PENDING = "Not yet deposited; author action pending"
CODE_INVENTORY = [("Category", "Script", "Role", "Output", "Repository status")]
for cat, script, role, out in [
    ("scRNA preprocessing (R)", "scripts_r/endo_cluster.R",
     "QC, normalization, harmony integration, clustering, UMAP", "08_scRNA/endothelial/"),
    ("scRNA preprocessing (R)", "scripts_r/endo_annotate.R",
     "Nine-state marker-based annotation and contaminant removal", "celltype_proportion_by_sample.csv"),
    ("scRNA preprocessing (R)", "scripts_r/endo_de.R",
     "State marker genes and patient-paired pseudobulk limma comparisons", "S2-S6"),
    ("scRNA preprocessing (R)", "scripts_r/endo_ligands.R",
     "Ligand expression ranking across endothelial states", "S8_ligand_rank_*"),
    ("scRNA preprocessing (R)", "scripts_r/42_export_reference.R",
     "Export reference labels and centroids for label transfer", "S9"),
    ("scRNA preprocessing (R)", "scripts_r/annot_check.R",
     "SingleR/celldex reference-based annotation check", "Figure S1"),
    ("Composition and statistics (Python)", "scripts/34_revise_composition.py",
     "Closed composition, BH-corrected group tests, CLR and Dirichlet-multinomial nulls",
     "S8_group_tests_BH, S8_compositional_null, S8_venous_tip_CLR"),
    ("Composition and statistics (Python)", "13_refcheck/71_stats_v3_5.py",
     "v3.5 independence-preserving 56-vs-4 exact tests, paired directions, "
     "60-vs-4 sensitivity, LOO BH families, bootstrap CIs, down-sampling", "S12"),
    ("Composition and statistics (Python)", "13_refcheck/63_sensitivity_v3_4.py",
     "Legacy pooled-design sensitivity analyses (S11)", "S11"),
    ("Diagnostic panel (Python)", "scripts/32_panel_nested_cv.py",
     "Nested patient-grouped CV with within-fold DEG and inner-fold L1 selection",
     "S8_nestedCV_folds, S8_gene_stability, S8_random_controls"),
    ("Diagnostic panel (Python)", "scripts/30_gse225178_validate.py",
     "Non-refitted descriptive checks in GSE20347, GSE45670 and GSE225178",
     "S8_external_validation, Figure S3"),
    ("Prognosis (Python)", "scripts/38_ph_and_stage.py",
     "Cox models, Schoenfeld test, stage association", "S8_stage_association_BH, S8_PH_test"),
    ("Prognosis (Python)", "scripts/48_tcga_esca_prognosis.py",
     "Cross-platform TCGA-ESCA ESCC prognostic re-check", "S10"),
    ("Correlation (Python)", "scripts/15_ec_t_correlation.py",
     "Nine-state versus T-cell-group Spearman correlations with BH correction",
     "S8_EC_vs_T_9x5_BH, S8_ExCD8_*"),
    ("External validation (Python)", "scripts/43_external_label_transfer.py",
     "Label transfer to GSE196756, GSE203067 and CELLxGENE atlas", "S9, Figure S4"),
    ("External validation (Python)", "scripts/47_paired_direction.py",
     "Two-sided paired sign tests across seven matched pairs", "S9"),
    ("CellChat (R)", "scripts_r/cellchat_run_v3_5.R",
     "v3.5 primary run with Other/unclassified T label; nboot = 50",
     "S7_CellChat_interactions, S7_tcell_group_counts_v3_5"),
    ("CellChat (R)", "scripts_r/cellchat_robust_v3_5.R",
     "v3.5 nboot = 200 repeat under two down-sampling seeds",
     "S8_CellChat_200perm, S8_CellChat_pathways_s1/s2"),
    ("CellChat (R)", "scripts_r/cellchat_run.R / cellchat_robust.R",
     "Archived v3.4 runs (marker-negative cells assigned to CD4 T)",
     "pre_S12 workbook backup"),
    ("Reproducibility (R)", "scripts_r/reproducibility_checks.R",
     "Software versions, seeds, cluster-resolution sensitivity", "S8_software, S8_seeds, S8_cluster_resolution"),
    ("Reproducibility (R)", "scripts_r/supp_tables.R",
     "Assembly of supplementary tables S1-S7", "S1-S7"),
    ("Figures (R)", "scripts_r/fig1_2.R, fig3_4.R, fig5_6.R, fig7.R, cellchat_figs.R",
     "Main figures 1-7 and Figure 5 with v3.5 display names and 56-vs-4 design",
     "12_submission/figures/"),
    ("Figures (Python)", "13_refcheck/72_patch_figure_scripts_v3_5.py, "
     "74_patch_figures_final.py, 73_replot_S4_v3_5.py",
     "v3.5 figure-script patches and Figure S4 replot", "figures/FigureS4.*"),
    ("Supplementary tables (Python)", "13_refcheck/77_merge_s12.py, "
     "78_split_and_verify.py",
     "S12 assembly, S7/S8 refresh, S1-S12 standalone files and consistency checks",
     "Supplementary_Tables_S1-S12_v3_5.xlsx, consistency report"),
    ("Manuscript (Python)", "13_refcheck/75_revise_en_v3_5.py, "
     "76_revise_zh_v3_5.py, 67/68/69_build_v3_4*.py",
     "English/Chinese source revisions and Frontiers-formatted docx generation",
     "manuscript_FI.md, manuscript_zh.md, three v3.5 docx files"),
]:
    CODE_INVENTORY.append((cat, script, role, out, PENDING))
write_sheet("S12_analysis_code_inventory", CODE_INVENTORY)

# ------------------------------------------------------------------ S7 re-run
S7_README = [
    ("Item", "Description"),
    ("Scope", "All ligand-receptor pairs retained after CellChat filtering in the v3.5 "
              "re-run (scripts_r/cellchat_run_v3_5.R), on cells pooled across tumours, "
              "adjacent mucosa and patients."),
    ("Why regenerated", "The v3.4 analysis assigned every marker-negative T cell to CD4 T "
                        "by default. In v3.5 these 2,109 cells form an explicit "
                        "Other/unclassified T group and are excluded from the "
                        "endothelial-to-T summaries. The v3.4 table (5,758 rows) is "
                        "archived in Supplementary_Tables_S1-S10_pre_S12_20261109.xlsx."),
    ("Columns", "prob is CellChat's dimensionless communication score (a summed triMean "
                "quantity in arbitrary units), NOT a probability. pval is the permutation "
                "P value from 50 permutations; with nboot = 50 its resolution is 0.02, so "
                "individual pairs cannot be declared statistically significant. All "
                "results are descriptive expression rankings, not evidence of "
                "cell-cell communication or antigen presentation."),
    ("S7_CellChat_interactions", "6,132 retained sender-receiver-resolved pairs; 1,556 "
                                 "connect an endothelial state to one of the five "
                                 "interpretable T-cell groups."),
    ("S7_tcell_group_counts_v3_5", "Mutually exclusive T-cell labels used in the CellChat "
                                   "object (8,000 sampled cells)."),
    ("Caution", "Tumour-derived and adjacent-derived cells from different patients were "
                "pooled, so a ligand and a receptor in one modelled pair need not originate "
                "from the same patient or the same tissue microenvironment. Do not "
                "interpret these rows as tumour-specific interactions."),
]
write_sheet("S7_README", S7_README, index=wb.sheetnames.index("S7_CellChat_interactions"))
write_sheet("S7_CellChat_interactions", csv_rows(os.path.join(CC, "cellchat_net.csv")))
write_sheet("S7_tcell_group_counts_v3_5",
            csv_rows(os.path.join(CC, "tcell_group_counts_v3_5.csv")))

# ------------------------------------------------------------------ S8 re-run
s1 = pd.read_csv(os.path.join(ROB, "summary_seed1.csv"))
s2 = pd.read_csv(os.path.join(ROB, "summary_seed2.csv"))
summ = pd.concat([s1, s2], ignore_index=True)
write_sheet("S8_CellChat_200perm", [list(summ.columns)] +
            [[native(v) for v in rec] for rec in
             summ.astype(object).where(pd.notna(summ), None).to_numpy().tolist()])

for seed in (1, 2):
    d = pd.read_csv(os.path.join(ROB, "pathway_ECtoT_seed%d.csv" % seed))
    d.insert(0, "seed", seed)
    write_sheet("S8_CellChat_pathways_s%d" % seed,
                [list(d.columns)] + [[native(v) for v in rec] for rec in
                                    d.astype(object).where(pd.notna(d), None).to_numpy().tolist()])

write_sheet("S8_tcell_groups_v3_5",
            csv_rows(os.path.join(ROB, "tcell_group_counts.csv")))

def set_readme_row(ws, key, value):
    """Idempotent: update an existing first-column key or append a new row."""
    for i in range(1, ws.max_row + 1):
        if ws.cell(row=i, column=1).value == key:
            ws.cell(row=i, column=2).value = value
            return
    ws.append([key, value])


set_readme_row(wb["S8_README"], "CellChat v3.5 re-run",
               "S8_CellChat_200perm, S8_CellChat_pathways_s1/s2 and S8_tcell_groups_v3_5 "
               "were regenerated with scripts_r/cellchat_robust_v3_5.R (nboot = 200, "
               "two down-sampling seeds). Marker-negative T cells are now labelled "
               "Other/unclassified T (2,109 cells) and excluded from endothelial-to-T "
               "summaries. Top three EC-to-T pathways (MHC-I, MIF, APP) are unchanged. "
               "The 1,556 EC-to-T pairs quoted in the manuscript come from the primary "
               "nboot = 50 run in S7, whereas these 200-permutation repeats give 1,558 "
               "and 1,551 pairs.")
CHANGED.add("S8_README")

set_readme_row(wb["S11_README"], "Note (v3.5)",
               "S11 was computed under the pooled 60-versus-4 design with asymptotic tests. "
               "It is retained as a robustness check. The v3.5 primary comparison respects "
               "sample independence (56 unmatched tumours vs 4 adjacent libraries, exact test) "
               "and is reported in Supplementary Table S12.")
CHANGED.add("S11_README")

# ------------------------------------------------------------------ top README
top = wb["README"]
updates = {
    "Title": "Supplementary Tables S1-S12",
    "File": "Supplementary_Tables_S1-S12_v3_5.xlsx (v3.5 delivery copy; identical to the "
            "historical Supplementary_Tables_S1-S10.xlsx workbook, whose filename is "
            "retained for script compatibility).",
    "S7": "Retained expressed ligand-receptor pairs from the v3.5 CellChat re-run "
          "(CellChat v2, CellChatDB.human, triMean, 50 permutations; "
          "Other/unclassified T excluded from endothelial-to-T summaries; see S7_README). "
          "Descriptive pooled transcript-level ranking only, not significant interactions.",
    "S11": "External-reviewer sensitivity analyses under the legacy pooled 60-versus-4 "
           "design: leave-one-normal-out, 529-cluster retention, exact tests, "
           "Hodges-Lehmann effect sizes and bootstrap CIs (7 sheets, S11_README first). "
           "Superseded as primary analysis by S12.",
}
for i in range(1, top.max_row + 1):
    key = top.cell(row=i, column=1).value
    if key in updates:
        top.cell(row=i, column=2).value = updates.pop(key)
assert not updates, "README keys not found: %s" % list(updates)
set_readme_row(top, "S12",
               "Independence-preserving sample-pairing re-analysis: exact primary "
               "56-unmatched-tumours-versus-4 comparison, four paired directions, "
               "pooled 60-versus-4 sensitivity, leave-one-out BH families, "
               "exact-versus-asymptotic contrast and down-sampling "
               "(9 sheets, S12_README first).")
set_readme_row(top, "Note on v3.5",
               "S7 and the CellChat sheets in S8 were regenerated after "
               "introducing the Other/unclassified T label; the archived "
               "pre-S12 workbook preserves the earlier outputs.")
CHANGED.add("README")

# Keep the workbook readable: README first, then S1..S12 blocks, preserving
# the original relative order inside each block.
def sort_key(name, pos):
    import re
    m = re.match(r"^S(\d+)", name)
    return (0, 0, pos) if name == "README" else ((1, int(m.group(1)), pos) if m else (2, 0, pos))

positions = {n: i for i, n in enumerate(wb.sheetnames)}
wb._sheets.sort(key=lambda ws: sort_key(ws.title, positions[ws.title]))
print("[OK] sheet order: %s ... %s" % (wb.sheetnames[0], wb.sheetnames[-1]))

wb.save(WB)
shutil.copy2(WB, DELIVERY)
print("[OK] saved:", WB)
print("[OK] delivery copy:", DELIVERY)

# ------------------------------------------------------------------ verification
after = openpyxl.load_workbook(DELIVERY, read_only=True, data_only=True)
print("[CHECK] sheets %d -> %d" % (len(old_names), len(after.sheetnames)))
assert len(after.sheetnames) == 69, "unexpected final sheet count %d" % len(after.sheetnames)
diffs = []
for name in old_names:
    if name in CHANGED:
        continue
    a = list(before[name].iter_rows(values_only=True))
    b = list(after[name].iter_rows(values_only=True))
    if a != b:
        diffs.append(name)
print("[CHECK] untouched original sheets verified: %d changed unexpectedly" % len(diffs))
assert not diffs, "unexpected cell changes: %s" % diffs

# key numbers
p = after["S12_primary_56v4"]
hdr = [c.value for c in next(p.iter_rows())]
rows = [dict(zip(hdr, r)) for r in p.iter_rows(min_row=2, values_only=True)]
by = {r["state"]: r for r in rows}
assert abs(by["Venous"]["med_normal_pct"] - 74.83766233766234) < 1e-9
assert abs(by["Venous"]["med_tumor_pct"] - 28.46065285114452) < 1e-9
assert abs(by["Venous"]["q_exact"] - 0.0004983235411732136) < 1e-12
assert by["Capillary"]["sig_exact"] in (False, "False", 0)
assert after["S7_CellChat_interactions"].max_row == 6133
assert after["S12_sample_pairing"].max_row == 65
print("[CHECK] key numbers PASSED")
before.close()
after.close()
print("[CHECK] workbook verification PASSED")
