# Analysis code for an ESCC endothelial-state atlas (manuscript v3.7)

**Software version:** 1.0.0
**Release date:** 2026-11-10
**Manuscript:** *Endothelial-state composition in esophageal squamous cell carcinoma: reduced venous representation, tip-like enrichment, and a candidate interferon-responsive state* (version 3.7, submitted to *Frontiers in Immunology*).
**Creators:** Liu Bin, Ni Zhen, Liang Runbin — Department of Thoracic Surgery,
Heyou Hospital, Foshan, Guangdong, China.
**Corresponding author ORCID:** 0000-0002-0617-3629.

## What is in this archive

Analysis source code only: 17 R scripts, 8 Python analysis scripts and 10 Python
scripts used to rebuild the supplementary tables and manuscript source. The
directory layout (`scripts_r/`, `scripts/`, `13_refcheck/`) preserves the paths
listed in Supplementary Table S12 (`S12_analysis_code_inventory`), which records
the 29 scripts used for the v3.5 re-analysis; the v3.6 and v3.7 additions are
listed under "Analysis order" below, so every number in the manuscript can be
traced to a script. `MANIFEST.tsv` gives the size and SHA-256 checksum of every
file. `PATHS.md` lists the machine-specific
absolute paths that must be edited before re-running the code.

## What is deliberately NOT included

- Raw or processed expression matrices (GSE160269, GSE53625, GSE20347,
  GSE45670, GSE225178, GSE196756, GSE203067, TCGA-ESCA). These are available
  from GEO (https://www.ncbi.nlm.nih.gov/geo/) and UCSC Xena
  (https://tcga.xenahubs.net) under their original accession numbers.
- Patient-level clinical tables beyond the public GEO supplementary files.
- Third-party full-text article PDFs, publisher material, figures and the
  manuscript files themselves. Redistribution rights were not granted.
- Intermediate objects (`.rds`, `.h5ad`, `.mtx.gz`), which are large and can be
  regenerated from the public sources.

## Public data sources

| Resource | Accession / location | Used for |
| --- | --- | --- |
| GSE160269 | GEO | 60 tumour + 4 adjacent non-malignant single-cell libraries |
| GSE53625, GSE20347, GSE45670 | GEO | bulk cohorts, diagnostic panel |
| GSE225178 | GEO | chemotherapy validation cohort |
| GSE196756, GSE203067 | GEO | external label transfer |
| CELLxGENE (Nowicki-Osuch gastrointestinal atlas) | CELLxGENE / API | external label transfer |
| TCGA-ESCA | UCSC Xena | cross-platform prognostic re-check |

## Software environment

R 4.4.3 with Seurat 5.5.0, harmony 1.2.4, limma 3.62.2, clusterProfiler 4.14.6,
CellChat 2.2.0.9001, ComplexHeatmap 2.22.0, SingleR 2.8.0 and celldex 1.16.0.
Python 3.13.14 with numpy 2.5.4, pandas 2.3.3, scipy 1.18.1, scikit-learn 1.9.1,
statsmodels 0.15.0, lifelines 0.30.3 and matplotlib 3.11.2. Exact sessionInfo
and the random-seed inventory are in Supplementary Table S8; the primary
independence-preserving re-analysis uses seed 20261109.

## Analysis order

1. `scripts_r/endo_cluster.R` -> `endo_annotate.R` -> `endo_de.R` ->
   `endo_ligands.R` -> `annot_check.R`: QC, integration, nine-state annotation,
   marker genes, reference-based check.
2. `scripts/34_revise_composition.py`: closed compositions, BH-corrected group
   tests, CLR and Dirichlet-multinomial nulls.
3. `scripts_r/cellchat_run_v3_5.R` and `cellchat_robust_v3_5.R`: primary
   (nboot = 50) and two-seed (nboot = 200) CellChat runs. Marker-negative T
   cells are labelled `Other/unclassified T`.
4. `13_refcheck/71_stats_v3_5.py`: the independence-preserving primary analysis — 56 unmatched
   tumours versus 4 adjacent non-malignant libraries, exact two-sided
   Mann-Whitney tests, BH over nine states, paired directions, pooled
   60-versus-4 sensitivity, leave-one-out BH families, bootstrap confidence
   intervals and equal-count down-sampling. Output reproduces Supplementary
   Table S12.
5. `scripts/32_panel_nested_cv.py`, `30_gse225178_validate.py`,
   `38_ph_and_stage.py`, `48_tcga_esca_prognosis.py`, `15_ec_t_correlation.py`,
   `43_external_label_transfer.py`, `47_paired_direction.py`: diagnostic,
   prognostic, correlation and external-validation analyses (S8-S10, Figure S4).
6. `scripts_r/fig*.R`, `cellchat_figs.R`, `reproducibility_checks.R`,
   `supp_tables.R`: figures and supplementary tables.
7. `13_refcheck/77_merge_s12.py`, `78_split_and_verify.py`, `94_merge_s13.py`,
   `101_merge_s14.py` and `107_verify_v3_7.py`: assemble the S1-S14 workbook,
   split stand-alone tables and run the numerical consistency checks against
   the manuscript (Tables S13 and S14 document tie-aware permutation tests,
   single-cell identity checks, seven-state stage tests and chemotherapy scores).
8. `13_refcheck/86_perm_v3_6.py`: complete tie-aware label-permutation enumeration
   (487,635 assignments for the primary design); `scripts_r/check_inflamed_identity.R`
   and `13_refcheck/99_v37_analyses.py`: Supplementary Table S14 analyses.

**Important:** the scripts contain absolute paths from the original Windows
analysis machine. They are research code, not a turn-key pipeline. Edit the
path constants listed in `PATHS.md` before re-running, and download the public
input data first. Every numerical result is also shipped as machine-readable
tables in Supplementary Tables S1-S14, so the manuscript can be audited
without re-running the full pipeline.

## Statistical design note

Four adjacent non-malignant libraries (P126N, P127N, P128N, P130N) share
patients with four tumours (P126T, P127T, P128T, P130T). Pooling all 60 tumours
would violate independence. The primary comparison therefore uses the 56
unmatched tumours; four matched pairs are reported as directions only; the
pooled 60-versus-4 comparison is a sensitivity analysis.

## Licence

Code: MIT (see `LICENSE`). Documentation in this archive: CC BY 4.0.
The underlying third-party datasets retain their original licences and are
not redistributed here.

## Citation

Once a DOI has been minted, cite this release as: Liu B, Ni Z, Liang R. (2026)
Analysis code for an ESCC endothelial-state atlas (manuscript v3.7). Zenodo.
Version 1.0.0. DOI.
