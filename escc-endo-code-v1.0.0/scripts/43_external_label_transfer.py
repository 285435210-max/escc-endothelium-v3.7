# -*- coding: utf-8 -*-
"""Step 1 (C-layer): external validation of the 9 endothelial states.

Stream-parses 10x mtx.gz (only ~198 target genes retained in memory, cached as
npz) + CELLxGENE h5ad, then evaluates three orthogonal readouts:
  1. author-supplied cell-type annotation of an independent normal EC atlas
  2. per-state marker expression / positive-cell fraction (no label transfer)
  3. exploratory label transfer by centroid rank correlation + module score

Outputs Table S9 sheets, Figure S4 and a run log.
"""
import os, gzip, json, warnings
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
BASE = r"D:/龙虾文档/2026-11-02-09-16-41/ESCC_endo"
RAW = os.path.join(BASE, "14_external_validation", "raw")
OUT = os.path.join(BASE, "14_external_validation", "out")
REF = r"D:/escc_sc/external"
G196 = r"D:/escc_sc/external/gse196756"
CACHE = r"D:/escc_sc/external/cache"
for d in [OUT, CACHE]:
    os.makedirs(d, exist_ok=True)

STATES = ["Arterial", "Capillary", "Venous", "Venous_inflamed",
          "Tip_angiogenic", "Matrix_EC", "IFNgamma_EC",
          "Proliferating", "Lymphatic"]
PAN_EC = ["PECAM1", "CDH5", "VWF", "CLDN5", "EGFL7", "KDR", "ERG",
          "FLT1", "ESAM", "CD34"]
LINEAGE = ["PTPRC", "EPCAM", "KRT5", "COL1A1", "DCN", "RGS5",
           "ACTA2", "PDGFRB", "CD79A", "LYZ"]

cent = pd.read_csv(os.path.join(REF, "ref_state_centroids.csv"))
mark = pd.read_csv(os.path.join(REF, "ref_state_markers_top20.csv"))
marker_by_state = {s: [g for g in mark.loc[mark.cluster == s, "gene"].astype(str)]
                   for s in STATES}
TARGET = sorted(set(PAN_EC + LINEAGE + [g for v in marker_by_state.values() for g in v]))
tpos = {g: i for i, g in enumerate(TARGET)}
cent_genes = set(cent["gene"].astype(str))
print("[init] target genes:", len(TARGET), "| states:", len(STATES), flush=True)


def load_mtx_stream(f_mat, f_feat, f_bc, key):
    """G x N log1p(CP10K) sparse matrix; cached after first parse."""
    npz = os.path.join(CACHE, key + ".npz")
    bcf = os.path.join(CACHE, key + "_barcodes.txt")
    if os.path.exists(npz) and os.path.exists(bcf):
        X = sparse.load_npz(npz)
        with open(bcf, encoding="utf-8") as h:
            barcodes = [l.strip() for l in h]
        print("    [cache]", key, X.shape, flush=True)
        return X, barcodes
    with gzip.open(f_feat, "rt") as h:
        genes = [l.rstrip("\n").split("\t")[1] for l in h]
    with gzip.open(f_bc, "rt") as h:
        barcodes = [l.strip() for l in h]
    ng, nc = len(genes), len(barcodes)
    rowmap = np.full(ng + 1, -1, dtype=np.int32)
    for i, g in enumerate(genes):
        if g in tpos:
            rowmap[i + 1] = tpos[g]
    total = np.zeros(nc, dtype=np.float64)
    rr, cc, vv = [], [], []
    with gzip.open(f_mat, "rt") as h:
        for line in h:
            if not line.startswith("%"):
                dims = line.split()
                break
        nrows, ncols, nnz = (int(x) for x in dims)
        assert nrows == ng and ncols == nc, (nrows, ncols, ng, nc)
        print("    [parse]", key, "cells", nc, "nnz", nnz, flush=True)
        for k, chunk in enumerate(pd.read_csv(h, sep=r"\s+", header=None,
                                              names=["g", "c", "v"],
                                              dtype=np.int32,
                                              chunksize=3_000_000, engine="c")):
            g = chunk["g"].to_numpy()
            c = chunk["c"].to_numpy() - 1
            v = chunk["v"].to_numpy().astype(np.float64, copy=False)
            total += np.bincount(c, weights=v, minlength=nc)
            p = rowmap[g]
            m = p >= 0
            if m.any():
                rr.append(p[m].astype(np.int32))
                cc.append(c[m].astype(np.int32))
                vv.append(v[m].astype(np.float32))
            if k % 10 == 0:
                print("      chunk", k, flush=True)
    X = sparse.coo_matrix((np.concatenate(vv),
                           (np.concatenate(rr), np.concatenate(cc))),
                          shape=(len(TARGET), nc)).tocsr()
    X.sum_duplicates()
    X = X @ sparse.diags((1e4 / np.maximum(total, 1)).astype(np.float32))
    X.data = np.log1p(X.data)
    sparse.save_npz(npz, X)
    with open(bcf, "w", encoding="utf-8") as h:
        h.write("\n".join(barcodes))
    print("    [done]", key, "target nnz", X.nnz, flush=True)
    return X, barcodes


def load_h5ad_ec():
    import anndata as ad
    npz = os.path.join(CACHE, "cxg_endo.npz")
    metaf = os.path.join(CACHE, "cxg_endo_meta.csv")
    if os.path.exists(npz) and os.path.exists(metaf):
        return sparse.load_npz(npz), pd.read_csv(metaf)
    a = ad.read_h5ad(os.path.join(RAW, "cxg_endo_3951.h5ad"))
    sym = a.var["feature_name"].astype(str).tolist()
    idx = [i for i, g in enumerate(sym) if g in tpos]
    sub = a.X[:, idx]
    sub = sub.T.tocsr() if sparse.issparse(sub) else sparse.csr_matrix(
        np.asarray(sub, dtype=np.float32).T)
    r = np.array([tpos[sym[i]] for i in idx], dtype=np.int32)
    co = sub.tocoo()
    M = sparse.coo_matrix((co.data.astype(np.float32), (r[co.row], co.col)),
                          shape=(len(TARGET), sub.shape[1])).tocsr()
    M.sum_duplicates()
    obs = a.obs[["donor_id", "tissue", "disease", "Patient_status",
                 "Detailed_Cell_Type", "Study"]].copy().reset_index(drop=True)
    sparse.save_npz(npz, M)
    obs.to_csv(metaf, index=False, encoding="utf-8-sig")
    return M, obs


def assign_states(X):
    gc = [g for g in TARGET if g in cent_genes]
    cidx = [tpos[g] for g in gc]
    C = cent.set_index("gene").loc[gc, STATES].to_numpy(dtype=np.float32)
    keep = C.std(axis=1) > 0
    gc = [g for g, k in zip(gc, keep) if k]
    cidx = [i for i, k in zip(cidx, keep) if k]
    C = C[keep]
    Xg = X[cidx, :].toarray().T.astype(np.float32)
    Rq = rankdata(Xg, axis=1).astype(np.float32)
    Rc = rankdata(C, axis=0).astype(np.float32)
    Rq -= Rq.mean(axis=1, keepdims=True)
    Rc -= Rc.mean(axis=0, keepdims=True)
    Rq /= (Rq.std(axis=1, keepdims=True) + 1e-6)
    Rc /= (Rc.std(axis=0, keepdims=True) + 1e-6)
    sim = (Rq @ Rc) / max(len(gc) - 1, 1)
    labA = np.array(STATES)[sim.argmax(axis=1)]
    ms = np.zeros((Xg.shape[0], len(STATES)), dtype=np.float32)
    g2col = {g: j for j, g in enumerate(gc)}
    for j, s in enumerate(STATES):
        cs = [g2col[g] for g in marker_by_state[s] if g in g2col]
        ms[:, j] = Xg[:, cs].mean(axis=1) if cs else 0.0
    ms = (ms - ms.mean(axis=0, keepdims=True)) / (ms.std(axis=0, keepdims=True) + 1e-6)
    labB = np.array(STATES)[ms.argmax(axis=1)]
    return labA, labB, float(np.mean(labA == labB)), len(gc)


def proportions(labels, sample, group, method, cohort):
    d = pd.DataFrame({"label": labels, "sample": np.asarray(sample),
                      "group": np.asarray(group)})
    n = d.groupby(["sample", "group", "label"]).size().unstack(fill_value=0)
    for s in STATES:
        if s not in n.columns:
            n[s] = 0
    n = n[STATES]
    tot = n.sum(axis=1)
    frac = n.div(tot.replace(0, np.nan), axis=0).fillna(0).reset_index()
    frac.insert(0, "cohort", cohort)
    frac.insert(1, "method", method)
    frac["n_cells"] = tot.values
    return frac


def marker_expression(X, sample, group, cohort):
    rows = []
    for s in STATES:
        gi = [tpos[g] for g in marker_by_state[s] if g in tpos]
        sub = X[gi, :]
        pos = (sub > 0).getnnz(axis=0)
        rows.append(dict(cohort=cohort, sample=sample, group=group, state=s,
                         n_markers=len(gi),
                         mean_log_expr=round(float(sub.mean()), 4),
                         frac_cells_ge1_marker=round(float((pos >= 1).mean()), 4),
                         frac_cells_ge2_markers=round(float((pos >= 2).mean()), 4),
                         n_cells=X.shape[1]))
    return rows


props, mex, agree, sens = [], [], [], []
pan_idx = [tpos[g] for g in PAN_EC if g in tpos]
lin_idx = [tpos[g] for g in LINEAGE if g in tpos]


def process(X, sample, group, cohort, filter_ec=True, det_thr=2):
    det = (X[pan_idx, :] > 0).getnnz(axis=0)
    for thr in [1, 2, 3, 4]:
        sens.append(dict(cohort=cohort, sample=str(sample[0])[:40],
                         threshold=thr, n_cells=int((det >= thr).sum()),
                         total_cells=X.shape[1]))
    keep = (det >= det_thr) if filter_ec else np.ones(X.shape[1], dtype=bool)
    print("    EC-like cells:", int(keep.sum()), "/", X.shape[1], flush=True)
    if keep.sum() < 10:
        return
    Xs = X[:, keep]
    s1 = np.asarray(sample)[keep]
    g1 = np.asarray(group)[keep]
    labA, labB, agr, ngene = assign_states(Xs)
    agree.append(dict(cohort=cohort, n_cells=int(keep.sum()),
                      agreement=round(agr, 4), scoring_genes=ngene))
    for meth, lab in [("centroid_rank", labA), ("marker_module", labB)]:
        props.append(proportions(lab, s1, g1, meth, cohort))
    for smp, grp in set(zip(s1, g1)):
        m = (s1 == smp) & (g1 == grp)
        mex.extend(marker_expression(Xs[:, m], smp, grp, cohort))


# ---- 1) CELLxGENE independent normal/benign endothelial atlas ----
print("\n[1] CELLxGENE endothelial atlas", flush=True)
M, obs = load_h5ad_ec()
grp = np.where(obs["disease"].astype(str) == "normal", "Normal_like",
               "Diseased_benign")
smp = np.array(["CXG_" + str(x) for x in obs["donor_id"].astype(str)])
process(M, smp, grp, "CELLxGENE_EC", filter_ec=False)
auth = obs.groupby(["disease", "Detailed_Cell_Type"]).size().reset_index(name="n")
auth["frac"] = auth["n"] / auth.groupby("disease")["n"].transform("sum")
auth.to_csv(os.path.join(OUT, "TableS9_author_annotation_CELLxGENE.csv"),
            index=False, encoding="utf-8-sig")
print("[CELLxGENE] author annotation:",
      obs["Detailed_Cell_Type"].value_counts().to_dict(), flush=True)

# ---- 2) GSE196756 paired tumour / adjacent normal ----
print("\n[2] GSE196756 (3 paired T/N)", flush=True)
for fm in sorted(f for f in os.listdir(G196) if f.endswith("_matrix.mtx.gz")):
    stem = fm.replace("_matrix.mtx.gz", "")
    g = "Tumor" if stem.endswith("T") else "Normal"
    X, bcs = load_mtx_stream(os.path.join(G196, fm),
                             os.path.join(G196, stem + "_features.tsv.gz"),
                             os.path.join(G196, stem + "_barcodes.tsv.gz"),
                             stem)
    process(X, np.full(len(bcs), stem), np.full(len(bcs), g), "GSE196756")

# ---- 3) GSE203067 paired N / T / LN ----
print("\n[3] GSE203067 (4 N / 4 T / 4 LN)", flush=True)
fmat = os.path.join(RAW, "GSE203067_matrix.mtx.gz")
X, bcs = load_mtx_stream(fmat,
                         os.path.join(RAW, "GSE203067_features.tsv.gz"),
                         os.path.join(RAW, "GSE203067_barcodes.tsv.gz"),
                         "GSE203067")
suf = np.array([int(b.rsplit("-", 1)[1]) for b in bcs])
grp = np.array(["Normal" if 5 <= k <= 8 else ("Tumor" if 9 <= k <= 12
               else "LymphNode") for k in suf])
smp = np.array(["N" + str(k - 4) if 5 <= k <= 8 else
                ("T" + str(k - 8) if 9 <= k <= 12 else "L" + str(k)) for k in suf])
process(X, smp, grp, "GSE203067")

# ---- 4) in-house GSE160269 reference (original annotation) ----
ref = pd.read_csv(os.path.join(BASE, "08_scRNA", "endothelial", "revised",
                               "celltype_proportion_closed.csv"))
refp = ref.pivot_table(index=["sample", "type"], columns="celltype",
                       values="frac", aggfunc="mean").reset_index()
refp = refp.rename(columns={"type": "group"})
refp.insert(0, "cohort", "GSE160269_reference")
refp.insert(1, "method", "original_annotation")
cnt = ref.groupby(["sample", "type"])["n"].sum().rename("n_cells").reset_index()
refp = refp.merge(cnt, left_on=["sample", "group"], right_on=["sample", "type"],
                  how="left").drop(columns=["type"])
props.append(refp[["cohort", "method", "sample", "group"] + STATES + ["n_cells"]])

allp = pd.concat(props, ignore_index=True)
allp.to_csv(os.path.join(OUT, "TableS9_state_proportions_by_sample.csv"),
            index=False, encoding="utf-8-sig")
mexdf = pd.DataFrame(mex)
mexdf.to_csv(os.path.join(OUT, "TableS9_marker_expression_by_sample.csv"),
             index=False, encoding="utf-8-sig")
pd.DataFrame(agree).to_csv(os.path.join(OUT, "TableS9_method_agreement.csv"),
                           index=False, encoding="utf-8-sig")
pd.DataFrame(sens).to_csv(os.path.join(OUT, "TableS9_EC_filter_sensitivity.csv"),
                          index=False, encoding="utf-8-sig")

rows = []
for (c, m, g), d in allp.groupby(["cohort", "method", "group"]):
    w = d["n_cells"].clip(lower=1)
    r = dict(cohort=c, method=m, group=g, n_samples=len(d), n_cells=int(d.n_cells.sum()))
    for s in STATES:
        r[s] = float(np.average(d[s], weights=w))
    rows.append(r)
summ = pd.DataFrame(rows)
summ.to_csv(os.path.join(OUT, "TableS9_cohort_summary.csv"),
            index=False, encoding="utf-8-sig")

# marker-based cohort summary (primary readout, no label transfer)
mrows = []
for (c, g, s), d in mexdf.groupby(["cohort", "group", "state"]):
    w = d["n_cells"].clip(lower=1)
    mrows.append(dict(cohort=c, group=g, state=s, n_samples=len(d),
                      n_cells=int(d.n_cells.sum()),
                      mean_log_expr=float(np.average(d.mean_log_expr, weights=w)),
                      frac_ge2=float(np.average(d.frac_cells_ge2_markers, weights=w))))
msumm = pd.DataFrame(mrows)
msumm.to_csv(os.path.join(OUT, "TableS9_marker_cohort_summary.csv"),
             index=False, encoding="utf-8-sig")

print("\n=== state proportions (centroid rank / original) ===", flush=True)
print(summ[summ.method.isin(["centroid_rank", "original_annotation"])]
      [["cohort", "group", "n_samples", "n_cells", "Venous",
        "Tip_angiogenic", "IFNgamma_EC"]].round(3).to_string(index=False), flush=True)
print("\n=== marker positive-cell fraction (primary evidence) ===", flush=True)
print(msumm[msumm.state.isin(["Tip_angiogenic", "IFNgamma_EC", "Venous"])]
      [["cohort", "group", "state", "n_cells", "mean_log_expr", "frac_ge2"]]
      .round(3).to_string(index=False), flush=True)
print("\n=== method agreement ===", flush=True)
print(pd.DataFrame(agree).to_string(index=False), flush=True)

# ---- Figure S4 ----
# 2x2, two-column width (7.1 inch at 600 dpi): matches Figures 1-7 / S1-S3.
# A single 1x4 row would be 19.5 inch wide and exceed the journal's page width.
fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.6))
ax = axes.ravel()
# A: author annotation of independent normal EC atlas
a = auth[auth.disease == "normal"].sort_values("n", ascending=False)
ax[0].bar(range(len(a)), a["frac"] * 100, color="#2e6da4", width=0.65)
ax[0].set_xticks(range(len(a)))
ax[0].set_xticklabels([x.replace("Arterial capillary", "Art. capillary")
                         for x in a["Detailed_Cell_Type"]], fontsize=9)
ax[0].set_ylabel("Endothelial cells (%)")
ax[0].set_title("(A) Independent normal EC atlas\n(author annotation)",
                  fontsize=10.5, fontweight="bold")
ax[0].set_ylim(0, 90)
ax[0].grid(axis="y", alpha=0.3, linewidth=0.6)
refv = float(refp.loc[refp.group == "Normal", "Venous"].median() * 100)
ax[0].axhline(refv, color="#b03a2e", linestyle="--", linewidth=1.3)
ax[0].text(1.45, refv + 2.5, "GSE160269 normal\nVenous %.1f%%" % refv,
             color="#b03a2e", fontsize=7.4, ha="center")

# B: marker expression, tip-like / interferon-responsive by cohort and group
b = msumm[(msumm.state.isin(["Tip_angiogenic", "IFNgamma_EC"])) &
          (msumm.group.isin(["Normal", "Normal_like", "Tumor"]))].copy()
order = ["CELLxGENE_EC", "GSE196756", "GSE203067"]
b["ord"] = b.cohort.map({c: i for i, c in enumerate(order)})
b = b.sort_values(["ord", "group", "state"])
lab = [c.replace("CELLxGENE_EC", "CXG") + " " +
       g.replace("Normal_like", "N").replace("Normal", "N").replace("Tumor", "T") +
       "\n" + s.replace("Tip_angiogenic", "Tip-like").replace("IFNgamma_EC", "Interferon-resp.")
       for c, g, s in zip(b.cohort, b.group, b.state)]
cols = ["#b03a2e" if s == "Tip_angiogenic" else "#e08b2c" for s in b.state]
ax[1].bar(range(len(b)), b["frac_ge2"] * 100, color=cols, width=0.68)
ax[1].set_xticks(range(len(b)))
ax[1].set_xticklabels(lab, fontsize=6.5, rotation=45, ha="right")
ax[1].set_ylabel("EC with >=2 state markers (%)")
ax[1].set_title("(B) State-marker expression\n(no label transfer)",
                  fontsize=10.5, fontweight="bold")
ax[1].grid(axis="y", alpha=0.3, linewidth=0.6)

# C: exploratory label transfer, venous fraction
c = summ[(summ.method.isin(["centroid_rank", "original_annotation"])) &
         (summ.group.isin(["Normal", "Normal_like", "Tumor"]))].copy()
c["ord"] = c.cohort.map({"GSE160269_reference": 0, "CELLxGENE_EC": 1,
                         "GSE196756": 2, "GSE203067": 3})
c = c.sort_values(["ord", "group"])
cl = [co.replace("GSE160269_reference", "GSE160269").replace("CELLxGENE_EC", "CXG") +
      " " + g.replace("Normal_like", "N").replace("Normal", "N").replace("Tumor", "T") +
      "\n(n=" + str(n) + ")" for co, g, n in zip(c.cohort, c.group, c.n_samples)]
ax[2].bar(range(len(c)), c["Venous"] * 100,
            color=["#2e6da4" if g != "Tumor" else "#7f8c8d" for g in c.group], width=0.68)
ax[2].set_xticks(range(len(c)))
ax[2].set_xticklabels(cl, fontsize=6.8, rotation=45, ha="right")
ax[2].set_ylabel("Venous-labelled EC (%)")
ax[2].set_title("(C) Exploratory label transfer\n(rank correlation, low agreement)",
                  fontsize=10.5, fontweight="bold")
ax[2].grid(axis="y", alpha=0.3, linewidth=0.6)

# D: patient-level paired direction (two independent ESCC cohorts, 7 pairs)
pt = pd.read_csv(os.path.join(OUT, "TableS9_paired_sign_test.csv"),
                 encoding="utf-8-sig")
pt = pt.set_index("state").loc[["Venous", "Tip_angiogenic", "IFNgamma_EC"]].reset_index()
dcol = ["#2e6da4", "#b03a2e", "#e08b2c"]
ax[3].bar(range(3), pt.pairs_concordant / pt.n_pairs * 100, color=dcol, width=0.65)
ax[3].set_xticks(range(3))
ax[3].set_xticklabels(["Venous\nN>T", "Tip-like\nN<T",
                       "Interferon-\nresponsive\nN<T"], fontsize=7.6)
ax[3].set_ylabel("Concordant pairs (%)")
ax[3].set_ylim(0, 122)
ax[3].set_title("(D) Paired direction, 7 pairs\n(2 independent ESCC cohorts)",
                  fontsize=10.5, fontweight="bold")
ax[3].grid(axis="y", alpha=0.3, linewidth=0.6)
for i, (k, q) in enumerate(zip(pt.pairs_concordant, pt.fdr_bh)):
    ax[3].text(i, k / 7 * 100 + 2, "%d/7\nq=%.3f" % (k, q),
               ha="center", va="bottom", fontsize=8)
# Frontiers: figure legends are submitted as text, not embedded in the image file.
# The full (A)-(D) legend is provided in the external-validation report (section 6).
plt.tight_layout(h_pad=1.8, w_pad=1.2)
for ext in ["png", "pdf"]:
    plt.savefig(os.path.join(OUT, "FigureS4_external_validation." + ext),
                dpi=600, bbox_inches="tight")
plt.close()
print("\nDONE ->", OUT, flush=True)
