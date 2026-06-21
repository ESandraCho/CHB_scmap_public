#!/usr/bin/env python3
"""
ca_accessibility_analysis.py — how much of the cardiomyocyte's Ca-handling machinery
is SURFACE (antibody-accessible) vs INTRACELLULAR (SR/ER, inaccessible), and how
the surface T/L subtype split differs between fetal and adult cells.

seed: 42

Rationale: a circulating antibody can only bind cell-surface channels. The
intracellular Ca machinery (RyR2/IP3R/SERCA on the SR/ER) carries the bulk of a
cardiomyocyte's Ca flux but is not reachable from outside the cell; only
surface channels (Cav, NCX, ORAI, TRP) are antibody-accessible. This script
quantifies, per group (fetal/adult conduction + working):
  - surface-accessible fraction of total Ca-machinery expression
  - within the accessible surface Cav, the T-type vs L-type split
on LINEAR expression, using curated subcellular annotation.

Outputs:
  outputs/tables/ca_accessibility.csv
  outputs/reports/ca_accessibility.md
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg
import config_ca_accessibility as acc
import config_channel_subtypes as sub

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("ca_access")

# reuse the conduction groups from the subtype overlay.
# NODE-RESOLVED groups matching F_node_channel_handling: SAN/AVN x pacemaker/working at
# each stage, plus an adult bulk-ventricular reference. The `col` field carries a sentinel
# string interpreted by group_mask:
#   "cell_type"/"cell_state" -> exact-match labset (fetal pace/work, adult pacemaker states)
#   "_nodework:<NODE>"        -> adult working CM dissected from that node region (SAN/AVN)
#   "_bulkvcm"                -> adult bulk ventricular CM (vCM* in LV/RV/SP/AX)
# Fetal working CM are taken per-node from each fetal (dissected SAN / AVN) object.
GROUPS = [
    ("fetal_SAN_pace",  "lim2024_fetal_san_qc.h5ad",  "cell_type",  {"pacemaker_CM"}),
    ("fetal_SAN_work",  "lim2024_fetal_san_qc.h5ad",  "cell_type",  {"working_CM"}),
    ("fetal_AVN_pace",  "lim2024_fetal_avn_qc.h5ad",  "cell_type",  {"pacemaker_CM"}),
    ("fetal_AVN_work",  "lim2024_fetal_avn_qc.h5ad",  "cell_type",  {"working_CM"}),
    ("adult_SAN_pace",  "kanemaru_conduction_qc.h5ad", "cell_state", {"SAN_P_cell"}),
    ("adult_SAN_work",  "kanemaru_conduction_qc.h5ad", "_nodework:SAN", None),
    ("adult_AVN_pace",  "kanemaru_conduction_qc.h5ad", "cell_state", {"AVN_P_cell"}),
    ("adult_AVN_work",  "kanemaru_conduction_qc.h5ad", "_nodework:AVN", None),
    ("adult_vCM_bulk",  "kanemaru_conduction_qc.h5ad", "_bulkvcm",     None),
]


def group_mask(a, col, labset):
    if isinstance(col, str) and col.startswith("_nodework:"):
        node = col.split(":", 1)[1]                          # SAN | AVN
        m = cfg.kanemaru_node_working_mask(a, node)          # node-region working CM
    elif col == "_bulkvcm":
        cs = a.obs["cell_state"].astype(str)
        m = cs.str.startswith("vCM").values & cfg.kanemaru_ventricular_mask(a)  # bulk ventricle
    else:
        m = a.obs[col].astype(str).isin(labset).values
    return m & cfg.kanemaru_platform_mask(a)    # single-platform (Multiome) for Kanemaru


def linear_sum_count(a, genes, mask):
    """Per-gene SUM of linear expression over masked cells, plus the cell count,
    so cell-weighted means can be pooled across files."""
    genes = [g for g in genes if g in a.var_names]
    n = int(mask.sum())
    if not genes or n == 0:
        return {g: 0.0 for g in genes}, n
    X = a[mask, genes].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lin = np.expm1(X)
    return {g: float(lin[:, i].sum()) for i, g in enumerate(genes)}, n


def _percell_surface_frac_one(a, mask, surf_genes, all_genes):
    """Per-cell surface-accessible fraction = surface / total Ca-machinery (linear),
    for masked cells with total > 0. Returns a 1-D array (one value per qualifying cell)."""
    surf_genes = [g for g in surf_genes if g in a.var_names]
    all_genes = [g for g in all_genes if g in a.var_names]
    if mask.sum() == 0:
        return np.array([])
    Xs = a[mask, surf_genes].X
    Xa = a[mask, all_genes].X
    Xs = Xs.toarray() if sparse.issparse(Xs) else np.asarray(Xs)
    Xa = Xa.toarray() if sparse.issparse(Xa) else np.asarray(Xa)
    surf = np.expm1(Xs).sum(1)
    tot = np.expm1(Xa).sum(1)
    ok = tot > 0
    return surf[ok] / tot[ok]


def _write_percell_surface_frac(all_genes, cache):
    surf_genes = list(acc.genes_by_access("SURFACE"))
    rows = []
    for label, fname, col, labset in GROUPS:
        fnames = (fname,) if isinstance(fname, str) else fname
        vals = []
        n_total = 0
        for fn in fnames:
            a = cache.get(fn)
            if a is None:
                a = ad.read_h5ad(cfg.PROC_DIR / fn); cache[fn] = a
            m = group_mask(a, col, labset)
            n_total += int(m.sum())
            vals.append(_percell_surface_frac_one(a, m, surf_genes, all_genes))
        vals = np.concatenate(vals) if vals else np.array([])
        if n_total < cfg.MIN_CELLS_PER_GROUP:
            continue
        for v in vals:
            rows.append({"group": label, "surface_frac": round(float(v), 4)})
    out = pd.DataFrame(rows)
    p = cfg.OUT_TABLES / "ca_accessibility_percell.csv"
    with open(p, "w") as f:
        f.write(cfg.header("PER-CELL surface-accessible Ca fraction (within-sample cell "
                           "dispersion for Fig 2A overlay; NOT donor-level error)") + "\n")
        out.to_csv(f, index=False)
    logger.info("wrote %s (%d cells across %d groups)", p, len(out), out.group.nunique())


def pooled_linear_mean(fnames, col, labset, genes, cache):
    """Cell-weighted mean linear expression pooled across one or more files."""
    if isinstance(fnames, str):
        fnames = (fnames,)
    sums, total = {g: 0.0 for g in genes}, 0
    for fname in fnames:
        a = cache.get(fname)
        if a is None:
            a = ad.read_h5ad(cfg.PROC_DIR / fname)
            cache[fname] = a
        s, n = linear_sum_count(a, genes, group_mask(a, col, labset))
        for g, v in s.items():
            sums[g] += v
        total += n
    means = {g: (sums[g] / total if total else 0.0) for g in genes}
    return means, total


def main():
    genes = list(acc.CA_MACHINERY.keys())
    cache, rows = {}, []
    for label, fname, col, labset in GROUPS:
        means, n = pooled_linear_mean(fname, col, labset, genes, cache)
        if n < cfg.MIN_CELLS_PER_GROUP:
            continue
        for g, m in means.items():
            rows.append({"group": label, "gene": g, "access": acc.access(g),
                         "compartment": acc.compartment(g), "expr": m,
                         "subtype": sub.subtype(g)})
        logger.info("%s: n=%d", label, n)
    df = pd.DataFrame(rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "ca_accessibility.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("Ca-machinery expression by subcellular accessibility "
                           "(surface=antibody-accessible vs intracellular SR/ER) per "
                           "conduction group; linear scale") + "\n")
        df.round(4).to_csv(f, index=False)

    # PER-CELL surface-accessible fraction (for the Fig 2A distribution overlay): the
    # group bar is a pooled per-cell-mean ratio; this records the cell-to-cell spread of
    # surface/(surface+intracellular+mixed) Ca-machinery within each cell, for cells whose
    # total Ca-machinery > 0. This is WITHIN-SAMPLE cell dispersion, NOT donor-level error.
    _write_percell_surface_frac(genes, cache)

    # per-group: accessible fraction of total Ca machinery, + surface T/L split
    groups = [g[0] for g in GROUPS if g[0] in df.group.unique()]
    summ = []
    for g in groups:
        d = df[df.group == g]
        total = d["expr"].sum()
        surf = d[d.access == "SURFACE"]["expr"].sum()
        intra = d[d.access == "INTRACELLULAR"]["expr"].sum()
        mixed = d[d.access == "MIXED"]["expr"].sum()
        # within surface Cav, T vs L
        cav = d[d.gene.isin(sub.genes_for_ion("Ca"))]
        L = cav[cav.subtype == "L-type"]["expr"].sum()
        T = cav[cav.subtype == "T-type"]["expr"].sum()
        summ.append({"group": g,
                     "surface_frac": surf / total if total else np.nan,
                     "intracellular_frac": intra / total if total else np.nan,
                     "mixed_frac": mixed / total if total else np.nan,
                     "surface_Cav_T_over_TplusL": T / (T + L) if (T + L) else np.nan,
                     "RYR2_expr": float(d[d.gene == "RYR2"]["expr"].iloc[0]) if (d.gene == "RYR2").any() else np.nan,
                     "total_Ca_machinery": total})
    sdf = pd.DataFrame(summ)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "ca_accessibility.md"
    with open(rp, "w") as f:
        f.write(cfg.header("Ca-machinery autoantibody accessibility") + "\n\n")
        f.write("# What fraction of the Ca machinery can an autoantibody even reach?\n\n")
        f.write("A circulating antibody binds only SURFACE channels. Intracellular Ca "
                "machinery (RyR2/IP3R/SERCA on SR/ER) carries the bulk of a "
                "cardiomyocyte's Ca flux but is INACCESSIBLE to antibody. The "
                "antibody-accessible Ca target set is therefore the surface slice. "
                "Linear expression; conduction groups.\n\n")
        f.write("## Surface vs intracellular Ca machinery, per group\n\n")
        f.write("| group | surface % | intracellular % | mixed % | surface-Cav T/(T+L) | RYR2 (intra, ref) |\n")
        f.write("|---|---|---|---|---|---|\n")
        for _, r in sdf.iterrows():
            f.write(f"| {r['group']} | {100*r['surface_frac']:.1f}% | "
                    f"{100*r['intracellular_frac']:.1f}% | {100*r['mixed_frac']:.1f}% | "
                    f"{r['surface_Cav_T_over_TplusL']:.3f} | {r['RYR2_expr']:.2f} |\n")
        f.write("\n## Interpretation\n\n")
        f.write("- The dominant Ca flux (RyR2, SR) is INTRACELLULAR and not "
                "antibody-accessible. The surface-accessible Ca channels are a "
                "minority of total Ca-machinery expression. Sarcolemmal Ca entry "
                "triggers SR release (Ca-induced Ca release), so surface channels "
                "are the entry point on which SR release depends. In the immature "
                "fetal node the SR is sparse and the cell relies more on sarcolemmal "
                "Ca entry and NCX (Seki 2003; Artman 1997).\n")
        f.write("- Within the surface-accessible Cav, the fetal node shows a higher "
                "T-type fraction than the adult node (see surface-Cav T/(T+L) column; "
                "cf. channel_subtype_conduction). The surface Ca subtype that differs "
                "most with development is the T-type (Cav3.1/CACNA1G).\n\n")
        f.write("**Caveats:** mRNA abundance is NOT protein localization — accessibility "
                "is assigned from CURATED UniProt subcellular annotation + channel "
                "biology, not measured here; transcript level is a proxy for channel "
                "availability. Cross-dataset, n=1 fetal/node, snRNA-seq.\n")

    logger.info("wrote %s + %s | surface frac range %.2f-%.2f",
                p1.name, rp.name, sdf.surface_frac.min(), sdf.surface_frac.max())


if __name__ == "__main__":
    main()
