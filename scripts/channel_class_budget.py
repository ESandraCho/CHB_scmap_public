#!/usr/bin/env python3
"""
channel_class_budget.py — within-cell ion-channel CLASS COMPOSITION per group, to
ask whether the ion-channel expression budget (and its surface-vs-SR Ca
dependency) shifts across developmental stage. Read SEPARATELY within each cell
type (AVN; ventricular working myocardium), not as a cell-type contrast.

seed: 42

Six functional classes (config_channel_classes.CLASS_ORDER): surface Ca (antibody-
accessible sarcolemmal Ca entry), SR Ca (intracellular SR Ca store), K, Na, HCN
(funny current), and transporter (Na/K-ATPase). For each group,
each class is reported as a FRACTION of the combined six-class pool (composition that
sums to 1.0 per group). As an internal within-group ratio this cancels a uniform
per-cell / per-dataset scaling and is robust to library-size offset. A
"dependency" ratio = surface_Ca / SR_Ca is also reported (how much the cell leans on antibody-
accessible surface Ca entry vs the buried SR store).

CAVEAT: a fetal-vs-adult composition is cross-dataset (fetal Lim vs adult Kanemaru)
and subject to abundance-dependent capture bias; it is reported direction-only.
The within-dataset numbers (the composition of any single group) are clean.

Outputs:
  results/tables/channel_class_budget.csv
  results/reports/channel_class_budget.md
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
import config_channel_classes as cc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("channel_class_budget")

# group -> (file(s), obs column, label set). fetal working vCM is pooled across the
# two fetal node dissections (cell-weighted), mirroring ca_accessibility_analysis.
# NODE-RESOLVED groups matching F_node_channel_handling (see ca_accessibility_analysis for
# the `col` sentinel convention: "_nodework:<NODE>" / "_bulkvcm").
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
    # Sim2021 whole-CM developmental trajectory (fetal -> child -> adult), single platform/
    # study across the three stages. Whole-CM (no nodal resolution); "_simcm:<stage>".
    ("sim_fetal_CM",    "sim2021_qc.h5ad", "_simcm:fetal", None),
    ("sim_child_CM",    "sim2021_qc.h5ad", "_simcm:child", None),
    ("sim_adult_CM",    "sim2021_qc.h5ad", "_simcm:adult", None),
]


def group_mask(a, col, labset):
    if isinstance(col, str) and col.startswith("_nodework:"):
        m = cfg.kanemaru_node_working_mask(a, col.split(":", 1)[1])   # node-region working CM
    elif col == "_bulkvcm":
        cs = a.obs["cell_state"].astype(str)
        m = cs.str.startswith("vCM").values & cfg.kanemaru_ventricular_mask(a)  # bulk ventricle
    elif isinstance(col, str) and col.startswith("_simcm:"):
        stage = col.split(":", 1)[1]                                  # fetal | child | adult
        return ((a.obs["cat"].astype(str) == "CM")
                & (a.obs["stage_bin"].astype(str) == stage)).values   # Sim has no modality mask
    else:
        m = a.obs[col].astype(str).isin(labset).values
    return m & cfg.kanemaru_platform_mask(a)    # single-platform (Multiome) for Kanemaru


def linear_sum_count(a, genes, mask):
    """Per-gene SUM of linear expression over masked cells + cell count, so
    cell-weighted means can be pooled across files."""
    genes = [g for g in genes if g in a.var_names]
    n = int(mask.sum())
    if not genes or n == 0:
        return {g: 0.0 for g in genes}, n
    X = a[mask, genes].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lin = np.expm1(X)
    return {g: float(lin[:, i].sum()) for i, g in enumerate(genes)}, n


def pooled_linear_means(fnames, col, labset, genes, cache):
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
    return {g: (sums[g] / total if total else 0.0) for g in genes}, total


def main():
    genes = list(cc.ALL_CLASS_GENES)
    cache, rows = {}, []
    for label, fname, col, labset in GROUPS:
        means, n = pooled_linear_means(fname, col, labset, genes, cache)
        if n < cfg.MIN_CELLS_PER_GROUP:
            logger.warning("%s: %d cells, skip", label, n)
            continue
        # class totals = sum of mean linear expression over the class's genes
        ctot = {c: sum(means.get(g, 0.0) for g in cc.CHANNEL_CLASSES[c])
                for c in cc.CLASS_ORDER}
        pool = sum(ctot.values())
        for c in cc.CLASS_ORDER:
            rows.append({"group": label, "cls": c,
                         "expr": ctot[c],
                         "frac": (ctot[c] / pool if pool else np.nan)})
        logger.info("%s: n=%d  pool=%.1f", label, n, pool)
    df = pd.DataFrame(rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "channel_class_budget.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("ion-channel functional-class composition per group: "
                           "summed mean linear expression and within-group fraction "
                           "of the six-class pool "
                           "(surface_Ca, SR_Ca, K, Na, HCN, transporter)") + "\n")
        df.round(4).to_csv(f, index=False)

    # wide fraction table + dependency ratio for the report
    groups = [g[0] for g in GROUPS if (df.group == g[0]).any()]
    wide = df.pivot(index="group", columns="cls", values="frac").reindex(groups)

    def dep(g):
        sub = df[df.group == g].set_index("cls")["expr"]
        s, r = sub.get("surface_Ca", np.nan), sub.get("SR_Ca", np.nan)
        return (s / r) if (r and r > 0) else np.nan

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "channel_class_budget.md"
    with open(rp, "w") as f:
        f.write(cfg.header("ion-channel class composition across dev stage") + "\n\n")
        f.write("# Ion-channel class composition across developmental stage\n\n")
        f.write("Within-cell composition of the six functional channel classes "
                "(surface Ca, SR Ca, K, Na, HCN, transporter), as a fraction of the "
                "combined six-class pool (the table below displays the surface Ca / SR "
                "Ca / K / Na columns; fractions are over all six classes). "
                "As an internal within-group ratio it cancels a uniform per-cell "
                "/ per-dataset scaling. Read SEPARATELY within each cell type: does the "
                "channel budget shift fetal -> adult in the AVN, and in the working "
                "ventricular myocardium (vCM)? The SAN is shown for completeness.\n\n")

        f.write("## Class composition (fraction of six-class pool), per group\n\n")
        f.write("| group | surface Ca | SR Ca | K | Na | surface_Ca/SR_Ca (dependency) |\n")
        f.write("|---|---|---|---|---|---|\n")
        for g in groups:
            r = wide.loc[g]
            f.write(f"| {g} | {r['surface_Ca']:.3f} | {r['SR_Ca']:.3f} | "
                    f"{r['K']:.3f} | {r['Na']:.3f} | {dep(g):.2f} |\n")

        # per-cell-type fetal->adult shift (the actual question)
        f.write("\n## Fetal -> adult shift, within each cell type "
                "(cross-dataset, direction-only)\n\n")
        pairs = [("AVN", "fetal_AVN_pace", "adult_AVN_pace"),
                 ("vCM", "fetal_work_vCM", "adult_work_vCM"),
                 ("SAN", "fetal_SAN_pace", "adult_SAN_pace")]
        f.write("| cell type | class | fetal frac | adult frac | direction |\n")
        f.write("|---|---|---|---|---|\n")
        for ct, fg, ag in pairs:
            if fg not in wide.index or ag not in wide.index:
                continue
            for c in cc.CLASS_ORDER:
                vf, va = wide.loc[fg, c], wide.loc[ag, c]
                d = ("rises adult" if va > vf * 1.1 else
                     "falls adult" if va < vf * 0.9 else "comparable")
                f.write(f"| {ct} | {cc.CLASS_LABELS[c]} | {vf:.3f} | {va:.3f} | {d} |\n")
            df_dep, da_dep = dep(fg), dep(ag)
            d = ("more surface-dependent fetal" if df_dep > da_dep * 1.1 else
                 "more surface-dependent adult" if df_dep < da_dep * 0.9 else "comparable")
            f.write(f"| {ct} | surface_Ca/SR_Ca | {df_dep:.2f} | {da_dep:.2f} | {d} |\n")

        f.write("\n## Interpretation\n\n")
        f.write("- The composition is an internal within-group ratio (robust to "
                "uniform library-size offset). Any single group's composition is "
                "clean; the *fetal-vs-adult* shift within a cell type is cross-dataset "
                "(fetal Lim vs adult Kanemaru) and is therefore direction-only, not a "
                "magnitude claim (abundance-dependent capture bias; see Methods).\n")
        f.write("- The surface_Ca/SR_Ca dependency ratio is the composition view of "
                "the same surface-vs-SR Ca-handling mode reported by NCX:RYR2 in the "
                "accessibility figure.\n\n")
        f.write("**Class gene sets (curated; neuronal/skeletal isoforms excluded):**\n\n")
        for c in cc.CLASS_ORDER:
            f.write(f"- {cc.CLASS_LABELS[c]}: {', '.join(cc.CHANNEL_CLASSES[c])}\n")
        f.write("\n**Caveats:** transcriptional composition, not channel protein / "
                "current density; class fractions depend on the curated gene sets; "
                "cross-stage = cross-dataset (direction-only); n=1 fetal donor/node; "
                "snRNA-seq.\n")

    logger.info("wrote %s + %s", p1.name, rp.name)


if __name__ == "__main__":
    main()
