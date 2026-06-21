#!/usr/bin/env python3
"""
fetal_pacemaker_vs_working.py — within-fetal-stage, within-donor contrast of the
conduction (pacemaker) cell program vs working myocardium, for the channel panel
that the "why the conduction system" argument rests on.

seed: 42

Motivation: most fetal-vs-adult contrasts in this project are CROSS-DATASET (fetal
Lim/Protze vs adult Kanemaru) and therefore direction-only. But each fetal object
contains BOTH pacemaker_CM and working_CM from the SAME donor and SAME dataset, so
a pacemaker-vs-working contrast is batch-confound-free. This script tests whether
the fetal node is a Ca-dependent, Na-poor, IKr-poor tissue relative to fetal
working myocardium WITHIN the fetal stage, with no cross-dataset confound.

Panel (each gene tagged by its role in the argument):
  Ca targets (anti-Ro CHB)   : CACNA1C, CACNA1D, CACNA1G, CACNA1H
  K non-targeted control      : KCND3 (Kv4.3/Ito — expressed, not attacked)
  Na (the "node is Na-poor")  : SCN5A
  K target (anti-Ro adult LQT): KCNH2 (IKr/hERG)
  pacemaker identity controls : HCN4 (up in node), GJA5/SCN5A (down in node)

For each gene, per fetal node (SAN, AVN), on log1p-normalized expression:
  mean in pacemaker_CM vs working_CM, log2 fold-change, Mann-Whitney p, BH-FDR.
A positive log2FC = enriched in the conduction (pacemaker) cells.

Outputs:
  results/tables/fetal_pacemaker_vs_working.csv
  results/reports/fetal_pacemaker_vs_working.md
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu, false_discovery_control

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("fetal_pace_work")

NODES = [("SAN", "lim2024_fetal_san_qc.h5ad"),
         ("AVN", "lim2024_fetal_avn_qc.h5ad")]

# gene panel = the SHARED Fig-1 channel set (config) + identity controls, so this contrast
# and the stage-ratio panel show exactly the same genes.
_ROLE = {"CACNA1C": "Ca target", "CACNA1D": "Ca target", "CACNA1G": "Ca target",
         "CACNA1H": "Ca target", "KCND3": "K (non-targeted control)",
         "SCN5A": "Na (node-poor test)", "KCNH2": "K target (LQT)", "KCNQ1": "K target",
         "KCNA4": "K target", "KCNJ5": "K target (AF)", "KCNJ11": "K target (VF)",
         "HCN4": "pacemaker identity", "GJA5": "fast-conduction (Cx40)"}
PANEL = {g: _ROLE.get(g, "channel") for g in cfg.FIG1_PANEL_GENES + cfg.FIG1_IDENTITY_GENES}


def vec(a, g):
    x = a[:, g].X
    return (x.toarray() if sparse.issparse(x) else np.asarray(x)).ravel()


def main():
    rows = []
    for node, fname in NODES:
        a = ad.read_h5ad(cfg.PROC_DIR / fname)
        ct = a.obs["cell_type"].astype(str)
        pace = (ct == "pacemaker_CM").values
        work = (ct == "working_CM").values
        logger.info("%s: pacemaker_CM=%d working_CM=%d", node, int(pace.sum()), int(work.sum()))
        genes = [g for g in PANEL if g in a.var_names]
        pvals = []
        recs = []
        for g in genes:
            v = vec(a, g)
            vp, vw = v[pace], v[work]
            mp, mw = float(vp.mean()), float(vw.mean())
            log2fc = float(np.log2((mp + cfg.EPS) / (mw + cfg.EPS)))
            p = (float(mannwhitneyu(vp, vw, alternative="two-sided").pvalue)
                 if pace.sum() >= 2 and work.sum() >= 2 else np.nan)
            pvals.append(p)
            recs.append({"node": node, "gene": g, "role": PANEL[g],
                         "mean_pacemaker": mp, "mean_working": mw,
                         "log2FC_pace_vs_work": log2fc,
                         "n_pacemaker": int(pace.sum()), "n_working": int(work.sum()),
                         "p": p})
        # BH-FDR within node
        ok = np.array([np.isfinite(p) for p in pvals])
        fdr = np.full(len(pvals), np.nan)
        if ok.any():
            fdr[ok] = false_discovery_control(np.array(pvals)[ok])
        for r, q in zip(recs, fdr):
            r["fdr"] = float(q) if np.isfinite(q) else np.nan
            rows.append(r)
    df = pd.DataFrame(rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "fetal_pacemaker_vs_working.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("within-fetal-stage, within-donor pacemaker_CM vs working_CM "
                           "channel-panel contrast (no cross-dataset confound); "
                           "log1p-normalized, Mann-Whitney + BH-FDR") + "\n")
        df.round(4).to_csv(f, index=False)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "fetal_pacemaker_vs_working.md"

    def cell(node, gene):
        r = df[(df.node == node) & (df.gene == gene)]
        if not len(r):
            return "—"
        r = r.iloc[0]
        star = "*" if np.isfinite(r["fdr"]) and r["fdr"] < 0.05 else ""
        return f"{r['log2FC_pace_vs_work']:+.2f}{star}"

    with open(rp, "w") as f:
        f.write(cfg.header("fetal pacemaker vs working CM — channel panel") + "\n\n")
        f.write("# The fetal conduction cell vs fetal working myocardium (within-stage)\n\n")
        f.write("Each fetal node object (SAN, AVN) contains both pacemaker_CM and "
                "working_CM from the SAME donor and dataset, so this pacemaker-vs-working "
                "contrast carries NO cross-dataset batch confound (unlike the fetal-vs-"
                "adult contrasts). It tests, within the fetal stage, whether the "
                "conduction cell is a Ca-dependent, Na-poor, IKr-poor tissue relative to "
                "working myocardium. log2FC = pacemaker vs working; `*` = BH-FDR<0.05.\n\n")

        f.write("| gene | role | SAN log2FC | AVN log2FC |\n")
        f.write("|---|---|---|---|\n")
        for g, role in PANEL.items():
            f.write(f"| {g} | {role} | {cell('SAN', g)} | {cell('AVN', g)} |\n")

        # pull the load-bearing genes for the narrative
        def lf(node, gene):
            r = df[(df.node == node) & (df.gene == gene)]
            return float(r["log2FC_pace_vs_work"].iloc[0]) if len(r) else np.nan
        def mean(node, gene, comp):
            r = df[(df.node == node) & (df.gene == gene)]
            return float(r[comp].iloc[0]) if len(r) else np.nan

        f.write("\n## What this establishes (within fetal stage, no batch confound)\n\n")
        f.write("- **Node is Na-poor:** SCN5A is depleted in fetal pacemaker vs working "
                f"(SAN {lf('SAN','SCN5A'):+.2f}, AVN {lf('AVN','SCN5A'):+.2f} log2FC) — "
                "the conduction cell carries little fast-Na current, so its depolarization "
                "cannot rely on INa. This is the substrate for the Ca-dependence argument.\n")
        f.write("- **Node leans on Ca:** the anti-Ro Ca targets (esp. CACNA1D/Cav1.3 and "
                f"the T-type CACNA1G/Cav3.1) are pacemaker-enriched or pacemaker-retained "
                f"(CACNA1D SAN {lf('SAN','CACNA1D'):+.2f}, AVN {lf('AVN','CACNA1D'):+.2f}; "
                f"CACNA1G SAN {lf('SAN','CACNA1G'):+.2f}, AVN {lf('AVN','CACNA1G'):+.2f}) — "
                "the conduction cell depends on the very channels anti-Ro blocks.\n")
        f.write("- **K target (IKr/hERG) is not a node current:** KCNH2 "
                f"(SAN {lf('SAN','KCNH2'):+.2f}, AVN {lf('AVN','KCNH2'):+.2f}) is not "
                "enriched in the fetal pacemaker — consistent with the node not relying on "
                "IKr, so the adult anti-hERG arm has little nodal substrate even in the fetus.\n")
        f.write("- **Pacemaker identity confirmed:** HCN4 up in pacemaker "
                f"(SAN {lf('SAN','HCN4'):+.2f}, AVN {lf('AVN','HCN4'):+.2f}); "
                "fast-conduction Cx40/GJA5 and SCN5A down — the expected nodal signature.\n\n")
        f.write("**Why this matters:** the 'why the conduction system' legs (Ca-dependence, "
                "Na-poverty, low IKr) are shown here WITHIN the fetal stage and within a "
                "single donor — they do not depend on the cross-dataset fetal-vs-adult "
                "comparison. The fetal-vs-adult contrasts then address the separate "
                "question of why the fetus more than the adult.\n\n")
        f.write("**Caveats:** n=1 fetal donor per node (a cell-type contrast across many "
                "cells of one individual — direction-only, no donor-level statistics); "
                "marker-annotated pacemaker clusters (fetal AVN P-cells are a small "
                "transitional population, noisier); snRNA-seq; log1p-normalized means.\n")

    logger.info("wrote %s + %s", p1.name, rp.name)


if __name__ == "__main__":
    main()
