#!/usr/bin/env python3
"""
adult_pacemaker_vs_working.py — the ADULT mirror of fetal_pacemaker_vs_working.py:
within-adult-dataset contrast of conduction (pacemaker) cells vs working myocardium for
the same channel panel, in the Kanemaru adult conduction atlas.

seed: 42

WHY this exists (the methodological point): fetal (Lim/Protze, GW19) and adult (Kanemaru)
nodal data are different donors / platforms / labs, AND there is NO developmentally-matched
shared cell type between them (the fetal "working CM" is GW19, the adult "working CM" is
20-75 y), so a DIRECT fetal-pacemaker vs adult-pacemaker expression comparison is confounded
by both batch and developmental stage and cannot be made cleanly. What CAN be compared is the
WITHIN-DATASET pacemaker/working contrast in each: each ratio is internal to one dataset and
therefore batch-free. Comparing the fetal log2FC to the adult log2FC then asks a well-posed
question — is a channel's pacemaker-SPECIFICITY a fetal feature, an adult feature, or shared?
— without ever subtracting two cross-dataset absolute levels. This script produces the adult
half so the two within-dataset contrasts can be placed side by side (Fig 1B,C vs the adult).

NOTE: the Kanemaru processed object carries only normalized/log expression (no raw counts);
the contrast is computed on that normalized X (expm1 -> linear mean), consistently for both
pacemaker and working cells of the SAME object, so the within-dataset ratio is valid. This is
a cell-type-SPECIFICITY contrast, not an absolute-level claim.

Outputs:
  results/tables/adult_pacemaker_vs_working.csv
  results/reports/adult_pacemaker_vs_working.md
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
logger = logging.getLogger("adult_pace_work")

ADULT = "kanemaru_conduction_qc.h5ad"
# adult node pacemaker cell_state -> labelled like the fetal nodes (SAN, AVN)
NODES = [("SAN", "SAN_P_cell"), ("AVN", "AVN_P_cell")]
WORKING_CELLTYPE = "Ventricular Cardiomyocyte"   # adult working myocardium reference

# same SHARED Fig-1 panel as the fetal script, for a like-for-like comparison
_ROLE = {"CACNA1C": "Ca target", "CACNA1D": "Ca target", "CACNA1G": "Ca target",
         "CACNA1H": "Ca target", "KCND3": "K (non-targeted control)",
         "SCN5A": "Na (node-poor test)", "KCNH2": "K target (LQT)", "KCNQ1": "K target",
         "KCNA4": "K target", "KCNJ5": "K target (AF)", "KCNJ11": "K target (VF)",
         "HCN4": "pacemaker identity", "GJA5": "fast-conduction (Cx40)"}
PANEL = {g: _ROLE.get(g, "channel") for g in cfg.FIG1_PANEL_GENES + cfg.FIG1_IDENTITY_GENES}


def linvec(a, g):
    """linear expression of gene g per cell: expm1 of the normalized/log X."""
    x = a[:, g].X
    x = (x.toarray() if sparse.issparse(x) else np.asarray(x)).ravel()
    return np.expm1(x)


def main():
    a = ad.read_h5ad(cfg.PROC_DIR / ADULT)
    cs = a.obs["cell_state"].astype(str).values
    # single-platform (Multiome) for Kanemaru: the SAN_P/AVN_P pacemaker states are
    # Multiome-only, so restrict to Multiome (no platform mix). The working reference is
    # NODE-MATCHED — working CMs dissected from the SAME node region as each pacemaker
    # (mirrors the fetal node-dissection working_CM), not bulk ventricle.
    plat = cfg.kanemaru_platform_mask(a)

    rows = []
    for node, state in NODES:
        pace = (cs == state) & plat
        work = cfg.kanemaru_node_working_mask(a, node) & plat   # working CM from this node region
        logger.info("%s: %s pacemaker=%d, node-working=%d", node, state, int(pace.sum()), int(work.sum()))
        if pace.sum() < 2:
            logger.warning("%s: too few pacemaker cells, skipping", node); continue
        genes = [g for g in PANEL if g in a.var_names]
        pvals, recs = [], []
        for g in genes:
            v = linvec(a, g)
            vp, vw = v[pace], v[work]
            mp, mw = float(vp.mean()), float(vw.mean())
            log2fc = float(np.log2((mp + cfg.EPS) / (mw + cfg.EPS)))
            p = float(mannwhitneyu(vp, vw, alternative="two-sided").pvalue)
            pvals.append(p)
            recs.append({"node": node, "gene": g, "role": PANEL[g],
                         "mean_pacemaker": mp, "mean_working": mw,
                         "log2FC_pace_vs_work": log2fc,
                         "n_pacemaker": int(pace.sum()), "n_working": int(work.sum()),
                         "p": p})
        ok = np.array([np.isfinite(p) for p in pvals])
        fdr = np.full(len(pvals), np.nan)
        if ok.any():
            fdr[ok] = false_discovery_control(np.array(pvals)[ok])
        for r, q in zip(recs, fdr):
            r["fdr"] = float(q) if np.isfinite(q) else np.nan
            rows.append(r)
    df = pd.DataFrame(rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "adult_pacemaker_vs_working.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("within-ADULT-dataset pacemaker (SAN_P/AVN_P) vs working (vCM) "
                           "channel-panel contrast (Kanemaru; no cross-dataset confound); "
                           "expm1 of normalized X, Mann-Whitney + BH-FDR. Cell-type-specificity "
                           "contrast — pairs with fetal_pacemaker_vs_working for fetal-vs-adult "
                           "ratio-of-ratios") + "\n")
        df.round(4).to_csv(f, index=False)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "adult_pacemaker_vs_working.md"

    def lf(node, gene):
        r = df[(df.node == node) & (df.gene == gene)]
        return float(r["log2FC_pace_vs_work"].iloc[0]) if len(r) else np.nan

    with open(rp, "w") as f:
        f.write(cfg.header("adult pacemaker vs working CM — channel panel (Kanemaru)") + "\n\n")
        f.write("# The adult conduction cell vs adult working myocardium (within-dataset)\n\n")
        f.write("Mirror of the fetal pacemaker-vs-working contrast, computed within the adult "
                "Kanemaru atlas (SAN_P / AVN_P pacemaker cells vs ventricular working CM). "
                "Because fetal (GW19) and adult nodal data share no developmentally-matched "
                "cell type, a direct fetal-vs-adult pacemaker comparison is confounded; the "
                "valid comparison is each dataset's INTERNAL pacemaker/working ratio, placed "
                "side by side. log2FC = pacemaker vs working; `*` = BH-FDR<0.05.\n\n")
        f.write("| gene | role | SAN log2FC | AVN log2FC |\n|---|---|---|---|\n")
        for g, role in PANEL.items():
            def cell(node):
                r = df[(df.node == node) & (df.gene == g)]
                if not len(r):
                    return "—"
                r = r.iloc[0]
                s = "*" if np.isfinite(r["fdr"]) and r["fdr"] < 0.05 else ""
                return f"{r['log2FC_pace_vs_work']:+.2f}{s}"
            f.write(f"| {g} | {role} | {cell('SAN')} | {cell('AVN')} |\n")
        f.write("\n**Interpretation (ratio-of-ratios with the fetal contrast):** a channel "
                "pacemaker-enriched in BOTH fetal and adult is a constitutive pacemaker channel; "
                "one enriched in fetal but NOT adult pacemaker is a fetal-specific conduction "
                "feature (the candidate vulnerability). This avoids ever comparing fetal and "
                "adult absolute levels directly.\n\n")
        f.write("**Caveats:** Kanemaru X is normalized/log (no raw counts) — contrast on expm1, "
                "valid as a within-dataset specificity ratio, not an absolute level; adult AVN "
                "P-cells n small; multi-donor adult vs single-donor fetal (the fetal side is the "
                "n=1 limit, not this side).\n")

    logger.info("wrote %s + %s", p1.name, rp.name)


if __name__ == "__main__":
    main()
