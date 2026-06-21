#!/usr/bin/env python3
"""
nodal_specificity_interaction.py — formal "constitutive vs fetal-specific" test for the
conduction-cell channel panel, via a cell-type x stage INTERACTION.

seed: 42

The problem this solves: fetal (GW19) and adult nodal data are different datasets with no
developmentally-matched shared cell type, so a direct fetal-vs-adult pacemaker comparison is
batch+stage confounded. But each dataset's pacemaker-vs-working contrast is internal and
batch-free. The well-posed fetal-vs-adult question is therefore the DIFFERENCE of the two
within-dataset log-ratios (a difference-in-differences / interaction):

    interaction = (fetal pacemaker - fetal working) - (adult pacemaker - adult working)

For each gene and node (SAN, AVN) the four cell groups are pooled and fit, on log1p expression,

    expr ~ C(celltype) * C(stage)        (celltype: pacemaker/working ; stage: fetal/adult)

The INTERACTION coefficient is exactly the difference-of-within-dataset-ratios above, and its
p-value tests whether a channel's pacemaker-SPECIFICITY differs between fetal and adult:
  interaction ~ 0  -> pacemaker-specificity is the SAME at both stages  -> CONSTITUTIVE nodal
  interaction > 0  -> MORE pacemaker-specific in fetal                  -> FETAL-specific
  interaction < 0  -> more pacemaker-specific in adult.
Crucially the STAGE MAIN EFFECT (confounded by batch/platform) is NOT interpreted; only the
interaction, which is batch-robust because each ratio is internal to one dataset.

HONEST CAVEAT (encoded in the report): the fetal side is n=1 donor per node, so cells are
pseudoreplicated on the fetal arm — the interaction p-values treat cells as independent and
therefore OVERSTATE significance. They are reported as a descriptive effect-size + nominal p,
NOT as donor-level inference. Direction/magnitude is the message.

Inputs (processed objects): fetal Lim/Protze (pacemaker_CM/working_CM), adult Kanemaru
(SAN_P_cell/AVN_P_cell vs Ventricular Cardiomyocyte). Adult X is normalized/log (no raw
counts) -> compared on the same log1p footing as fetal.

Outputs:
  results/tables/nodal_specificity_interaction.csv
  results/reports/nodal_specificity_interaction.md
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("nodal_interaction")

PANEL = ["CACNA1C", "CACNA1D", "CACNA1G", "CACNA1H", "CACNA1E", "SCN5A", "KCNH2", "HCN4", "GJA5"]

# per node: (fetal file, fetal pace label, fetal work label), adult uses Kanemaru
FETAL = {"SAN": "lim2024_fetal_san_qc.h5ad", "AVN": "lim2024_fetal_avn_qc.h5ad"}
ADULT = "kanemaru_conduction_qc.h5ad"
ADULT_PACE = {"SAN": "SAN_P_cell", "AVN": "AVN_P_cell"}
ADULT_WORK_CTYPE = "Ventricular Cardiomyocyte"


def col_log1p(a, g):
    """log1p-normalized expression vector for gene g (fetal X and Kanemaru X are both
    log1p-normalized, so they are on the same footing)."""
    x = a[:, g].X
    return (x.toarray() if sparse.issparse(x) else np.asarray(x)).ravel()


def main():
    fetal_obj = {n: ad.read_h5ad(cfg.PROC_DIR / f) for n, f in FETAL.items()}
    adult = ad.read_h5ad(cfg.PROC_DIR / ADULT)
    a_cs = adult.obs["cell_state"].astype(str).values
    # single-platform (Multiome) for Kanemaru: nodal states are Multiome-only, so restrict
    # to Multiome to avoid a platform main effect. The adult working reference is
    # NODE-MATCHED (working CMs from the SAME node region as the pacemaker), mirroring the
    # fetal node-dissection working_CM — so the interaction is not confounded by comparing
    # fetal node-adjacent working CM against adult bulk ventricle.
    a_plat = cfg.kanemaru_platform_mask(adult)

    rows = []
    for node in ("SAN", "AVN"):
        fa = fetal_obj[node]
        f_ct = fa.obs["cell_type"].astype(str).values
        f_pace = f_ct == "pacemaker_CM"
        f_work = f_ct == "working_CM"
        a_pace = (a_cs == ADULT_PACE[node]) & a_plat
        a_work = cfg.kanemaru_node_working_mask(adult, node) & a_plat   # node-region working
        if f_pace.sum() < 2 or a_pace.sum() < 2:
            logger.warning("%s: too few pacemaker cells", node); continue
        genes = [g for g in PANEL if g in fa.var_names and g in adult.var_names]
        for g in genes:
            # assemble the 4-group long frame on log1p expression
            fp, fw = col_log1p(fa, g)[f_pace], col_log1p(fa, g)[f_work]
            ap, aw = col_log1p(adult, g)[a_pace], col_log1p(adult, g)[a_work]
            df = pd.DataFrame({
                "expr": np.concatenate([fp, fw, ap, aw]),
                "celltype": (["pace"] * len(fp) + ["work"] * len(fw)
                             + ["pace"] * len(ap) + ["work"] * len(aw)),
                "stage": (["fetal"] * (len(fp) + len(fw)) + ["adult"] * (len(ap) + len(aw))),
            })
            # baseline = work, adult; interaction = pace:fetal
            m = smf.ols("expr ~ C(celltype, Treatment('work')) * C(stage, Treatment('adult'))",
                        data=df).fit()
            iname = [t for t in m.params.index if ":" in t]
            inter = float(m.params[iname[0]]) if iname else np.nan
            ip = float(m.pvalues[iname[0]]) if iname else np.nan
            fetal_ratio = float(fp.mean() - fw.mean())     # log1p-mean within-fetal pace-work
            adult_ratio = float(ap.mean() - aw.mean())     # within-adult
            if not np.isfinite(inter):
                verdict = "n/a"
            elif ip >= 0.05:
                verdict = "constitutive (both stages)"
            elif inter > 0:
                verdict = "fetal-specific"
            else:
                verdict = "adult-specific"
            rows.append({"node": node, "gene": g,
                         "fetal_pace_minus_work": round(fetal_ratio, 4),
                         "adult_pace_minus_work": round(adult_ratio, 4),
                         "interaction_coef": round(inter, 4), "interaction_p": ip,
                         "verdict": verdict,
                         "n_fetal_pace": int(f_pace.sum()), "n_adult_pace": int(a_pace.sum())})
        logger.info("%s done", node)
    out = pd.DataFrame(rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "nodal_specificity_interaction.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("celltype x stage interaction per channel: (fetal pace-work) - "
                           "(adult pace-work) on log1p expression; interaction = batch-robust "
                           "fetal-vs-adult pacemaker-specificity difference. STAGE MAIN EFFECT "
                           "NOT interpreted (batch-confounded). n=1 fetal donor -> cells "
                           "pseudoreplicated, p nominal/descriptive only") + "\n")
        out.round(4).to_csv(f, index=False)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "nodal_specificity_interaction.md"
    with open(rp, "w") as f:
        f.write(cfg.header("nodal channel specificity: constitutive vs fetal-specific") + "\n\n")
        f.write("# Is each conduction-cell channel constitutively nodal, or fetal-specific?\n\n")
        f.write("Fetal and adult nodal datasets share no developmentally-matched cell type, so "
                "a direct fetal-vs-adult pacemaker comparison is batch+stage confounded. The "
                "well-posed test is the DIFFERENCE of the two within-dataset pacemaker/working "
                "log-ratios (a celltype x stage interaction). The interaction is batch-robust "
                "(each ratio internal to one dataset); the stage MAIN effect is confounded and "
                "is not interpreted.\n\n")
        f.write("| node | gene | fetal pace-work | adult pace-work | interaction | p | verdict |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, r in out.iterrows():
            f.write(f"| {r['node']} | {r['gene']} | {r['fetal_pace_minus_work']:+.2f} | "
                    f"{r['adult_pace_minus_work']:+.2f} | {r['interaction_coef']:+.2f} | "
                    f"{r['interaction_p']:.1e} | {r['verdict']} |\n")
        f.write("\n## Reading\n\n")
        f.write("- **constitutive (both stages):** pacemaker-enriched (or depleted) to the same "
                "degree in fetal and adult — a permanent property of the conduction cell. The "
                "anti-Ro Ca targets that are constitutive mean the antibody always has a nodal "
                "target; fetal vulnerability then comes from the developmental LEVEL axis "
                "(higher fetal CM expression, Fig 3) + accessibility (Fig 2), not from the "
                "channel being nodal only in the fetus.\n")
        f.write("- **fetal-specific:** more pacemaker-specific in fetal than adult — a candidate "
                "for a fetal-only conduction vulnerability.\n\n")
        f.write("**Caveat (important):** the fetal side is n=1 donor per node, so cells are "
                "pseudoreplicated and the interaction p-values OVERSTATE significance (they treat "
                "cells as independent). They are descriptive effect-sizes with nominal p, not "
                "donor-level inference; the direction and magnitude carry the message. Adult X is "
                "log-normalized (no raw counts); fetal and adult are compared on the same log1p "
                "footing.\n")

    logger.info("wrote %s + %s", p1.name, rp.name)


if __name__ == "__main__":
    main()
