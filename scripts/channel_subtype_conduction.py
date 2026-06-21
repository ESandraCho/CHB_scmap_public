#!/usr/bin/env python3
"""
channel_subtype_conduction.py — conduction-system overlay of the channel-subtype-
ratio method, restricted to pacemaker cells.

seed: 42

Working myocardium switches the T-type Ca component off early in development, while
the cardiac conduction system retains a T-type component. This script applies the
subtype-fraction method (linear scale, fraction-of-ion-family) to PACEMAKER cells
specifically, node-matched fetal-vs-adult:
  - fetal SAN pacemaker (Lim 2024) vs adult SAN_P (Kanemaru)
  - fetal AVN pacemaker (Protze GSE297072, dissected tissue) vs adult AVN_P (Kanemaru)
and pacemaker-vs-working within each, to compare the T-type-leaning of the fetal
node Ca program with the adult node (and the K/Na/HCN subtype programs).

Outputs:
  outputs/tables/channel_subtype_conduction.csv
  outputs/reports/channel_subtype_conduction.md
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
import config_channel_subtypes as sub

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("subtype_cond")

# (group_label, dataset_file, obs_col, set-of-labels)
GROUPS = [
    ("fetal_SAN_pace",  "lim2024_fetal_san_qc.h5ad",  "cell_type",  {"pacemaker_CM"}),
    ("fetal_SAN_work",  "lim2024_fetal_san_qc.h5ad",  "cell_type",  {"working_CM"}),
    ("fetal_AVN_pace",  "lim2024_fetal_avn_qc.h5ad",  "cell_type",  {"pacemaker_CM"}),
    ("fetal_AVN_work",  "lim2024_fetal_avn_qc.h5ad",  "cell_type",  {"working_CM"}),
    ("adult_SAN_pace",  "kanemaru_conduction_qc.h5ad", "cell_state", {"SAN_P_cell"}),
    ("adult_AVN_pace",  "kanemaru_conduction_qc.h5ad", "cell_state", {"AVN_P_cell"}),
    ("adult_work_vCM",  "kanemaru_conduction_qc.h5ad", "cell_state", None),  # vCM* prefix
]


def linear_mean(a, genes, mask):
    genes = [g for g in genes if g in a.var_names]
    if not genes or mask.sum() == 0:
        return {g: 0.0 for g in genes}
    X = a[mask, genes].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lin = np.expm1(X)
    return {g: float(lin[:, i].mean()) for i, g in enumerate(genes)}


def subtype_fracs(means):
    rows = [{"gene": g, "ion": sub.ion(g), "subtype": sub.subtype(g), "expr": m}
            for g, m in means.items()]
    df = pd.DataFrame(rows)
    df = df[df["ion"].notna()]
    st = df.groupby(["ion", "subtype"], as_index=False)["expr"].sum()
    st["ion_total"] = st.groupby("ion")["expr"].transform("sum")
    st["fraction"] = np.where(st["ion_total"] > 0, st["expr"] / st["ion_total"], np.nan)
    return st


def main():
    genes = list(sub.SUBTYPES.keys())
    cache = {}
    rows = []
    for label, fname, col, labset in GROUPS:
        a = cache.get(fname) or ad.read_h5ad(cfg.PROC_DIR / fname)
        cache[fname] = a
        c = a.obs[col].astype(str)
        if labset is None:                            # node-matched working reference
            # working CMs from the node regions (SAN+AVN), mirroring fetal node-dissection
            # working_CM — NOT bulk ventricle.
            mask = cfg.kanemaru_node_working_mask(a)
        else:
            mask = c.isin(labset).values
        mask = mask & cfg.kanemaru_platform_mask(a)   # single-platform (Multiome) for Kanemaru
        n = int(mask.sum())
        if n < cfg.MIN_CELLS_PER_GROUP:
            logger.warning("%s: %d cells, skip", label, n); continue
        st = subtype_fracs(linear_mean(a, genes, mask))
        st["group"] = label; st["n"] = n
        rows.append(st)
        logger.info("%s: n=%d", label, n)
    allst = pd.concat(rows, ignore_index=True)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "channel_subtype_conduction.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("conduction-cell channel SUBTYPE fractions, node-matched "
                           "fetal-vs-adult pacemaker; linear scale") + "\n")
        allst[["group", "ion", "subtype", "expr", "fraction", "n"]].round(4).to_csv(
            f, index=False)

    # pivot: subtype fraction per group, for the report
    def frac(group, ion, subtype):
        r = allst[(allst.group == group) & (allst.ion == ion)
                  & (allst.subtype == subtype)]["fraction"]
        return float(r.iloc[0]) if len(r) else np.nan

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "channel_subtype_conduction.md"
    groups_order = [g[0] for g in GROUPS if g[0] in allst.group.unique()]
    with open(rp, "w") as f:
        f.write(cfg.header("conduction-cell channel subtype ratios") + "\n\n")
        f.write("# Channel subtype dominance in CONDUCTION cells — fetal vs adult node\n\n")
        f.write("Subtype-fraction method restricted to PACEMAKER cells, node-matched "
                "(fetal SAN vs adult SAN_P; fetal AVN vs adult AVN_P), with working CM "
                "for contrast. Compares the T-type-leaning of the fetal node Ca program "
                "with the adult node. The working-myocardium T-type switch does not apply "
                "here because nodal cells retain a T-type component.\n\n")

        # --- Ca T-type fraction ---
        f.write("## Calcium: T-type fraction of the Ca-channel program\n\n")
        f.write("| group | n | L-type | T-type | T/(T+L) leaning |\n")
        f.write("|---|---|---|---|---|\n")
        for g in groups_order:
            n = int(allst[allst.group == g]["n"].iloc[0])
            L = frac(g, "Ca", "L-type"); T = frac(g, "Ca", "T-type")
            lean = T / (T + L) if (T + L) > 0 else np.nan
            f.write(f"| {g} | {n} | {L:.3f} | {T:.3f} | {lean:.3f} |\n")
        # explicit fetal-vs-adult node contrast on T-type leaning
        def lean(g):
            L = frac(g, "Ca", "L-type"); T = frac(g, "Ca", "T-type")
            return T / (T + L) if (T + L) > 0 else np.nan
        f.write("\n**Node-matched T-type leaning (T/(T+L)):**\n")
        for node, fg, ag in [("SAN", "fetal_SAN_pace", "adult_SAN_pace"),
                             ("AVN", "fetal_AVN_pace", "adult_AVN_pace")]:
            if fg in groups_order and ag in groups_order:
                fl, al = lean(fg), lean(ag)
                verdict = ("fetal MORE T-type-leaning" if fl > al * 1.1
                           else "fetal LESS T-type-leaning" if fl < al * 0.9
                           else "comparable")
                f.write(f"- {node}: fetal {fl:.3f} vs adult {al:.3f} — {verdict}\n")
        if "adult_work_vCM" in groups_order:
            f.write(f"- (adult working vCM reference: {lean('adult_work_vCM'):.3f} — "
                    "T-type switched off)\n")

        # --- other families, compact ---
        for ion, subs, title in [
            ("K", ["Kv_rapid", "Kir", "Kv_Ito", "Kir_GIRK"],
             "Potassium (IKr/hERG, IK1, Ito, GIRK-vagal)"),
            ("Na", ["Nav_cardiac"], "Sodium (Nav1.5)"),
            ("cation", ["HCN"], "HCN funny current")]:
            f.write(f"\n## {title}\n\n")
            f.write("| group | " + " | ".join(subs) + " |\n")
            f.write("|---" * (len(subs) + 1) + "|\n")
            for g in groups_order:
                vals = [f"{frac(g, ion, s):.3f}" if not np.isnan(frac(g, ion, s)) else "—"
                        for s in subs]
                f.write(f"| {g} | " + " | ".join(vals) + " |\n")

        f.write("\n## Summary\n\n")
        f.write("Pacemaker cells retain a T-type Ca component that working myocardium "
                "lacks (compare nodal groups vs adult_work_vCM). The node-matched "
                "fetal-vs-adult contrast above reports the T-type fraction of the nodal "
                "Ca program in fetal versus adult cells — i.e. the share of the nodal Ca "
                "program contributed by the T-type subtype (Cav3.1, CACNA1G). Read "
                "alongside the per-gene fetal/adult means in the companion "
                "fetal-vs-adult tables (CACNA1G specifically).\n\n")
        f.write("**Caveats:** cross-dataset (fetal Lim/Protze vs adult Kanemaru — "
                "direction-only, different prep/object); single fetal donor per node; "
                "marker-annotated pacemaker clusters (AVN P-cells small/noisier); "
                "T-type absolute levels low so fractions noisier; snRNA-seq.\n")

    logger.info("wrote %s + %s", p1.name, rp.name)


if __name__ == "__main__":
    main()
