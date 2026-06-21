#!/usr/bin/env python3
"""
node_channel_handling_expression.py — node-resolved fetal-vs-adult expression of
cardiac ion channels and the Ca-handling apparatus, by SAN / AVN / vCM.

seed: 42

Computes per-cell-group mean expression (linear TP10K) and detection rate for a
curated panel of Ca / K / Na channels plus the Ca-handling machinery (NCX1, SERCA2a,
RyR2, PLN, CASQ2, ...) implicated in AV-block pathogenesis (CIRCEP.115.003432).

Groups are computed WITHIN each dataset (no cross-platform integration), because
fetal nodal identity exists only in Lim/Protze (cell_type) and adult nodal identity
only in Kanemaru (cell_state); cross-dataset numbers are descriptive and platform is
reported per group:
  - Lim 2024 fetal SAN  : pacemaker_CM / working_CM        (10x snRNA)
  - Protze   fetal AVN  : pacemaker_CM / working_CM        (10x snRNA)
  - Kanemaru adult      : SAN_P / AVN_P / vCM   (Multiome-RNA nuclei ONLY; nodal cells
                          are 100% Multiome nuclei, so vCM is restricted to match)
  - Sim 2021            : fetal / child / adult CM   (snRNA whole-CM ladder, NO node)

Outputs (CSV public, MD report private):
  results/tables/node_channel_handling_expression.csv   (long: group x gene)
  results/reports/node_channel_handling_expression.md   (pivot tables by family)
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
import config_chb_handling_panel as panel

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("node_expr")

GENES = panel.PANEL_GENES


def _group_mask(a, obs_col, matcher):
    """Return a boolean mask of cells in the group, resolving Sim CM-by-stage and
    Kanemaru Multiome-nucleus selectors on the fly from cell_type + stage_bin."""
    if obs_col.startswith("_sim_cm_"):
        stage = obs_col.rsplit("_", 1)[1]            # fetal | child | adult
        return (a.obs["cat"].astype(str) == "CM") & (a.obs["stage_bin"].astype(str) == stage)
    if obs_col.startswith("_kan_multiome_"):
        # Adult pacemaker state (SAN_P / AVN_P), Multiome nuclei.
        which = obs_col.replace("_kan_multiome_", "")    # SAN_P | AVN_P
        mod = (a.obs["modality"].astype(str) == "Multiome-RNA").to_numpy()
        sel = (a.obs["cell_state"].astype(str) == f"{which}_cell").to_numpy()
        return mod & sel
    if obs_col.startswith("_kan_nodework_"):
        # NODE-MATCHED adult working reference: working CMs dissected from the SAME node
        # region as the pacemaker (mirrors the fetal node-adjacent working_CM), Multiome.
        node = obs_col.replace("_kan_nodework_", "")     # SAN | AVN
        mod = (a.obs["modality"].astype(str) == "Multiome-RNA").to_numpy()
        return mod & cfg.kanemaru_node_working_mask(a, node)
    if obs_col == "_kan_bulk_vCM":
        # Separate BULK ventricular reference: vCM* from true ventricular regions
        # (LV/RV/SP/AX), away from the node, Multiome.
        mod = (a.obs["modality"].astype(str) == "Multiome-RNA").to_numpy()
        cs = a.obs["cell_state"].astype(str).str.startswith("vCM").to_numpy()
        return mod & cs & cfg.kanemaru_ventricular_mask(a)
    c = a.obs[obs_col].astype(str)
    kind, val = matcher
    if kind == "prefix":
        return c.str.startswith(val)
    return c.isin(val)


def _expr_stats(a, genes, mask):
    """mean linear TP10K and detection rate (frac cells > 0) for each gene."""
    genes = [g for g in genes if g in a.var_names]
    idx = np.where(mask.to_numpy() if hasattr(mask, "to_numpy") else np.asarray(mask))[0]
    sub = a[idx, genes]
    X = sub.X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lin = np.expm1(X)                       # back to linear TP10K
    mean_lin = lin.mean(axis=0)
    det = (X > 0).mean(axis=0)              # detection fraction (log scale > 0 == count > 0)
    return {g: (float(mean_lin[i]), float(det[i])) for i, g in enumerate(genes)}


def main():
    cache = {}
    rows = []
    for key, disp, fname, obs_col, matcher in panel.GROUPS:
        a = cache.get(fname)
        if a is None:
            a = ad.read_h5ad(cfg.PROC_DIR / fname)
            cache[fname] = a
        mask = _group_mask(a, obs_col, matcher)
        n = int(np.asarray(mask).sum())
        if n < cfg.MIN_CELLS_PER_GROUP:
            logger.warning("%s: only %d cells (<%d), skipping", key, n, cfg.MIN_CELLS_PER_GROUP)
            continue
        stats = _expr_stats(a, GENES, mask)
        for g in GENES:
            if g not in stats:
                continue
            mean_lin, det = stats[g]
            rows.append({
                "group": key, "group_label": disp, "dataset": fname.replace("_qc.h5ad", ""),
                "n_cells": n, "gene": g, "label": panel.LABEL[g],
                "family": panel.FAMILY[g], "mean_tp10k": mean_lin, "detect_frac": det,
            })
        logger.info("%s (%s): n=%d", key, disp, n)

    df = pd.DataFrame(rows)
    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    csv_path = cfg.OUT_TABLES / "node_channel_handling_expression.csv"
    with open(csv_path, "w") as f:
        f.write(cfg.header("node-resolved fetal-vs-adult channel + Ca-handling "
                           "expression by SAN/AVN/vCM; linear TP10K, within-dataset") + "\n")
        df.round(4).to_csv(f, index=False)
    logger.info("wrote %s", csv_path)

    _write_report(df)


def _write_report(df):
    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "node_channel_handling_expression.md"
    # column order = group order as defined, restricted to those present
    group_order = [(k, d) for k, d, *_ in panel.GROUPS if k in set(df.group)]
    n_by_group = df.drop_duplicates("group").set_index("group")["n_cells"].to_dict()

    def pivot(value):
        p = df.pivot_table(index=["family", "label", "gene"], columns="group",
                           values=value, sort=False)
        return p.reindex(columns=[k for k, _ in group_order])

    mean_p = pivot("mean_tp10k")
    det_p = pivot("detect_frac")

    with open(rp, "w") as f:
        f.write(cfg.header("node-resolved channel + Ca-handling expression") + "\n\n")
        f.write("# Cardiac ion channels and Ca-handling apparatus: fetal vs adult, "
                "by SAN / AVN / vCM\n\n")
        f.write("Per-cell-group **mean expression (linear TP10K)** computed WITHIN each "
                "dataset; no cross-platform integration. Fetal nodal identity from Lim 2024 "
                "(SAN, GSE279630) and Protze (AVN, GSE297072) `cell_type`; adult nodal "
                "identity from Kanemaru `cell_state` (SAN_P_cell / AVN_P_cell / vCM*). "
                "Sim 2021 whole-CM (fetal / child / adult) is a NON-nodal developmental "
                "ladder, single platform across the three stages. Cross-dataset numbers are "
                "descriptive (platforms differ: all 10x-family, TP10K-log1p; Lim/Protze "
                "snRNA, Kanemaru mixed sc/sn + multiome, Sim 10x).\n\n")

        f.write("## Cell groups (n cells)\n\n")
        f.write("| group | label | dataset | n |\n|---|---|---|---|\n")
        for k, d in group_order:
            sub = df[df.group == k].iloc[0]
            f.write(f"| {k} | {d} | {sub['dataset']} | {n_by_group[k]} |\n")
        f.write("\n")

        for tbl, title, note in [
            (mean_p, "Mean expression (linear TP10K)",
             "Mean over cells in the group of expm1(log1p TP10K)."),
            (det_p, "Detection rate (fraction of cells with non-zero counts)",
             "Fraction of cells in the group expressing the gene; controls for the "
             "snRNA dropout that depresses mean expression of large transcripts."),
        ]:
            f.write(f"## {title}\n\n{note}\n\n")
            cols = [k for k, _ in group_order]
            hdr = "| family | gene | " + " | ".join(dict(group_order)[c] for c in cols) + " |\n"
            f.write(hdr)
            f.write("|" + "---|" * (2 + len(cols)) + "\n")
            last_fam = None
            for (fam, label, gene), r in tbl.iterrows():
                fam_cell = fam if fam != last_fam else ""
                last_fam = fam
                vals = " | ".join("." if pd.isna(r[c]) else f"{r[c]:.2f}"
                                  if tbl is mean_p else f"{r[c]:.2f}" for c in cols)
                f.write(f"| {fam_cell} | {label} | {vals} |\n")
            f.write("\n")

        _write_gradient_section(f, df)

        f.write("---\n*Gene roles and sources: see `config_chb_handling_panel.py` "
                "(Ca-handling additions per CIRCEP.115.003432).*\n")
    logger.info("wrote %s", rp)


def _write_gradient_section(f, df):
    """Per STAGE BLOCK: which cell type expresses each gene highest, and the
    cell-type gradient (mean TP10K ordered high->low). Within-stage, cross-cell-type
    for fetal/adult; cross-dev-stage only for the Sim block."""
    present = set(df.group)
    label_by_group = dict(zip(df.group, df.group_label))
    mean_by = df.set_index(["group", "gene"])["mean_tp10k"].to_dict()

    f.write("## Highest-expressing cell type per gene (within stage)\n\n")
    f.write("For each gene, the cell type with the **highest mean TP10K** in that stage "
            "block, then the full gradient (group: value, high -> low). Fetal and adult "
            "blocks compare CELL TYPES within one stage and platform; the Sim block is the "
            "only across-DEV-stage gradient (whole-CM, no nodal split).\n\n")

    for block, groups in panel.STAGE_BLOCKS.items():
        groups = [g for g in groups if g in present]
        if not groups:
            continue
        f.write(f"### {block}\n\n")
        f.write("| family | gene | highest | gradient (mean TP10K, high -> low) |\n")
        f.write("|---|---|---|---|\n")
        last_fam = None
        for g in panel.PANEL_GENES:
            vals = [(grp, mean_by.get((grp, g), float("nan"))) for grp in groups]
            vals = [(grp, v) for grp, v in vals if not np.isnan(v)]
            if not vals:
                continue
            vals.sort(key=lambda kv: kv[1], reverse=True)
            top_grp, top_val = vals[0]
            # only call a "winner" if it is meaningfully above the next group
            second = vals[1][1] if len(vals) > 1 else 0.0
            winner = label_by_group[top_grp] if top_val >= 0.05 else "(all ~0)"
            if top_val >= 0.05 and second > 0 and top_val / max(second, 1e-9) < 1.25:
                winner += " ≈"  # near-tie flag
            grad = ", ".join(f"{label_by_group[grp].replace('Fetal ','').replace('Adult ','').replace('Sim ','')}: {v:.2f}"
                             for grp, v in vals)
            fam = panel.FAMILY[g]
            fam_cell = fam if fam != last_fam else ""
            last_fam = fam
            f.write(f"| {fam_cell} | {panel.LABEL[g]} | {winner} | {grad} |\n")
        f.write("\n")


if __name__ == "__main__":
    main()
