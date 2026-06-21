#!/usr/bin/env python3
"""
node_datasets_methods_table.py — methods-section datasets & cell-types table for the
node-resolved channel + Ca-handling expression analysis.

seed: 42

Emits two tables, with cell counts computed LIVE from the processed objects (never
hard-coded) so the methods text cannot drift from the data actually analysed:

  1. DATASET table  — one row per dataset: accession, stage, chemistry, the platform
     ACTUALLY used by this analysis, citation, and live total CM/used-cell counts.
  2. CELL-GROUP table — one row per analysis group (the GROUPS in the panel config):
     dataset, cell-type label / obs selector, platform, and live n cells.

Provenance fields come from config_chb_handling_panel.DATASET_PROVENANCE (transcribed
from DATA_SOURCES.md, the authoritative source). Counts come from the objects.

Outputs:
  results/tables/node_datasets_methods.csv          (cell-group rows)
  results/tables/node_datasets_methods_summary.csv  (dataset rows)
  results/reports/node_datasets_methods.md          (both, formatted for methods)
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg
import config_chb_handling_panel as panel
from node_channel_handling_expression import _group_mask

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("methods_tbl")

# selector text shown in the methods table (human-readable obs criterion per group)
SELECTOR_TEXT = {
    "fetal_SAN_pace": "cell_type == pacemaker_CM",
    "fetal_SAN_work": "cell_type == working_CM",
    "fetal_AVN_pace": "cell_type == pacemaker_CM",
    "fetal_AVN_work": "cell_type == working_CM",
    "adult_SAN_pace": "cell_state == SAN_P_cell, modality == Multiome-RNA",
    "adult_SAN_work": "working CM (vCM*/aCM*), region == SAN, modality == Multiome-RNA",
    "adult_AVN_pace": "cell_state == AVN_P_cell, modality == Multiome-RNA",
    "adult_AVN_work": "working CM (vCM*/aCM*), region == AVN, modality == Multiome-RNA",
    "adult_vCM":      "cell_state startswith vCM, region in {LV,RV,SP,AX}, modality == Multiome-RNA",
    "sim_fetal_CM":   "cat == CM & stage_bin == fetal",
    "sim_child_CM":   "cat == CM & stage_bin == child",
    "sim_adult_CM":   "cat == CM & stage_bin == adult",
}


def main():
    cache = {}
    group_rows = []
    used_by_file = {}
    for key, disp, fname, obs_col, matcher in panel.GROUPS:
        a = cache.get(fname)
        if a is None:
            a = ad.read_h5ad(cfg.PROC_DIR / fname)
            cache[fname] = a
        n = int(np.asarray(_group_mask(a, obs_col, matcher)).sum())
        used_by_file[fname] = used_by_file.get(fname, 0) + n
        prov = panel.DATASET_PROVENANCE[fname]
        group_rows.append({
            "group": key, "group_label": disp, "dataset": prov["name"],
            "object_file": fname, "selector": SELECTOR_TEXT.get(key, ""),
            "platform_used": prov["platform_used"], "n_cells": n,
        })
        logger.info("%s: n=%d", key, n)

    gdf = pd.DataFrame(group_rows)

    # dataset-level summary — only the objects this analysis actually loads (the four
    # h5ad-backed node/dev datasets in panel.GROUPS), for which object/used counts are
    # computed live. The non-h5ad corroboration datasets (Cui/Wang/Lazar) carry no counts
    # here and are covered by supp_dataset_analysis_table.py / build_supplementary_tables.py.
    ds_rows = []
    for fname in cache:
        prov = panel.DATASET_PROVENANCE[fname]
        ds_rows.append({
            "dataset": prov["name"], "object_file": fname, "accession": prov["accession"],
            "tissue_stage": prov["tissue_stage"], "chemistry": prov["chemistry"],
            "platform_used": prov["platform_used"], "donors": prov["donors"],
            "citation": prov["citation"],
            "n_cells_object": int(a.n_obs) if (a := cache[fname]) is not None else 0,
            "n_cells_used": used_by_file.get(fname, 0),
        })
    ddf = pd.DataFrame(ds_rows)

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p_groups = cfg.OUT_TABLES / "node_datasets_methods.csv"
    p_summary = cfg.OUT_TABLES / "node_datasets_methods_summary.csv"
    hdr = cfg.header("datasets & cell-types used by node_channel_handling_expression; "
                     "counts live from objects, provenance from DATA_SOURCES.md")
    for path, frame in [(p_groups, gdf), (p_summary, ddf)]:
        with open(path, "w") as f:
            f.write(hdr + "\n")
            frame.to_csv(f, index=False)
        logger.info("wrote %s", path)

    _write_report(ddf, gdf)


def _write_report(ddf, gdf):
    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "node_datasets_methods.md"
    with open(rp, "w") as f:
        f.write(cfg.header("datasets & cell-types — methods table") + "\n\n")
        f.write("# Datasets and cell types used (methods)\n\n")
        f.write("Single-cell/-nucleus datasets and the cell groups drawn from each for the "
                "node-resolved fetal-vs-adult channel + Ca-handling expression analysis. "
                "All objects are TP10K-normalized, log1p; means are computed on linear "
                "(expm1) values. Cell counts are computed live from the objects; provenance "
                "(accession / stage / chemistry / citation) is transcribed from "
                "`DATA_SOURCES.md`.\n\n")

        f.write("## Datasets\n\n")
        f.write("| Dataset | Accession | Tissue / stage | Chemistry (object) | "
                "Platform USED here | Donors | n cells (object) | n cells (used) | Citation |\n")
        f.write("|" + "---|" * 9 + "\n")
        for _, r in ddf.iterrows():
            f.write(f"| {r['dataset']} | {r['accession']} | {r['tissue_stage']} | "
                    f"{r['chemistry']} | {r['platform_used']} | {r['donors']} | "
                    f"{int(r['n_cells_object'])} | {int(r['n_cells_used'])} | {r['citation']} |\n")
        f.write("\n")

        f.write("## Cell groups (analysis units)\n\n")
        f.write("Each row is one group compared in the analysis. Fetal and adult blocks are "
                "WITHIN-stage cross-CELL-TYPE; the Sim block is the only cross-dev-stage "
                "comparison (whole-CM, no nodal resolution).\n\n")
        f.write("| Block | Group | Dataset | Cell-type selector | Platform | n cells |\n")
        f.write("|---|---|---|---|---|---|\n")
        block_of = {g: b for b, gs in panel.STAGE_BLOCKS.items() for g in gs}
        for _, r in gdf.iterrows():
            f.write(f"| {block_of.get(r['group'], '')} | {r['group_label']} | {r['dataset']} | "
                    f"`{r['selector']}` | {r['platform_used'].split(';')[0].split(':')[0]} | "
                    f"{int(r['n_cells'])} |\n")
        f.write("\n")

        f.write("## Platform-consistency note (for methods/limitations)\n\n")
        f.write("The adult conduction object (Kanemaru) is a multi-platform atlas (10x snRNA "
                "+ Multiome-RNA + scRNA). The dissected SAN and AVN regions were sequenced "
                "**only by Multiome-RNA single nuclei**, so the adult pacemaker states "
                "(`SAN_P_cell`, `AVN_P_cell`) exist only in the Multiome data (0 snRNA cells). "
                "To keep the within-adult node-vs-working comparison on a single platform, the "
                "adult ventricular working group (`vCM`) is therefore RESTRICTED to Multiome-RNA "
                "nuclei as well. snRNA cannot be used for the adult node because the node is not "
                "represented in the snRNA fraction.\n")
    logger.info("wrote %s", rp)


if __name__ == "__main__":
    main()
