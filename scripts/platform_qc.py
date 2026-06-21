#!/usr/bin/env python3
"""
platform_qc.py — platform-suitability QC for the CHB analyses, addressing the critique
that droplet 10x under-samples large cardiomyocytes (cell-size / RNA-content bias).

Two honest defenses, quantified:

  (A) Per-cell transcriptome complexity by platform. The Smart-seq-family (STRT) datasets
      (Cui 2019, Wang 2020) recover MORE genes per cell than droplet 10x — the expected,
      and favourable, contrast — so the cross-platform corroboration of the developmental
      direction (S3) is not a case of one weak platform agreeing with another. Reported as
      median genes/cell per dataset.

  (B) Depth-robustness of the nodal within-donor contrast. The headline nodal result
      (pacemaker vs working CM, Fig 1A/B) is computed WITHIN one donor and one chemistry,
      so any 10x size/depth bias is shared by both populations and largely cancels. This is
      tested directly: recompute the pacemaker/working log2 fold-change for the anti-Ro targets
      on depth-matched cells (pacemaker and working cells downsampled to a common
      total-count quantile range) and confirm the sign/direction is unchanged. If the result
      survived only because pacemaker and working cells differ in depth, depth-matching would
      remove it; it does not.

NOTE ON SCOPE: there is NO Smart-seq dataset of dissected human SAN/AVN conduction tissue in
existence; pacemaker cells are resolvable only in droplet single-nucleus data, and the ADULT
node specifically only in 10x Multiome nuclei (the SAN/AVN regions were Multiome-sequenced;
SAN_P/AVN_P have zero 3'-snRNA cells). The Smart-seq corroboration (Cui/Wang) therefore
supports the WHOLE-CM developmental arm, not the nodal arm; the nodal arm's platform defense
is (B) plus the use of single-nucleus (size-bias-robust) data, not cross-platform replication.
Multiome is reported in the complexity panel (A) for completeness — it is droplet snRNA-family,
not an orthogonal technology, so it does not serve as cross-platform corroboration.

Env: structural_epitope.
Outputs:
  results/tables/platform_qc.json
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config_chb as cfg

OUT = cfg.OUT_TABLES / "platform_qc.json"
PROC = cfg.PROC_DIR
RAW = cfg.DATA_DIR
TARGETS = list(cfg.TARGET_PHENOTYPE)


def _genes_per_cell_10x(fname, cm_only_col=None):
    """median genes/cell and counts/cell for a processed 10x object."""
    import anndata as ad
    a = ad.read_h5ad(PROC / fname)
    gcol = "n_genes_by_counts" if "n_genes_by_counts" in a.obs else "n_genes"
    g = a.obs[gcol].to_numpy(float)
    c = a.obs["total_counts"].to_numpy(float) if "total_counts" in a.obs else np.full(len(g), np.nan)
    return {"median_genes_per_cell": float(np.median(g)),
            "median_counts_per_cell": float(np.median(c)),
            "n_cells": int(len(g))}


def qc_10x_datasets():
    out = {}
    for tag, fname in [("Sim2021_10x", "sim2021_qc.h5ad"),
                       ("Lim_fetalSAN_10x", "lim2024_fetal_san_qc.h5ad"),
                       ("Protze_fetalAVN_10x", "lim2024_fetal_avn_qc.h5ad")]:
        try:
            out[tag] = {**_genes_per_cell_10x(fname), "platform": "10x snRNA"}
        except Exception as e:
            out[tag] = {"error": str(e)}
    return out


def qc_kanemaru_multiome_nodal():
    """Per-cell complexity of the adult NODAL cells actually used (Kanemaru SAN_P/AVN_P),
    which are 100% 10x Multiome single nuclei. Reported so the complexity panel (S1) shows
    where the Multiome platform carrying the headline nodal result sits — NOT as a
    cross-platform corroboration (Multiome is droplet snRNA-family, not an orthogonal
    technology; the orthogonal Smart-seq corroboration is S3)."""
    import anndata as ad
    try:
        a = ad.read_h5ad(PROC / "kanemaru_conduction_qc.h5ad")
    except Exception as e:
        return {"error": str(e)}
    cs = a.obs["cell_state"].astype(str)
    mod = a.obs["modality"].astype(str)
    mask = mod.eq("Multiome-RNA").to_numpy() & cs.isin(["SAN_P_cell", "AVN_P_cell"]).to_numpy()
    gcol = "n_genes_by_counts" if "n_genes_by_counts" in a.obs else "n_genes"
    g = a.obs[gcol].to_numpy(float)[mask]
    c = (a.obs["total_counts"].to_numpy(float)[mask]
         if "total_counts" in a.obs else np.full(g.shape, np.nan))
    if g.size == 0:
        return {"error": "no Multiome nodal cells"}
    return {"median_genes_per_cell": float(np.median(g)),
            "median_counts_per_cell": float(np.median(c)),
            "n_cells": int(g.size), "platform": "10x Multiome (snRNA-family)"}


def qc_wang_strt():
    """Wang 2020 STRT: genes/cell from the cell_info table (Distinct.Genes...)."""
    base = RAW / "wang2020_gse109816"
    info = pd.read_csv(base / "GSE109816_normal_heart_cell_info.txt.gz", sep="\t", index_col=0)
    cm = info[info["Type"].isin(["N_LA_CM", "N_LV_CM"])]
    gcol = [c for c in info.columns if "Distinct.Genes" in c][0]
    ucol = [c for c in info.columns if c == "Distinct.UMIs"]
    res = {"platform": "STRT (non-10x)",
           "median_genes_per_cell": float(cm[gcol].median()),
           "n_cells": int(len(cm))}
    if ucol:
        res["median_counts_per_cell"] = float(cm[ucol[0]].median())
    return res


def qc_cui_strt():
    """Cui 2019 STRT: TPM matrix -> genes/cell = non-zero TPM entries per column."""
    p = RAW / "cui2019_gse106118" / "GSE106118_UMI_count_merge.txt.gz"
    with gzip.open(p, "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    genes_per_cell = (df > 0).sum(axis=0)
    return {"platform": "STRT (non-10x)",
            "median_genes_per_cell": float(genes_per_cell.median()),
            "n_cells": int(df.shape[1])}


def nodal_depth_robustness():
    """For each fetal node, recompute the pacemaker/working log2 FC of the anti-Ro
    targets on DEPTH-MATCHED cells (both populations restricted to the overlapping
    central total-count band), and compare to the full-data FC. Direction-stable =
    the contrast is not a depth artifact."""
    import anndata as ad
    res = {}
    for fname, lab in [("lim2024_fetal_san_qc.h5ad", "fetal_SAN"),
                       ("lim2024_fetal_avn_qc.h5ad", "fetal_AVN")]:
        a = ad.read_h5ad(PROC / fname)
        counts = a.layers["counts"] if "counts" in a.layers else a.X
        tot = np.asarray(counts.sum(1)).ravel().astype(float)
        ct = a.obs["cell_type"].astype(str).values
        pace = ct == "pacemaker_CM"
        work = ct == "working_CM"

        # depth-match: keep cells in the overlap of the two populations' 10-90 pctile depth
        lo = max(np.percentile(tot[pace], 10), np.percentile(tot[work], 10))
        hi = min(np.percentile(tot[pace], 90), np.percentile(tot[work], 90))
        band = (tot >= lo) & (tot <= hi)

        def log2fc(mask_a, mask_b, gi):
            va = _tp10k_vec(counts, tot, gi, mask_a)
            vb = _tp10k_vec(counts, tot, gi, mask_b)
            return float(np.log2((np.nanmean(va) + cfg.EPS) / (np.nanmean(vb) + cfg.EPS)))

        genes = ["CACNA1D", "CACNA1G", "KCNH2", "SCN5A", "HCN4"]
        full, matched = {}, {}
        for g in genes:
            if g not in a.var_names:
                continue
            gi = a.var_names.get_loc(g)
            full[g] = round(log2fc(pace, work, gi), 3)
            matched[g] = round(log2fc(pace & band, work & band, gi), 3)
        res[lab] = {
            "depth_band_counts": [float(lo), float(hi)],
            "n_pace_full": int(pace.sum()), "n_work_full": int(work.sum()),
            "n_pace_matched": int((pace & band).sum()), "n_work_matched": int((work & band).sum()),
            "log2fc_full": full, "log2fc_depth_matched": matched,
            "direction_stable": {g: bool(np.sign(full[g]) == np.sign(matched[g])) for g in full},
        }
    return res


def _tp10k_vec(counts, tot, gi, mask):
    v = counts[mask][:, gi]
    v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
    t = tot[mask].copy(); t[t == 0] = np.nan
    return v / t * cfg.TARGET_SUM


def main():
    print("  10x datasets complexity ...")
    tenx = qc_10x_datasets()
    print("  Wang STRT complexity ...")
    wang = qc_wang_strt()
    print("  Cui STRT complexity ...")
    cui = qc_cui_strt()
    print("  Kanemaru Multiome nodal complexity ...")
    kan_nodal = qc_kanemaru_multiome_nodal()
    print("  nodal depth-robustness ...")
    nodal = nodal_depth_robustness()

    complexity = {**tenx, "Kanemaru_adultNodal_Multiome": kan_nodal,
                  "Wang_adult_STRT": wang, "Cui_fetal_STRT": cui}
    all_stable = all(all(v["direction_stable"].values()) for v in nodal.values())

    out = {
        "purpose": "Platform-suitability QC vs the '10x under-samples cardiomyocytes' critique.",
        "scope_note": "No Smart-seq/snRNA dataset of dissected human SAN/AVN exists; pacemaker "
                      "cells are resolvable only in droplet single-nucleus data, and the adult "
                      "node specifically only in 10x MULTIOME nuclei (the SAN/AVN regions were "
                      "Multiome-sequenced; SAN_P/AVN_P have zero 3'-snRNA cells). Cui/Wang STRT "
                      "corroborate the WHOLE-CM developmental arm; the nodal arm rests on "
                      "single-nucleus data (size-bias-robust) + the within-donor depth-robustness "
                      "test below. NOTE: Multiome is droplet snRNA-family, NOT an orthogonal "
                      "platform — it is reported in the complexity panel for completeness, not as "
                      "cross-platform corroboration (that role is Smart-seq, S3).",
        "A_per_cell_complexity": complexity,
        "B_nodal_depth_robustness": nodal,
        "nodal_direction_stable_all": bool(all_stable),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {OUT.name}")
    print("  per-cell genes (median):",
          {k: v.get("median_genes_per_cell") for k, v in complexity.items()})
    print("  nodal direction stable after depth-matching:", all_stable)


if __name__ == "__main__":
    main()
