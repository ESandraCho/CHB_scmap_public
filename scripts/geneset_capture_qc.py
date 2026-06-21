#!/usr/bin/env python3
"""
geneset_capture_qc.py — method-support QC addressing the criticism that snRNA-seq
under-captures low-abundance genes (electrophysiology genes in particular, which are
lower-abundance than structural sarcomere genes).

Design: the groups that OVERLAP between the snRNA and Smart-seq platforms are WORKING
cardiomyocytes, matched by developmental stage. snRNA vs Smart-seq capture is compared for
three gene sets — electrophysiology (ion channels/transporters), structural (sarcomere),
housekeeping — in matched working-CM groups:

  FETAL  : snRNA  = Lim fetal-SAN working_CM, Protze fetal-AVN working_CM, Sim fetal CM
           Smart-seq = Cui 2019 late-fetal (17-25 wk) ventricular CM (STRT)
  ADULT  : snRNA  = Sim adult CM, Kanemaru ventricular CM
           Smart-seq = Wang 2020 adult ventricular CM (STRT)

Part 1 — gene capture (per gene set, per group):
  (a) DETECTION: % of the gene set detected (>0 in >=1% of cells) — robust to units.
  (b) LEVEL: per-cell expression FRACTION (gene counts / per-cell total), so TPM (Cui),
      UMI (Wang) and counts (snRNA) are unit-matched; report the gene-set median.
  (c) PATTERN: do the gene-set members keep the same relative expression between snRNA and
      Smart-seq — Spearman correlation of per-gene mean fraction, snRNA vs Smart-seq, on the
      shared detected genes (a module-level concordance).

Part 2 — subtype-ratio concordance: for each multi-subtype ion-channel family (Ca L/T,
  Kv rapid/slow/Ito, ...), the within-family subtype RATIO computed on snRNA vs Smart-seq at
  matched stage should agree.

CAVEAT handled in-code: the Kanemaru processed object carries only normalized/log X (no raw
counts), so its per-cell fraction is approximated from expm1(X)/rowsum and flagged; its
detection % is exact.

Env: structural_epitope.
Outputs: results/tables/geneset_capture_qc.json
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
import scipy.sparse as sp

import config_chb as cfg
import config_channel_subtypes as cst

OUT = cfg.OUT_TABLES / "geneset_capture_qc.json"
PROC = cfg.PROC_DIR
RAW = cfg.DATA_DIR
GENESETS = cfg.DATA_DIR / "geneset_cache" / "genesets.json"
DETECT_MIN_FRAC_CELLS = 0.01      # gene "detected" if >0 in >=1% of the group's cells


# ----------------------------------------------------------------- loaders
# Each loader returns a (genes x cells) dense-ish accessor as a dict:
#   {"mat": 2D ndarray genes-in-rows? no -> per-gene helpers}, simplified:
# Returns (var_index: dict gene->row, frac_matrix: cells x genes sparse/array of per-cell
# fractions, n_cells). To keep memory bounded, per-gene stats are computed on the fly per gene set.

def _anndata_group(fname, mask_fn, use_counts=True):
    import anndata as ad
    a = ad.read_h5ad(PROC / fname)
    mask = mask_fn(a)
    if use_counts and "counts" in a.layers:
        X = a.layers["counts"][mask]
        approx = False
    else:
        X = a.X[mask]                      # normalized/log (Kanemaru) -> expm1 approx
        X = sp.csr_matrix(np.expm1(X.toarray())) if sp.issparse(X) else sp.csr_matrix(np.expm1(np.asarray(X)))
        approx = True
    X = sp.csr_matrix(X) if not sp.issparse(X) else X
    var = {g: i for i, g in enumerate(a.var_names)}
    return X, var, int(mask.sum()), approx


def group_snrna(tag):
    if tag == "Lim_fetalSAN":
        return _anndata_group("lim2024_fetal_san_qc.h5ad",
                              lambda a: a.obs["cell_type"].astype(str).values == "working_CM")
    if tag == "Protze_fetalAVN":
        return _anndata_group("lim2024_fetal_avn_qc.h5ad",
                              lambda a: a.obs["cell_type"].astype(str).values == "working_CM")
    if tag == "Sim_fetal":
        return _anndata_group("sim2021_qc.h5ad",
                              lambda a: (a.obs["cat"].astype(str).values == "CM") &
                                        (a.obs["stage_bin"].astype(str).values == "fetal"))
    if tag == "Sim_adult":
        return _anndata_group("sim2021_qc.h5ad",
                              lambda a: (a.obs["cat"].astype(str).values == "CM") &
                                        (a.obs["stage_bin"].astype(str).values == "adult"))
    if tag == "Kanemaru_adult":
        # Restrict to Multiome-RNA single nuclei (cfg.kanemaru_platform_mask) so this QC
        # group is genuinely single-nucleus (the "snRNA" family). NOTE: unlike the nodal
        # FIGURES (which use a node-matched working reference), this QC group is restricted
        # to TRUE VENTRICULAR regions (LV/RV/SP/AX) on purpose — the Smart-seq comparator
        # (Wang 2020) is adult VENTRICULAR CM, so the capture comparison must be ventricle-
        # vs-ventricle, not node-region working CM. Excludes AVN-dissected vCM.
        import config_chb as _cfg
        return _anndata_group("kanemaru_conduction_qc.h5ad",
                              lambda a: (a.obs["cell_type"].astype(str).values == "Ventricular Cardiomyocyte")
                                        & _cfg.kanemaru_platform_mask(a)
                                        & _cfg.kanemaru_ventricular_mask(a),   # ventricle-vs-ventricle vs Wang
                              use_counts=False)
    raise ValueError(tag)


def _df_group_to_X(df):
    """genes(rows) x cells(cols) DataFrame -> (csr cells x genes, var dict, n)."""
    X = sp.csr_matrix(df.to_numpy().T.astype(float))
    var = {g: i for i, g in enumerate(df.index.astype(str))}
    return X, var, df.shape[1], False


def group_cui_late_fetal():
    """Cui STRT, late-fetal (17-25 wk) ventricular CM (LV/RV columns)."""
    p = RAW / "cui2019_gse106118" / "GSE106118_UMI_count_merge.txt.gz"
    with gzip.open(p, "rt") as f:
        df = pd.read_csv(f, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    import re
    # columns are HE{week}W_{donor}_{chamber}.{n}; keep late-fetal (17-25 wk) ventricular (LV/RV)
    keep = [c for c in df.columns
            if (m := re.match(r"^HE(\d+)W_\d+_(LV|RV)", str(c))) and 17 <= int(m.group(1)) <= 25]
    return _df_group_to_X(df[keep])


def group_wang_adult():
    """Wang STRT, adult ventricular CM (N_LV_CM)."""
    base = RAW / "wang2020_gse109816"
    info = pd.read_csv(base / "GSE109816_normal_heart_cell_info.txt.gz", sep="\t", index_col=0)
    ids = set(info.index[info["Type"] == "N_LV_CM"])
    with gzip.open(base / "GSE109816_normal_heart_umi_matrix.csv.gz", "rt") as f:
        df = pd.read_csv(f, sep=",", index_col=0)
    df.index = df.index.astype(str)
    cols = [c for c in df.columns if c in ids]
    return _df_group_to_X(df[cols])


# ----------------------------------------------------------------- metrics
def per_cell_fraction_stats(X, var, genes):
    """For a cells x genes count(-like) matrix, per gene: detection frac and mean per-cell
    fraction (gene / per-cell total). Returns dict gene -> (detect_frac, mean_cellfrac)."""
    tot = np.asarray(X.sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
    out = {}
    for g in genes:
        if g not in var:
            out[g] = (0.0, np.nan); continue
        col = X[:, var[g]]
        col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        detect = float(np.mean(col > 0))
        frac = np.nanmean(col / tot)
        out[g] = (detect, float(frac))
    return out


def geneset_summary(stats, geneset):
    present = [stats[g] for g in geneset if g in stats]
    det = [d for d, _ in present]
    detected = [(g, stats[g]) for g in geneset if g in stats and stats[g][0] >= DETECT_MIN_FRAC_CELLS]
    return {
        "n_in_set": len(geneset),
        "pct_detected": round(100 * np.mean([d >= DETECT_MIN_FRAC_CELLS for d in det]), 1) if det else 0.0,
        "median_detection_frac": round(float(np.median(det)), 3) if det else 0.0,
        "median_cellfrac": float(np.median([f for _, (_, f) in detected if np.isfinite(f)])) if detected else np.nan,
        "_per_gene_frac": {g: s[1] for g, s in detected},   # for cross-platform pattern corr
    }


def pattern_corr(snrna_set, smartseq_set):
    """Spearman of per-gene mean cellfrac, snRNA vs Smart-seq, on shared detected genes."""
    from scipy.stats import spearmanr
    a, b = snrna_set["_per_gene_frac"], smartseq_set["_per_gene_frac"]
    shared = [g for g in a if g in b and np.isfinite(a[g]) and np.isfinite(b[g])]
    if len(shared) < 5:
        return {"rho": None, "n_shared": len(shared)}
    rho, p = spearmanr([a[g] for g in shared], [b[g] for g in shared])
    return {"rho": round(float(rho), 3), "p": float(p), "n_shared": len(shared)}


def cardiac_channel_capture(X, var):
    """Per-channel ABSOLUTE capture for every cardiac-expressed channel subtype: detection
    fraction (% cells >0) and mean per-cell expression fraction. This is the honest
    per-channel comparison — unlike a within-family RATIO, it is not distorted by the
    volatile per-cell-fraction of an abundant denominator channel (e.g. CACNA1C). Each
    channel speaks for itself across platforms."""
    tot = np.asarray(X.sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
    out = {}
    for g in sorted(cst.CARDIAC_SUBTYPE_GENES):
        ion, sub = cst.ion(g), cst.subtype(g)
        if g not in var:
            out[g] = {"ion": ion, "subtype": sub, "detect_pct": 0.0, "cellfrac": np.nan}
            continue
        col = X[:, var[g]]
        col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        out[g] = {"ion": ion, "subtype": sub,
                  "detect_pct": round(100 * float(np.mean(col > 0)), 2),
                  "cellfrac": float(np.nanmean(col / tot))}
    return out


def subtype_ratios(X, var):
    """Within-family subtype RATIO from per-cell-fraction means: for each ion family with
    >1 subtype, the fraction of family signal carried by each subtype. Restricted to
    CARDIAC-relevant channels — trace neuronal/skeletal isoforms (config NON_CARDIAC_*) sit at
    the snRNA detection floor and would make the ratio reflect capture luck, not biology.
    NOTE: retained for completeness but NOT used in the figure — a within-family ratio divides
    by the abundant L-type denominator (CACNA1C), whose per-cell fraction differs markedly
    across platforms, so the ratio reflects the denominator not the subtype's own capture. The
    per-channel absolute capture (cardiac_channel_capture) is the figure metric instead."""
    tot = np.asarray(X.sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
    fam_sub = {}
    for g, (ion, sub, *_ ) in cst.SUBTYPES.items():
        if g not in var or ion is None or not cst.is_cardiac_subtype(g):
            continue
        col = X[:, var[g]]
        col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        val = float(np.nanmean(col / tot))
        fam_sub.setdefault(ion, {}).setdefault(sub, 0.0)
        fam_sub[ion][sub] += val
    # ratios within each family
    ratios = {}
    for ion, subs in fam_sub.items():
        s = sum(subs.values())
        if s > 0 and len(subs) > 1:
            ratios[ion] = {sub: round(v / s, 4) for sub, v in subs.items()}
    return ratios


# ----------------------------------------------------------------- driver
def main():
    gs = json.loads(GENESETS.read_text())
    sets = gs["genesets"]
    all_set_genes = sorted({g for v in sets.values() for g in v} |
                           set(cst.SUBTYPES.keys()))

    groups = {
        "fetal": {
            "snRNA": ["Lim_fetalSAN", "Protze_fetalAVN", "Sim_fetal"],
            "smartseq": "Cui_late_fetal",
        },
        "adult": {
            "snRNA": ["Sim_adult", "Kanemaru_adult"],
            "smartseq": "Wang_adult",
        },
    }
    loaders = {
        "Cui_late_fetal": group_cui_late_fetal,
        "Wang_adult": group_wang_adult,
    }

    result = {"detection_min_frac_cells": DETECT_MIN_FRAC_CELLS,
              "geneset_sources": gs["sources"], "geneset_counts": gs["counts"],
              "stages": {}}

    for stage, spec in groups.items():
        print(f"  [{stage}] loading groups ...")
        loaded = {}
        for tag in spec["snRNA"]:
            X, var, n, approx = group_snrna(tag)
            loaded[tag] = (X, var, n, approx, "snRNA")
        sm = spec["smartseq"]
        X, var, n, approx = loaders[sm]()
        loaded[sm] = (X, var, n, approx, "Smart-seq")

        # part 1: gene-set capture per group
        per_group = {}
        setstats = {}
        for tag, (X, var, n, approx, plat) in loaded.items():
            stats = per_cell_fraction_stats(X, var, all_set_genes)
            setstats[tag] = {name: geneset_summary(stats, sets[name]) for name in sets}
            per_group[tag] = {"platform": plat, "n_cells": n, "counts_approx": approx,
                              "genesets": {name: {k: v for k, v in setstats[tag][name].items()
                                                  if not k.startswith("_")} for name in sets}}

        # part 1c: snRNA-vs-Smartseq pattern correlation per gene set (each snRNA group vs the
        # Smart-seq group of the same stage)
        patt = {}
        for tag in spec["snRNA"]:
            patt[tag] = {name: pattern_corr(setstats[tag][name], setstats[sm][name]) for name in sets}

        # part 2: per-channel ABSOLUTE capture (figure metric) + subtype ratios (completeness)
        chan_cap = {tag: cardiac_channel_capture(loaded[tag][0], loaded[tag][1]) for tag in loaded}
        subr = {tag: subtype_ratios(loaded[tag][0], loaded[tag][1]) for tag in loaded}

        result["stages"][stage] = {
            "groups": per_group,
            "snrna_vs_smartseq_pattern_rho": patt,
            "cardiac_channel_capture": chan_cap,
            "subtype_ratios": subr,
            "smartseq_ref": sm,
        }
        print(f"  [{stage}] done: {[ (t, per_group[t]['n_cells']) for t in per_group ]}")

    # nodal pacemaker-vs-working per-channel capture in the dissected fetal nodes (all snRNA):
    # shows the F1D-relevant channels (Cav3.1/Cav1.3/HCN4) are detected AND pacemaker-enriched,
    # i.e. the within-platform contrasts are not a capture floor.
    print("  [nodal] pacemaker-vs-working per-channel capture ...")
    import anndata as ad
    nodal = {}
    for fname, lab in [("lim2024_fetal_san_qc.h5ad", "fetal_SAN"),
                       ("lim2024_fetal_avn_qc.h5ad", "fetal_AVN")]:
        a = ad.read_h5ad(PROC / fname)
        X = a.layers["counts"] if "counts" in a.layers else a.X
        var = {g: i for i, g in enumerate(a.var_names)}
        ct = a.obs["cell_type"].astype(str).values
        grp = {}
        for cell in ("pacemaker_CM", "working_CM"):
            m = ct == cell
            if m.sum() >= 20:
                grp[cell] = {"n_cells": int(m.sum()),
                             "channels": cardiac_channel_capture(sp.csr_matrix(X)[m], var)}
        nodal[lab] = grp
    result["nodal_pacemaker_vs_working"] = nodal

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2))
    print(f"\nWrote {OUT.name}")


if __name__ == "__main__":
    main()
