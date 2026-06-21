#!/usr/bin/env python3
"""
cav_subtype_pace_work_stats.py — pacemaker-vs-working statistics for the Cav-subtype
COMPOSITION shown in Fig 1B, per node x stage pair.

seed: 42

Question (Sandra): within each node x stage group, is the per-cell FRACTION of the
Cav pool contributed by CACNA1D (Cav1.3) and by CACNA1G (Cav3.1) significantly HIGHER
in pacemaker cells than in their matched (same node, same stage) working cardiomyocytes?

Design:
  - Cav pool = {CACNA1C, CACNA1D, CACNA1G, CACNA1H} (the four Fig-1B subtypes), linear TP10K.
  - Per cell with Cav-pool > 0: fraction_i = Cav_i / sum(Cav pool).
  - Per subtype, per pair: one-sided Mann-Whitney U (pacemaker > working) on the per-cell
    fractions — matching the per-cell Mann-Whitney convention of Fig 1D/E (fetal/adult
    pacemaker-vs-working), here on the composition fraction rather than raw expression.
  - BH-FDR across all (pair x subtype) tests.

Pairs (pacemaker vs node-matched working): fetal SAN, fetal AVN, adult SAN, adult AVN.
The two headline channels are CACNA1D and CACNA1G; CACNA1C/CACNA1H are reported for context.

Caveat (same as Fig 1): the fetal groups are n=1 donor, so cells are pseudoreplicated —
the p-values are within-sample (cell-level), descriptive, not donor-level inference.

Outputs:
  results/tables/cav_subtype_pace_work_stats.csv
  results/reports/cav_subtype_pace_work_stats.md
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
from node_channel_handling_expression import _group_mask

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("cav_pw_stats")

CAV4 = ["CACNA1C", "CACNA1D", "CACNA1G", "CACNA1H"]
LABEL = {"CACNA1C": "Cav1.2 (L)", "CACNA1D": "Cav1.3 (L)",
         "CACNA1G": "Cav3.1 (T)", "CACNA1H": "Cav3.2 (T)"}
HEADLINE = {"CACNA1D", "CACNA1G"}

# (pair_label, node, stage, file, pace obs_col/matcher, work obs_col/matcher)
PAIRS = [
    ("fetal SAN", "lim2024_fetal_san_qc.h5ad",
     ("cell_type", ("isin", {"pacemaker_CM"})), ("cell_type", ("isin", {"working_CM"}))),
    ("fetal AVN", "lim2024_fetal_avn_qc.h5ad",
     ("cell_type", ("isin", {"pacemaker_CM"})), ("cell_type", ("isin", {"working_CM"}))),
    ("adult SAN", "kanemaru_conduction_qc.h5ad",
     ("_kan_multiome_SAN_P", ("isin", {"__TRUE__"})), ("_kan_nodework_SAN", ("isin", {"__TRUE__"}))),
    ("adult AVN", "kanemaru_conduction_qc.h5ad",
     ("_kan_multiome_AVN_P", ("isin", {"__TRUE__"})), ("_kan_nodework_AVN", ("isin", {"__TRUE__"}))),
]


def _percell_cav_fractions(a, mask):
    """cells-in-group x 4 array of per-cell Cav-subtype fractions (Cav_i / Cav-pool),
    keeping only cells whose Cav pool > 0. Returns (frac_array, n_cells_group, n_qualifying)."""
    genes = [g for g in CAV4 if g in a.var_names]
    idx = np.where(np.asarray(mask))[0]
    X = a[idx, genes].X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    lin = np.expm1(X)
    pool = lin.sum(1)
    ok = pool > 0
    frac = np.full((lin.shape[0], len(CAV4)), np.nan)
    col = {g: i for i, g in enumerate(CAV4)}
    for j, g in enumerate(genes):
        frac[ok, col[g]] = lin[ok, j] / pool[ok]
    return frac, int(len(idx)), int(ok.sum())


def main():
    cache = {}
    rows = []
    for pair, fname, (pcol, pmatch), (wcol, wmatch) in PAIRS:
        a = cache.get(fname)
        if a is None:
            a = ad.read_h5ad(cfg.PROC_DIR / fname); cache[fname] = a
        pmask = _group_mask(a, pcol, pmatch)
        wmask = _group_mask(a, wcol, wmatch)
        pf, n_p, nq_p = _percell_cav_fractions(a, pmask)
        wf, n_w, nq_w = _percell_cav_fractions(a, wmask)
        logger.info("%s: pacemaker n=%d (Cav+ %d), working n=%d (Cav+ %d)",
                    pair, n_p, nq_p, n_w, nq_w)
        for j, g in enumerate(CAV4):
            p_vals = pf[:, j][~np.isnan(pf[:, j])]
            w_vals = wf[:, j][~np.isnan(wf[:, j])]
            mean_p = float(p_vals.mean()) if p_vals.size else np.nan
            mean_w = float(w_vals.mean()) if w_vals.size else np.nan
            # one-sided Mann-Whitney: pacemaker fraction > working fraction
            if p_vals.size >= 3 and w_vals.size >= 3 and (p_vals.var() + w_vals.var()) > 0:
                pval = float(mannwhitneyu(p_vals, w_vals, alternative="greater").pvalue)
            else:
                pval = np.nan
            rows.append({
                "pair": pair, "gene": g, "subtype": LABEL[g], "headline": g in HEADLINE,
                "n_pace_cavpos": int(p_vals.size), "n_work_cavpos": int(w_vals.size),
                "mean_frac_pace": mean_p, "mean_frac_work": mean_w,
                "delta_pace_minus_work": (mean_p - mean_w) if np.isfinite(mean_p) and np.isfinite(mean_w) else np.nan,
                "p_one_sided_pace_gt_work": pval,
            })

    df = pd.DataFrame(rows)
    ok = df["p_one_sided_pace_gt_work"].notna().to_numpy()
    fdr = np.full(len(df), np.nan)
    if ok.any():
        fdr[ok] = false_discovery_control(df.loc[ok, "p_one_sided_pace_gt_work"].to_numpy())
    df["fdr_bh"] = fdr
    df["sig_fdr_0.05"] = df["fdr_bh"] < 0.05

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    csv = cfg.OUT_TABLES / "cav_subtype_pace_work_stats.csv"
    with open(csv, "w") as f:
        f.write(cfg.header("Cav-subtype composition pacemaker-vs-working stats (Fig 1B): "
                           "per-cell Cav-pool fraction, one-sided MWU pace>work, BH-FDR") + "\n")
        df.round(5).to_csv(f, index=False)
    logger.info("wrote %s", csv)

    _write_report(df)


def _write_report(df):
    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "cav_subtype_pace_work_stats.md"
    with open(rp, "w") as f:
        f.write(cfg.header("Cav-subtype pacemaker-vs-working composition stats") + "\n\n")
        f.write("# Is the Cav1.3 / Cav3.1 share of the Cav pool higher in pacemaker than "
                "matched working cardiomyocytes? (Fig 1B)\n\n")
        f.write("Per node x stage pair, per Cav subtype: per-cell fraction of the four-Cav pool "
                "(Cav_i / [Cav1.2+Cav1.3+Cav3.1+Cav3.2], cells with Cav-pool > 0). One-sided "
                "Mann-Whitney U (pacemaker > working); BH-FDR across all pair x subtype tests. "
                "Headline channels **CACNA1D (Cav1.3)** and **CACNA1G (Cav3.1)** in bold. The "
                "fetal pairs are n=1 donor (cells pseudoreplicated): p-values are within-sample "
                "cell-level, descriptive, not donor-level inference.\n\n")
        f.write("| pair | subtype | mean frac pace | mean frac work | Δ (pace−work) | "
                "p (1-sided, pace>work) | BH-FDR | sig |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for _, r in df.iterrows():
            name = f"**{r['subtype']}**" if r["headline"] else r["subtype"]
            sig = "✓" if r["sig_fdr_0.05"] else ""
            def fmt(x): return "n/a" if pd.isna(x) else (f"{x:.2e}" if x < 1e-3 else f"{x:.4f}")
            f.write(f"| {r['pair']} | {name} | {r['mean_frac_pace']:.3f} | "
                    f"{r['mean_frac_work']:.3f} | {r['delta_pace_minus_work']:+.3f} | "
                    f"{fmt(r['p_one_sided_pace_gt_work'])} | {fmt(r['fdr_bh'])} | {sig} |\n")
        f.write("\n*Source: `cav_subtype_pace_work_stats.py`. ✓ = BH-FDR < 0.05.*\n")
    logger.info("wrote %s", rp)


if __name__ == "__main__":
    main()
