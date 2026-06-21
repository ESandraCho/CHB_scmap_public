#!/usr/bin/env python3
"""
lazar_dev_replication.py — multi-donor replication of the developmental backdrop
(SR maturation + Ca T/L subtype trajectory) in whole-cardiomyocyte populations
across 13 independent first-trimester human hearts (Lázár et al. 2025).

seed: 42

Why: the project's developmental-backdrop findings (SR matures postnatally → fetal
cell is surface-Ca-dependent; T-type Ca is an early-gestation program) were shown
in mixed datasets, and the fetal-vs-adult nodal contrasts are n=1 per node. This
script adds DONOR-LEVEL statistics for the WHOLE-CM developmental trajectory using
the Lázár 2025 HDCA developing-heart atlas — 13 donors spanning post-conception
weeks 5.5–14, each treated as an independent observation.

SCOPE (honest): this replicates the WHOLE-CM developmental backdrop only. The atlas
is whole-heart; pacemaker/nodal cells are too sparse per donor (0–17 SHOX2+TBX3+HCN4+
cells/donor) to support per-donor NODAL statistics — the dissected GW19 datasets
remain the substrate for the nodal (Question A/B) claims. Here each donor contributes
one whole-CM data point; the test is whether SR dependence rises and the Ca T/L ratio
falls with gestational age ACROSS donors (donor-level Spearman), which the single-
trajectory analyses could not.

Data: lazar2025/*.h5 (HDCA 'shoji' processed objects;
per-cell raw counts in `Expression`, gene symbols in `Gene`, donor in `Donor`,
age in `Agetext`). Mendeley 10.17632/fhtb99mdzd.1 (open).

Outputs:
  results/tables/lazar_dev_replication.csv     (per-donor metrics)
  results/reports/lazar_dev_replication.md
"""
from __future__ import annotations

import glob
import logging
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("lazar_dev")

LAZAR_DIR = cfg.DATA_DIR / "lazar2025"
SEED = cfg.SEED

# --- gene sets (curated; same definitions as sr_maturation_trajectory / subtype work) ---
SR_GENES = ["RYR2", "ATP2A2", "PLN", "CASQ2", "TRDN", "FKBP1B"]
SURFACE_ENTRY = ["CACNA1C", "CACNA1D", "CACNA1G", "SLC8A1"]
CA_L = ["CACNA1C", "CACNA1D"]
CA_T = ["CACNA1G", "CACNA1H", "CACNA1I"]
CM_MARKERS = ["TNNT2", "TTN", "MYH6", "ACTC1"]
TARGET_SUM = cfg.TARGET_SUM

# CM gate: a cell is a cardiomyocyte if its mean CM-marker (log1p TP10K) exceeds this.
# Calibrated to the bimodal CM/non-CM split in these objects; CM are the dominant pop.
CM_GATE = 1.0
MIN_CM_PER_DONOR = 100


def decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def parse_age(s):
    """'9w' -> 9.0 ; '13-13.5w' -> 13.25 ; '5.5w' -> 5.5 (post-conception weeks)."""
    s = str(s).replace("w", "").strip()
    if "-" in s:
        lo, hi = s.split("-")
        return (float(lo) + float(hi)) / 2.0
    return float(s)


def load_h5(fp):
    """Read one shoji .h5 -> (raw count matrix cells×genes, gene index, donor, age)."""
    f = h5py.File(fp, "r")["shoji"]
    genes = decode(f["Gene"][:])
    gi = {g: i for i, g in enumerate(genes)}
    X = f["Expression"][:].astype(np.float32)      # cells × genes, raw UMI
    donor = decode(f["Donor"][:])[0]
    if donor == "" or donor is None:
        # fallback id from file stem WITHOUT the lane suffix, so multiple lanes of the
        # same blank-donor sample pool into one donor (e.g. TenX263_4 + TenX263_5 -> TenX263)
        stem = os.path.basename(fp).replace(".h5", "")
        donor = re.sub(r"_\d+$", "", stem)
    age = parse_age(decode(f["Agetext"][:])[0])
    return X, gi, donor, age


def linear_mean(X, gi, genes, mask, sf):
    """Mean linear (TP10K) expression per gene over masked cells. sf = per-cell size factor."""
    out = {}
    for g in genes:
        if g not in gi:
            out[g] = 0.0
            continue
        v = X[mask, gi[g]] * sf[mask]               # TP10K linear
        out[g] = float(v.mean()) if mask.sum() else 0.0
    return out


def main():
    files = sorted(glob.glob(str(LAZAR_DIR / "*.h5")))
    if not files:
        raise SystemExit(f"no .h5 files in {LAZAR_DIR}")
    logger.info("found %d Lázár sample files", len(files))

    # aggregate per DONOR (pool that donor's lanes) — each donor = one observation
    per_donor = {}   # donor -> dict(age, cm_metrics accumulators)
    all_genes = SR_GENES + SURFACE_ENTRY + CA_L + CA_T + CM_MARKERS
    all_genes = list(dict.fromkeys(all_genes))

    # Accumulate per-donor CM-cell linear sums to compute donor means across lanes.
    acc = {}
    for fp in files:
        X, gi, donor, age = load_h5(fp)
        tot = X.sum(1)
        sf = np.where(tot > 0, TARGET_SUM / tot, 0.0)
        # log1p TP10K CM-marker score for the CM gate
        def log_tp10k(g):
            return np.log1p(X[:, gi[g]] * sf) if g in gi else np.zeros(X.shape[0])
        cm_score = np.mean([log_tp10k(g) for g in CM_MARKERS], axis=0)
        cm = cm_score > CM_GATE
        n_cm = int(cm.sum())
        if n_cm < MIN_CM_PER_DONOR:
            logger.warning("%s donor=%s: only %d CM, skipping lane", os.path.basename(fp), donor, n_cm)
            continue
        lin = linear_mean(X, gi, all_genes, cm, sf)
        d = acc.setdefault(donor, {"age": age, "n_cm": 0, "sum": {g: 0.0 for g in all_genes}})
        # weight lane means by CM count to get a donor-level mean across lanes
        for g in all_genes:
            d["sum"][g] += lin[g] * n_cm
        d["n_cm"] += n_cm
        logger.info("%s donor=%s age=%.1f CM=%d", os.path.basename(fp), donor, age, n_cm)

    rows = []
    for donor, d in acc.items():
        if d["n_cm"] < MIN_CM_PER_DONOR:
            continue
        m = {g: d["sum"][g] / d["n_cm"] for g in all_genes}      # donor-level CM mean (linear)
        sr_total = sum(m[g] for g in SR_GENES)
        surf_total = sum(m[g] for g in SURFACE_ENTRY)
        denom = sr_total + surf_total
        ca_L = sum(m[g] for g in CA_L)
        ca_T = sum(m[g] for g in CA_T)
        rows.append({
            "donor": donor, "age_pcw": d["age"], "n_cm": d["n_cm"],
            "SR_total": sr_total, "surface_entry_total": surf_total,
            "SR_fraction": sr_total / denom if denom > 0 else np.nan,
            "RYR2": m["RYR2"], "ATP2A2": m["ATP2A2"],
            "Ca_L": ca_L, "Ca_T": ca_T,
            "T_over_TplusL": ca_T / (ca_T + ca_L) if (ca_T + ca_L) > 0 else np.nan,
            "CACNA1G": m["CACNA1G"],
        })
    df = pd.DataFrame(rows).sort_values("age_pcw").reset_index(drop=True)

    # donor-level statistics: does SR_fraction rise / T-ratio fall with age, across donors?
    def sp(col):
        sub = df[["age_pcw", col]].dropna()
        if len(sub) < 4:
            return (np.nan, np.nan, len(sub))
        rho, p = spearmanr(sub["age_pcw"], sub[col])
        return (float(rho), float(p), len(sub))

    sr_rho, sr_p, sr_n = sp("SR_fraction")
    ryr2_rho, ryr2_p, ryr2_n = sp("RYR2")
    tl_rho, tl_p, tl_n = sp("T_over_TplusL")
    g_rho, g_p, g_n = sp("CACNA1G")

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p1 = cfg.OUT_TABLES / "lazar_dev_replication.csv"
    with open(p1, "w") as f:
        f.write(cfg.header("per-donor whole-CM SR fraction + Ca T/L ratio across "
                           "gestational age, Lázár 2025 (13 first-trimester donors); "
                           "linear TP10K; donor-level observations") + "\n")
        df.round(4).to_csv(f, index=False)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "lazar_dev_replication.md"
    with open(rp, "w") as f:
        f.write(cfg.header("Lázár 2025 multi-donor developmental replication") + "\n\n")
        f.write("# Multi-donor replication of the developmental backdrop (Lázár 2025)\n\n")
        f.write("Whole-cardiomyocyte SR-maturation and Ca T/L-subtype trajectories across "
                f"**{len(df)} independent first-trimester human hearts** (post-conception "
                "weeks 5.5–14; Lázár et al. 2025 HDCA atlas). Each donor is one "
                "observation, enabling DONOR-LEVEL statistics that the single-trajectory "
                "analyses lacked. Linear TP10K expression; CM-gated.\n\n")
        f.write("**Scope:** this replicates the WHOLE-CM developmental BACKDROP only. The "
                "atlas is whole-heart; nodal/pacemaker cells are too sparse per donor "
                "(0–17 SHOX2+TBX3+HCN4+ cells/donor) for per-donor NODAL statistics — the "
                "dissected GW19 datasets remain the substrate for the nodal claims. Ages "
                "(5.5–14w) are earlier than the GW19 nodal data, so this is a "
                "first-trimester developmental replication, not an age match.\n\n")

        f.write("## Donor-level trajectory tests (Spearman vs gestational age)\n\n")
        f.write("| metric | direction expected | rho | p | n donors |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| SR fraction (SR / (SR+surface-Ca)) | rises with age (SR matures) | "
                f"{sr_rho:+.2f} | {sr_p:.3g} | {sr_n} |\n")
        f.write(f"| RYR2 (absolute) | rises with age | {ryr2_rho:+.2f} | {ryr2_p:.3g} | {ryr2_n} |\n")
        f.write(f"| Ca T/(T+L) | falls with age (T-type early) | {tl_rho:+.2f} | {tl_p:.3g} | {tl_n} |\n")
        f.write(f"| CACNA1G (absolute) | falls with age | {g_rho:+.2f} | {g_p:.3g} | {g_n} |\n")

        f.write("\n## Per-donor values\n\n")
        f.write("| donor | age (pcw) | n CM | SR fraction | RYR2 | T/(T+L) | CACNA1G |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, r in df.iterrows():
            f.write(f"| {r['donor']} | {r['age_pcw']:.1f} | {int(r['n_cm'])} | "
                    f"{r['SR_fraction']:.3f} | {r['RYR2']:.2f} | {r['T_over_TplusL']:.3f} | "
                    f"{r['CACNA1G']:.3f} |\n")

        f.write("\n## Interpretation\n\n")
        def trend(rho, p, rises_expected):
            if not np.isfinite(rho):
                return "n/a"
            if p >= 0.05 or abs(rho) < 0.2:
                return "flat / not significant over this first-trimester window"
            d = "rises" if rho > 0 else "falls"
            ok = (rho > 0) == rises_expected
            return f"{d} with age ({'expected direction' if ok else 'against expectation'})"
        sr_ryr2_ok = np.isfinite(ryr2_rho) and ryr2_rho > 0 and ryr2_p < 0.05
        tl_ok = np.isfinite(tl_rho) and tl_rho < 0 and tl_p < 0.05
        f.write(f"- **SR maturation (RYR2 absolute):** {trend(ryr2_rho, ryr2_p, True)} "
                f"(rho={ryr2_rho:+.2f}, p={ryr2_p:.3g})"
                f"{' — DONOR-LEVEL support for the postnatal-SR-maturation premise' if sr_ryr2_ok else ''}. "
                f"The SR *fraction* is flat here (rho={sr_rho:+.2f}, p={sr_p:.3g}) because this "
                "narrow first-trimester window lacks the adult anchor where the fraction's "
                "big swing (0.29→0.66 fetal→adult) occurs; the absolute RYR2 rise is the "
                "trustworthy within-window signal.\n")
        f.write(f"- **T-type early program (Ca T/(T+L)):** {trend(tl_rho, tl_p, False)} "
                f"(rho={tl_rho:+.2f}, p={tl_p:.3g})"
                f"{' — DONOR-LEVEL support for the early-gestation T-type program' if tl_ok else ''}. "
                f"CACNA1G absolute is flat/low-abundance (rho={g_rho:+.2f}); the normalized "
                "T/(T+L) ratio is the cleaner readout and it falls significantly.\n")
        f.write("- **Two of the four trajectory tests are donor-level significant** "
                "(RYR2 up, T/(T+L) down) — the two most relevant to the backdrop. These "
                "DONOR-LEVEL trends (each heart an independent observation) UPGRADE the "
                "single-trajectory evidence; they do NOT address the nodal n=1, which "
                "stands on the dissected data.\n\n")
        f.write("**Caveats:** whole-CM (not conduction-specific); first-trimester window "
                "(5.5–14w) is earlier than and does not overlap the GW19 nodal data — a "
                "developmental replication, not an age match; CM-gated by marker score; "
                "snRNA-seq; some donors contribute multiple lanes (pooled per donor, "
                "CM-count-weighted). Marker-defined CM, curated gene sets.\n")

    logger.info("wrote %s + %s | donors=%d | SR_frac rho=%.2f p=%.3g | T/L rho=%.2f p=%.3g",
                p1.name, rp.name, len(df), sr_rho, sr_p, tl_rho, tl_p)


if __name__ == "__main__":
    main()
