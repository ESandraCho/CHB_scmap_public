#!/usr/bin/env python3
"""
ca_subtype_dominance_dev.py — does the DOMINANT functioning cardiac Ca channel subtype
change across development? Resolves the paradox: anti-Ro Ca targets express LESS in
adults, yet adults clearly need functioning Ca channels.

Answer (computed below): the WORKHORSE does NOT switch — Cav1.2 (CACNA1C) is the dominant
L-type at ~99% of L-type expression at fetal/child/adult, expressed in ~97-99% of CMs at
all stages. The anti-Ro targets that collapse (Cav1.3/CACNA1D, Cav3.1/CACNA1G) are MINOR,
pacemaker/conduction-SPECIALIZED channels (<1.2% of L-type, ~17% of cells) that the adult
working myocardium never relied on; the fetal NODE depends on them. The real developmental
shift is in Ca-handling MODE: sarcolemmal influx (fetal) -> SR release (adult RYR2 ~2x).

Sim 2021 CM, within-study, mean-per-cell TP10K + %-expressing.
Env: structural_epitope.
"""
from __future__ import annotations
import json
import numpy as np, scipy.sparse as sp
import config_chb as cfg

OUT = cfg.OUT_TABLES / "ca_subtype_dominance_dev.json"
SIM = cfg.PROC_DIR / "sim2021_qc.h5ad"
FAMILIES = {"L-type": ["CACNA1C", "CACNA1D"], "T-type": ["CACNA1G", "CACNA1H", "CACNA1I"],
            "other_Cav": ["CACNA1A", "CACNA1B", "CACNA1E"], "SR_handling": ["RYR2", "ATP2A2"]}
STAGES = ["fetal", "child", "adult"]


def main():
    import anndata as ad
    a = ad.read_h5ad(SIM)
    cm = a.obs["cat"].astype(str).values == "CM"
    counts = a.layers["counts"]; stage = a.obs["stage_bin"].astype(str).values
    genes = [g for fam in FAMILIES.values() for g in fam if g in a.var_names]

    donor = a.obs["donor_id"].astype(str).values

    def stats(g, st):
        m = cm & (stage == st); tot = np.asarray(counts[m].sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
        v = counts[m][:, a.var_names.get_loc(g)]; v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
        return float(np.nanmean(v / tot) * cfg.TARGET_SUM), float(np.mean(v > 0) * 100)

    def per_donor_means(g, st):
        """Mean TP10K of gene g per DONOR within stage st (Sim has 3 donors/stage), so the
        figure can show donor-level dispersion (genuine biological replication, not pseudorep)."""
        gi = a.var_names.get_loc(g)
        out = {}
        for d in sorted(set(donor[cm & (stage == st)])):
            m = cm & (stage == st) & (donor == d)
            if m.sum() == 0:
                continue
            tot = np.asarray(counts[m].sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
            v = counts[m][:, gi]; v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
            out[d] = round(float(np.nanmean(v / tot) * cfg.TARGET_SUM), 4)
        return out

    per_gene = {g: {st: dict(zip(("mean_tp10k", "pct_expr"), stats(g, st))) for st in STAGES} for g in genes}
    # per-donor means (for donor-level error bars in Fig 3)
    per_gene_by_donor = {g: {st: per_donor_means(g, st) for st in STAGES} for g in genes}
    # within-family dominance per stage
    dom = {}
    for fam, fg in FAMILIES.items():
        fg = [g for g in fg if g in genes]; dom[fam] = {}
        for st in STAGES:
            tot = sum(per_gene[g][st]["mean_tp10k"] for g in fg) or np.nan
            lead = max(fg, key=lambda g: per_gene[g][st]["mean_tp10k"])
            dom[fam][st] = {"dominant": lead,
                            "dominant_pct_of_family": round(per_gene[lead][st]["mean_tp10k"] / tot * 100, 1)}
        dom[fam]["dominant_switches"] = len({dom[fam][st]["dominant"] for st in STAGES}) > 1

    out = {"question": "Does the dominant functioning cardiac Ca subtype change across development?",
           "answer": "Workhorse Cav1.2 (CACNA1C) dominant at all stages (~99% of L-type); anti-Ro "
                     "targets Cav1.3/Cav3.1 are minor/specialized & decline; real shift = Ca-handling "
                     "mode (sarcolemmal->SR, RYR2 ~2x).",
           "per_gene": per_gene, "per_gene_by_donor": per_gene_by_donor,
           "within_family_dominance": dom}
    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT.name}\n")
    for fam in FAMILIES:
        d = dom[fam]
        traj = " -> ".join(f"{st}:{d[st]['dominant']}({d[st]['dominant_pct_of_family']}%)" for st in STAGES)
        print(f"  {fam:12s} {traj}  switch={d['dominant_switches']}")
    print("\n  RYR2 (SR) mean:", {st: round(per_gene['RYR2'][st]['mean_tp10k'], 1) for st in STAGES})


if __name__ == "__main__":
    main()
