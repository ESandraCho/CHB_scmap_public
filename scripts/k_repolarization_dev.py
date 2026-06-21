#!/usr/bin/env python3
"""
k_repolarization_dev.py — does the working myocardium become more IKr-DEPENDENT for
repolarization across development? Tests the substrate premise behind the fetal-vs-adult
QT asymmetry: the adult QT phenotype is an IKr/hERG (KCNH2) disease, yet the fetus is
QT-spared (Brucato 2004; magnetocardiography). If the fetal ventricle does not lean on
IKr for repolarization, an anti-hERG antibody has little functional consequence there.

Two transcript proxies, reported jointly (neither alone proves functional dependence):
  1. IKr SHARE of the repolarizing-K pool = KCNH2 / sum(repolarizing-K). A rising share
     with maturation is a transcript signature of leaning on IKr for repolarization.
  2. Repolarization RESERVE = N_eff (effective number of distinct K currents) =
     exp(Shannon entropy) of the per-CURRENT K composition. High reserve = many redundant
     repolarizing currents buffering a single-channel (IKr) block; low reserve = vulnerable.
  Joint readout: IKr-DEPENDENCE = high IKr share WITH low reserve.

HONEST SCOPE: transcript composition is NOT a measurement of functional IKr dependence
(that is an electrophysiology property); this is a proxy. Computed on Sim 2021 CM,
within-study, fetal/child/adult (the same-platform, confound-free stage axis), whole-CM
(working-dominated — appropriate for QT, a distributed-ventricle readout). The Sim atlas
has no dissected node, so this is the ventricular/working program, not the node.

Sim 2021 CM, within-study, mean-per-cell TP10K. Env: structural_epitope.
"""
from __future__ import annotations
import json
import numpy as np, scipy.sparse as sp
import config_chb as cfg

OUT = cfg.OUT_TABLES / "k_repolarization_dev.json"
SIM = cfg.PROC_DIR / "sim2021_qc.h5ad"
STAGES = ["fetal", "child", "adult"]

# repolarizing-K channels grouped by the distinct cardiac CURRENT they carry (so the reserve
# N_eff counts independent currents, not redundant paralogs). Curated (config_channel_subtypes
# / Nerbonne & Kass 2005). IKr (KCNH2) is the QT/hERG target; the others are the reserve.
K_CURRENTS = {
    "IKr":  ["KCNH2"],                       # rapid delayed rectifier — hERG; the QT/adult target
    "IKs":  ["KCNQ1"],                       # slow delayed rectifier
    "IK1":  ["KCNJ2", "KCNJ4", "KCNJ12"],    # inward rectifier (resting/terminal repol)
    "Ito":  ["KCND3", "KCNA4"],              # transient outward
    "IKur": ["KCNA5"],                       # ultrarapid delayed rectifier
    "IKATP":["KCNJ11", "KCNJ8"],             # K-ATP
    "IKACh":["KCNJ5", "KCNJ3"],              # G-protein-gated
}
ALL_K = [g for genes in K_CURRENTS.values() for g in genes]

# EAD (early-afterdepolarization) SUBSTRATE: when IKr/hERG is lost, the plateau prolongs into the
# window where reactivatable INWARD currents reopen and drive EADs -> torsades. The EAD-driving
# inward currents are the bulk L-type ICaL (Cav1.2/CACNA1C; Cav1.3/CACNA1D is a minor low-voltage
# L-type, reported but ~100x smaller) and the late/window Na current (Nav1.5/SCN5A). The reported
# quantity is the transcriptional BALANCE of these reactivatable inward channels relative to the outward
# repolarizing-K reserve. IMPORTANT framing: this ratio is a COMPOSITION of channel transcripts
# (a proxy for which channels the cell is PROVISIONED for / prefers), NOT a current/flux magnitude.
# A higher fetal ratio means the fetal cell carries a larger PROPORTION of L-type-Ca + Na channels
# that COULD be potentiated under a hERG block — it does NOT mean ~2x the current, nor that EADs
# occur (that is downstream functional electrophysiology not measured here).
EAD_INWARD_CAV12 = ["CACNA1C"]            # bulk ICaL — the main EAD driver (NOT a fetal-enriched target)
EAD_INWARD_CAV13 = ["CACNA1D"]            # Cav1.3 L-type (low-voltage; minor; a fetal-enriched anti-Ro target)
EAD_INWARD_LATENA = ["SCN5A"]             # late/window Na


def main():
    import anndata as ad
    a = ad.read_h5ad(SIM)
    cm = a.obs["cat"].astype(str).values == "CM"
    counts = a.layers["counts"]; stage = a.obs["stage_bin"].astype(str).values
    genes = [g for g in ALL_K if g in a.var_names]

    donor = a.obs["donor_id"].astype(str).values

    def mean_tp10k(g, st, don=None):
        m = cm & (stage == st)
        if don is not None:
            m = m & (donor == don)
        if m.sum() == 0:
            return np.nan
        tot = np.asarray(counts[m].sum(1)).ravel().astype(float); tot[tot == 0] = np.nan
        v = counts[m][:, a.var_names.get_loc(g)]
        v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
        return float(np.nanmean(v / tot) * cfg.TARGET_SUM)

    per_gene = {g: {st: mean_tp10k(g, st) for st in STAGES} for g in genes}

    # EAD inward-current genes (loaded separately; not part of the K pool)
    ead_genes = [g for g in (EAD_INWARD_CAV12 + EAD_INWARD_CAV13 + EAD_INWARD_LATENA)
                 if g in a.var_names]
    ead_tp10k = {g: {st: mean_tp10k(g, st) for st in STAGES} for g in ead_genes}

    # per-current totals, IKr share, and reserve N_eff per stage
    per_stage = {}
    for st in STAGES:
        cur = {c: sum(per_gene[g][st] for g in gs if g in per_gene)
               for c, gs in K_CURRENTS.items()}
        pool = sum(cur.values()) or np.nan
        ikr_share = cur["IKr"] / pool if pool else np.nan
        # reserve = effective number of distinct K currents (exp Shannon entropy over currents)
        ps = np.array([cur[c] for c in K_CURRENTS]) / pool
        ps = ps[ps > 0]
        neff = float(np.exp(-np.sum(ps * np.log(ps)))) if len(ps) else np.nan
        # reserve EXCLUDING IKr: how much non-IKr repolarization backs up a hERG block
        nonikr = {c: cur[c] for c in K_CURRENTS if c != "IKr"}
        npool = sum(nonikr.values()) or np.nan
        # EAD substrate: reactivatable INWARD (ICaL Cav1.2 [+Cav1.3] + late-Na) over OUTWARD
        # repolarizing reserve (full K pool). Higher => more EAD/torsades-prone under a hERG block.
        def g(gene): return ead_tp10k.get(gene, {}).get(st, 0.0)
        cav12 = sum(g(x) for x in EAD_INWARD_CAV12)
        cav13 = sum(g(x) for x in EAD_INWARD_CAV13)
        latena = sum(g(x) for x in EAD_INWARD_LATENA)
        inward_cav12 = cav12 + latena              # conservative: bulk ICaL only (NOT target-driven)
        inward_full = cav12 + cav13 + latena        # full L-type (incl. minor fetal-enriched Cav1.3)
        per_stage[st] = {
            "current_tp10k": {c: round(cur[c], 4) for c in K_CURRENTS},
            "K_pool_tp10k": round(pool, 4),
            "IKr_share": round(ikr_share, 4),
            "reserve_Neff": round(neff, 4),
            "nonIKr_pool_tp10k": round(npool, 4),
            "IKr_over_nonIKr": round(cur["IKr"] / npool, 4) if npool else None,
            "EAD_substrate": {
                "Cav12_ICaL_tp10k": round(cav12, 4),
                "Cav13_tp10k": round(cav13, 4),
                "lateNa_SCN5A_tp10k": round(latena, 4),
                "Krepol_outward_tp10k": round(pool, 4),
                "inward_over_outward_Cav12only": round(inward_cav12 / pool, 4) if pool else None,
                "inward_over_outward_Cav12plus13": round(inward_full / pool, 4) if pool else None,
            },
        }

    # PER-DONOR inward:outward ratio (for donor-level error bars in Fig 3C). Sim has 3 donors
    # per stage — genuine biological replication. Same definition as inward_over_outward_Cav12only
    # ((Cav1.2 ICaL + late-Na) / full K pool), computed within each donor.
    inward_over_outward_by_donor = {}
    for st in STAGES:
        dons = sorted(set(donor[cm & (stage == st)]))
        per_d = {}
        for d in dons:
            kpool = sum(mean_tp10k(g, st, d) for gs in K_CURRENTS.values() for g in gs
                        if g in per_gene)
            inw = sum(mean_tp10k(g, st, d) for g in (EAD_INWARD_CAV12 + EAD_INWARD_LATENA)
                      if g in ead_tp10k)
            if kpool and not np.isnan(kpool):
                per_d[d] = round(inw / kpool, 4)
        inward_over_outward_by_donor[st] = per_d

    out = {
        "question": "(1) Does the working myocardium become more IKr-dependent across development? "
                    "(2) Is the fetal CM a stronger EAD/torsades substrate under a hERG block?",
        "proxies": "IKr share of repolarizing-K pool + reserve N_eff (IKr dependence/reserve); "
                   "EAD substrate = reactivatable inward (ICaL Cav1.2[+Cav1.3] + late-Na SCN5A) over "
                   "outward repolarizing-K reserve. Transcript proxies, NOT measurements of current "
                   "density or of EADs themselves.",
        "scope": "Sim 2021 within-study, whole-CM (working-dominated), fetal/child/adult; no node.",
        "stages": STAGES,
        "per_current_genes": K_CURRENTS,
        "ead_inward_genes": {"Cav1.2_ICaL": EAD_INWARD_CAV12, "Cav1.3": EAD_INWARD_CAV13,
                             "lateNa": EAD_INWARD_LATENA},
        "per_gene_tp10k": per_gene,
        "ead_inward_tp10k": ead_tp10k,
        "per_stage": per_stage,
        "inward_over_outward_by_donor": inward_over_outward_by_donor,
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT.name}\n")
    print(f"{'stage':8s} {'IKr share':>10s} {'reserve Neff':>13s} {'IKr/nonIKr':>11s}  KCNH2 tp10k")
    for st in STAGES:
        s = per_stage[st]
        print(f"{st:8s} {s['IKr_share']:>10.3f} {s['reserve_Neff']:>13.3f} "
              f"{(s['IKr_over_nonIKr'] or float('nan')):>11.3f}  {per_gene.get('KCNH2', {}).get(st, float('nan')):.4f}")
    print(f"\n{'stage':8s} {'Cav1.2':>8s} {'SCN5A':>7s} {'K-repol':>8s} "
          f"{'in:out(1.2)':>12s} {'in:out(1.2+1.3)':>16s}")
    for st in STAGES:
        e = per_stage[st]["EAD_substrate"]
        print(f"{st:8s} {e['Cav12_ICaL_tp10k']:>8.2f} {e['lateNa_SCN5A_tp10k']:>7.3f} "
              f"{e['Krepol_outward_tp10k']:>8.3f} {e['inward_over_outward_Cav12only']:>12.2f} "
              f"{e['inward_over_outward_Cav12plus13']:>16.2f}")
    f0 = per_stage["fetal"]["EAD_substrate"]["inward_over_outward_Cav12only"]
    a0 = per_stage["adult"]["EAD_substrate"]["inward_over_outward_Cav12only"]
    print(f"\n  Fetal inward:outward balance is ~{f0 / a0:.1f}x the adult's (Cav1.2-only). NOTE: this is a "
          f"COMPOSITION/BALANCE of channel transcripts, NOT a current/flux measurement. It means the "
          f"fetal cell is provisioned with a larger PROPORTION of L-type Ca + Na channels that COULD be "
          f"reactivated/potentiated under a hERG block, vs its repolarizing-K reserve. Substrate "
          f"provision, not demonstrated EADs and not a doubling of current.")


if __name__ == "__main__":
    main()
