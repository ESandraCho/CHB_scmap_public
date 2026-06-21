#!/usr/bin/env python3
"""
build_channelopathy_classification.py — classification of cardiac ion-channel (and
related) targets by AUTOIMMUNE channelopathy type (congenital/fetal vs adult) and
clinical SYMPTOM (AV block, QT prolongation, etc.), curated from the cited literature.

This is a CURATED LITERATURE TABLE (not computed from data) — the evidence backbone for
the fetal-vs-adult vulnerability analysis. Each row: target gene/protein, channel/current,
autoimmune antibody, channelopathy class, life-stage, clinical symptom/ECG, EP mechanism,
and the supporting reference(s). A separate column flags the corresponding INHERITED
(genetic) channelopathy for the same channel (from Kline 2019), so the autoimmune and
genetic phenotypes can be compared.

Sources (cited literature):
  Lazzerini 2017 Nat Rev Cardiol (nrcardio.2017.61) — master Table 1 (Ca/K/Na)
  Qu 2019 Front Cardiovasc Med (fcvm-06-00054) — Ca channelopathies, Table 1
  Jin Li 2020 Curr Cardiol Rep (11886_2020) — atrial/nodal/ventricular Table 1
  Capecchi 2019 Heart Rhythm (S1547527119301390) — K+ channelopathies
  Cao/Lu 2025 Front Immunol (fimmu-16-1561061) — anti-Ro adult arrhythmias
  Lazzerini & Boutjdir 2025 Heart Rhythm (PIIS1547527125021010) — contemporary review
  Lazzerini 2023 JACC EP (S2405500X23002050) + Keefe 2023 editorial — anti-Ro Ca-block AVB
  Benjamin 2025 Lancet Rheumatol (S266599132500092X) — anti-ATP1A1 fetal CHB biomarker
  Kline 2019 Med Clin N Am (S0025712519300495) — INHERITED channelopathy genotype-phenotype

Env: any python3 (stdlib only — pure curation -> JSON + TSV).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_TABLES = HERE.parent / "results" / "tables"
OUT_JSON = OUT_TABLES / "channelopathy_classification.json"
OUT_TSV = OUT_TABLES / "channelopathy_classification.tsv"

# Each entry: curated from the cited reference Table(s). `arrhythmia_class` (brady | tachy)
# follows the Lazzerini & Boutjdir 2025 Fig-1 "classification by arrhythmic phenotype":
# Bradyarrhythmias = AVB / sinus bradycardia / SA block; Tachyarrhythmias = LQTS / SQTS /
# BrS / AF / VF-CA. A channel that appears in BOTH branches (Cav1.2, Nav1.5) gets one row
# per branch, so its overall type resolves to "both".
# fields: gene, channel_current, antibody, autoimmune_class, arrhythmia_class, life_stage,
#         symptom_ecg, ep_mechanism, inherited_counterpart, refs
CLASSIFICATION = [
    # ===== CONGENITAL / FETAL autoimmune (anti-Ro/SSA + anti-La/SSB) — BRADY (AVB/SB) =====
    {"gene": "CACNA1D", "channel_current": "Cav1.3 L-type ICaL",
     "antibody": "anti-Ro/SSA (52kD)", "autoimmune_class": "autoimmune congenital heart block",
     "arrhythmia_class": "brady",
     "life_stage": "fetal/neonatal", "symptom_ecg": "AV block, sinus bradycardia",
     "ep_mechanism": "inhibits ICaL -> nodal AP inhibition (SAN/AVN pacemaking)",
     "inherited_counterpart": "SANDD (sinoatrial node dysfunction + deafness)",
     "refs": ["Lazzerini2017_T1", "Qu2019_T1", "JinLi2020_T1", "Lazzerini2025"]},
    {"gene": "CACNA1C", "channel_current": "Cav1.2 L-type ICaL",
     "antibody": "anti-Ro/SSA (52kD)", "autoimmune_class": "autoimmune congenital heart block",
     "arrhythmia_class": "brady",
     "life_stage": "fetal/neonatal (+ adult acquired)", "symptom_ecg": "AV block, sinus bradycardia",
     "ep_mechanism": "inhibits ICaL at AVN node cells (S5-S6 pore loop domain I+III); ubiquitous channel",
     "inherited_counterpart": "Timothy syndrome / LQT8; Brugada 3",
     "refs": ["Lazzerini2017_T1", "Qu2019_T1", "Lazzerini2023_JACCEP", "Keefe2023", "Lazzerini2025"]},
    {"gene": "CACNA1G", "channel_current": "Cav3.1 T-type ICaT",
     "antibody": "anti-Ro/SSA", "autoimmune_class": "autoimmune congenital heart block",
     "arrhythmia_class": "brady",
     "life_stage": "fetal/neonatal", "symptom_ecg": "sinus bradycardia, AV block",
     "ep_mechanism": "inhibits ICaT -> nodal AP inhibition; mapped fetal epitope p305 (S5-S6 repeat I)",
     "inherited_counterpart": "(not a classic inherited channelopathy gene)",
     "refs": ["Lazzerini2017_T1", "Qu2019_T1", "Strandberg2013", "Lazzerini2025"]},
    {"gene": "CACNA1H", "channel_current": "Cav3.2 T-type ICaT",
     "antibody": "anti-Ro/SSA", "autoimmune_class": "autoimmune congenital heart block",
     "arrhythmia_class": "brady",
     "life_stage": "fetal (embryonic-predominant)", "symptom_ecg": "sinus bradycardia, AV block",
     "ep_mechanism": "inhibits ICaT (weaker, ~19%); predominant embryonic T-type isoform",
     "inherited_counterpart": "(epilepsy/other; minor cardiac)",
     "refs": ["Lazzerini2017_T1", "Qu2019_T1", "Lazzerini2025"]},
    {"gene": "ATP1A1", "channel_current": "Na+/K+-ATPase alpha-1 (transporter, NOT a channel)",
     "antibody": "anti-ATP1A1 (maternal)", "autoimmune_class": "autoimmune congenital heart block (biomarker)",
     "arrhythmia_class": "brady",
     "life_stage": "fetal/neonatal", "symptom_ecg": "AV block (predictive biomarker)",
     "ep_mechanism": "novel fetal cardiac protein target; 100% of CHB-affected, absent unaffected",
     "inherited_counterpart": "(ATP1A1 mutations -> neuro/renal, not classic CHB)",
     "refs": ["Benjamin2025"]},

    # ===== ADULT autoimmune — TACHY (LQTS / SQTS / AF / VF-CA) =====
    {"gene": "KCNH2", "channel_current": "Kv11.1 / hERG IKr",
     "antibody": "anti-Ro/SSA (52kD)", "autoimmune_class": "autoimmune long-QT syndrome",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "QTc prolongation, complex ventricular arrhythmias, TdP",
     "ep_mechanism": "inhibits IKr (acute pore block +/or chronic internalization) -> APD prolongation; S5-pore turret epitope",
     "inherited_counterpart": "LQT2 (LoF) / SQT1 (GoF)",
     "refs": ["Lazzerini2017_T1", "Capecchi2019", "Cao2025", "JinLi2020_T1", "Lazzerini2025"]},
    {"gene": "KCNQ1", "channel_current": "Kv7.1 IKs",
     "antibody": "agonist-like anti-KCNQ1", "autoimmune_class": "autoimmune short-QT / DCM-associated",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "QTc shortening, ventricular tachycardia",
     "ep_mechanism": "increases IKs (agonist-like) -> APD shortening",
     "inherited_counterpart": "LQT1 (LoF) / SQT2 (GoF) / familial AF",
     "refs": ["Lazzerini2017_T1", "Capecchi2019", "JinLi2020_T1", "Lazzerini2025"]},
    {"gene": "KCNA4", "channel_current": "Kv1.4 Ito",
     "antibody": "anti-Kv1.4 (myasthenia gravis)", "autoimmune_class": "autoimmune QT prolongation",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "QTc prolongation, TdP",
     "ep_mechanism": "inhibiting (mechanism NA in table)",
     "inherited_counterpart": "(not a classic isolated cardiac channelopathy)",
     "refs": ["Lazzerini2017_T1", "Lazzerini2025"]},
    {"gene": "KCNJ5", "channel_current": "Kir3.4 / GIRK4 IK(ACh)",
     "antibody": "anti-Kir3.4", "autoimmune_class": "autoimmune atrial fibrillation",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "atrial fibrillation",
     "ep_mechanism": "anti-Kir3.4 autoantibody associated with AF (Lazzerini 2025 Fig 1)",
     "inherited_counterpart": "(KCNJ5 -> primary aldosteronism; LQT13)",
     "refs": ["Lazzerini2025"]},
    {"gene": "KCNJ11", "channel_current": "Kir6.2 / KATP",
     "antibody": "anti-Kir6.2", "autoimmune_class": "autoimmune VF / cardiac arrest (short-coupled VF)",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "ventricular fibrillation / cardiac arrest",
     "ep_mechanism": "anti-Kir6.2 (KATP) autoantibody associated with VF/CA (Lazzerini 2025 Fig 1)",
     "inherited_counterpart": "(KCNJ11 -> neonatal diabetes; Cantu)",
     "refs": ["Lazzerini2025"]},

    # ===== channels appearing in BOTH brady and tachy branches (-> 'both') =====
    {"gene": "SCN5A", "channel_current": "Nav1.5 INa",
     "antibody": "anti-Nav1.5 (SEPARATE entity, NOT anti-Ro)", "autoimmune_class": "autoimmune AV block / SA block (adult)",
     "arrhythmia_class": "brady",
     "life_stage": "adult", "symptom_ecg": "AV block / sinoatrial block (idiopathic, adult)",
     "ep_mechanism": "reduces INa (downregulation); AVN vulnerable due to low Na density",
     "inherited_counterpart": "LQT3 / Brugada / PCCD / SCN5A-DCM",
     "refs": ["Lazzerini2017_T1", "JinLi2020_T1", "Lazzerini2025"]},
    {"gene": "SCN5A", "channel_current": "Nav1.5 INa",
     "antibody": "anti-Nav1.5 (SEPARATE entity, NOT anti-Ro)", "autoimmune_class": "autoimmune Brugada-pattern / VA",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "Brugada-type ECG, ventricular arrhythmia",
     "ep_mechanism": "anti-Nav1.5 associated with Brugada-syndrome phenotype (Lazzerini 2025 Fig 1, tachy branch)",
     "inherited_counterpart": "LQT3 / Brugada / PCCD / SCN5A-DCM",
     "refs": ["Lazzerini2025"]},
    {"gene": "CACNA1C", "channel_current": "Cav1.2 L-type ICaL (agonist-like)",
     "antibody": "agonist-like anti-Cav1.2 (anti-alpha1C)", "autoimmune_class": "DCM-associated ventricular arrhythmia / VF-CA",
     "arrhythmia_class": "tachy",
     "life_stage": "adult", "symptom_ecg": "ventricular tachycardia, VF/cardiac arrest, sudden cardiac death",
     "ep_mechanism": "increases ICaL in ventricular CMs (agonist-like, N-terminus D1) -> APD prolongation + EADs",
     "inherited_counterpart": "Timothy/LQT8",
     "refs": ["Lazzerini2017_T1", "Qu2019_T1", "Lazzerini2025"]},
]

REFERENCES = {
    "Lazzerini2017_T1": "Lazzerini et al. 2017 Nat Rev Cardiol 14:521 (Table 1) — nrcardio.2017.61.pdf",
    "Qu2019_T1": "Qu et al. 2019 Front Cardiovasc Med 6:54 (Table 1) — fcvm-06-00054.pdf",
    "JinLi2020_T1": "Jin Li 2020 Curr Cardiol Rep 23:3 (Table 1) — 11886_2020_Article_1430.pdf",
    "Capecchi2019": "Capecchi et al. 2019 Heart Rhythm (K+ channelopathies) — 1-s2.0-S1547527119301390-main.pdf",
    "Cao2025": "Cao, Lu et al. 2025 Front Immunol 16:1561061 — fimmu-16-1561061.pdf",
    "Lazzerini2025": "Lazzerini & Boutjdir 2025 Heart Rhythm (review) — PIIS1547527125021010.pdf",
    "Lazzerini2023_JACCEP": "Lazzerini et al. 2023 JACC Clin EP 9 — 1-s2.0-S2405500X23002050-main.pdf",
    "Keefe2023": "Keefe & Wehrens 2023 JACC Clin EP editorial — keefe-wehrens-2023-...pdf",
    "Strandberg2013": "Strandberg et al. 2013 PLoS One (Cav3.1 p305 fetal epitope)",
    "Benjamin2025": "Benjamin et al. 2025 Lancet Rheumatol 7:e554 (anti-ATP1A1) — 1-s2.0-S266599132500092X-main.pdf",
    "Kline2019": "Kline & Costantini 2019 Med Clin N Am 103:809 (inherited) — 1-s2.0-S0025712519300495.pdf",
}


def main():
    OUT_TABLES.mkdir(parents=True, exist_ok=True)
    payload = {
        "description": "Cardiac ion-channel targets classified by autoimmune channelopathy "
                       "type (congenital/fetal vs adult) and clinical symptom. Curated from "
                       "the cited literature. Inherited (genetic) counterpart listed for comparison.",
        "references": REFERENCES,
        "classification": CLASSIFICATION,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    cols = ["gene", "channel_current", "antibody", "autoimmune_class", "arrhythmia_class",
            "life_stage", "symptom_ecg", "ep_mechanism", "inherited_counterpart", "refs"]
    lines = ["\t".join(cols)]
    for r in CLASSIFICATION:
        lines.append("\t".join(str(r[c]) if c != "refs" else ";".join(r["refs"]) for c in cols))
    OUT_TSV.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT_JSON.name} + {OUT_TSV.name}\n")

    # console summary grouped by class+stage
    print("=== CONGENITAL / FETAL autoimmune channelopathies ===")
    for r in CLASSIFICATION:
        if "congenital" in r["autoimmune_class"]:
            print(f"  {r['gene']:8s} {r['channel_current']:28s} -> {r['symptom_ecg']}")
    print("\n=== ADULT autoimmune channelopathies ===")
    for r in CLASSIFICATION:
        if "congenital" not in r["autoimmune_class"]:
            print(f"  {r['gene']:8s} {r['channel_current']:32s} -> {r['symptom_ecg']}")


if __name__ == "__main__":
    main()
