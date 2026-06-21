"""
config_chb.py — paths, thresholds, gene panels, and ion-channel classification for
the cardiac ion-channel developmental-expression analysis. Self-contained (this
public analysis directory has no external project dependency).

All values here are curated inputs (gene panels, label vocabularies, QC
thresholds) or standard references, not results. Sources for channel
classification: standard cardiac electrophysiology (Nerbonne & Kass 2005 Physiol
Rev; Catterall 2011; Amin et al. 2010).
"""
from __future__ import annotations

from datetime import date
from pathlib import Path

# --- paths ---
HERE = Path(__file__).resolve().parent           # public/scripts
PUBLIC = HERE.parent                              # public
DATA_DIR = PUBLIC.parent / "private" / "data"     # heavy raw/processed data (private)
PROC_DIR = DATA_DIR / "processed"
OUT_TABLES = PUBLIC / "results" / "tables"
# Reports are interpretive narrative (.md): in the full repo they are written PRIVATE
# (only computed tables and figures are public). When the public repo is checked out
# standalone (no sibling private/ dir), reports fall back to a local results/reports/
# so every script runs self-contained.
_PRIVATE = PUBLIC.parent / "private"
OUT_REPORTS = (_PRIVATE / "reports") if _PRIVATE.is_dir() else (PUBLIC / "results" / "reports")

# --- analysis parameters ---
SEED = 42
EPS = 1e-3
TARGET_SUM = 1e4                  # TP10K normalization
MIN_GENES_PER_CELL = 200
MAX_PCT_MT = 20.0
MIN_CELLS_PER_GROUP = 20

_VERSION = "CHB v1 (standalone)"


def header(input_desc: str) -> str:
    return f"# version: {_VERSION} | generated: {date.today().isoformat()} | seed: {SEED} | input: {input_desc}"


# ---------------------------------------------------------------------------
# PLATFORM CONSISTENCY for the adult conduction object (Kanemaru).
# The Kanemaru conduction atlas is multi-platform (10x snRNA + Multiome-RNA + scRNA).
# The dissected SAN/AVN regions were sequenced ONLY by Multiome-RNA single nuclei, so
# the adult pacemaker states (SAN_P_cell / AVN_P_cell) exist only in the Multiome data
# (0 snRNA cells). To keep within-adult node-vs-working comparisons on a single platform,
# every Kanemaru group (including ventricular working vCM, which otherwise pools ~75%
# snRNA) is restricted to Multiome-RNA. This helper applies that restriction; it is a
# no-op on objects without a `modality` column (i.e. the fetal Lim/Protze and Sim
# objects), so callers can apply it uniformly. CURATED METHODOLOGICAL CONSTRAINT.
KANEMARU_PLATFORM = "Multiome-RNA"


def kanemaru_platform_mask(a):
    """Boolean mask selecting the platform-consistent (Multiome-RNA) Kanemaru cells.
    Returns all-True for objects with no `modality` column (non-Kanemaru)."""
    import numpy as np
    if "modality" not in a.obs.columns:
        return np.ones(a.n_obs, dtype=bool)
    return (a.obs["modality"].astype(str).values == KANEMARU_PLATFORM)


# ---------------------------------------------------------------------------
# REGION PURITY for the adult ventricular working reference (Kanemaru).
# Kanemaru's `vCM*` cell_state clustering includes working cardiomyocytes physically
# dissected from the AV-NODE tissue block (region == "AVN"): ~29% of the Multiome vCM
# cells come from the AVN region, not from ventricle. Pooling these into the working-CM
# reference contaminates it with node-adjacent myocytes (inflates nodal markers TBX3/
# HCN1/HCN4, deflates RYR2), artificially narrowing the node-vs-working contrast. The
# clean ventricular reference is restricted to the true ventricular dissection regions.
# (The SAN region's working CMs are labelled aCM*, i.e. atrial, so they do not enter a
# vCM-prefixed reference; only the AVN region leaks into vCM.) CURATED CONSTRAINT.
KANEMARU_VENTRICULAR_REGIONS = ("LV", "RV", "SP", "AX")   # true ventricle: LV/RV/septum/apex


def kanemaru_ventricular_mask(a):
    """Boolean mask restricting to true ventricular regions (LV/RV/SP/AX), excluding the
    AVN-dissected cells that Kanemaru labels vCM. No-op (all-True) on objects without a
    `region` column (non-Kanemaru). Used for the BULK-VENTRICLE reference (node vs distal
    myocardium). For the node-matched working reference use kanemaru_node_working_mask."""
    import numpy as np
    if "region" not in a.obs.columns:
        return np.ones(a.n_obs, dtype=bool)
    return np.isin(a.obs["region"].astype(str).values, KANEMARU_VENTRICULAR_REGIONS)


# ---------------------------------------------------------------------------
# NODE-MATCHED adult working reference. The FETAL datasets are dissected SAN / AVN
# tissue, so the fetal `working_CM` is working myocardium captured WITHIN the node
# block (node-adjacent / transitional working CM), not bulk ventricle. To make the
# adult pacemaker-vs-working contrast symmetric with the fetal one, the adult working
# reference is the working CMs dissected from the SAME node region:
#   - AVN comparison: working CMs (vCM*/aCM*) with region == "AVN"
#   - SAN comparison: working CMs (aCM*) with region == "SAN"  (the SAN sits in atrial
#     tissue, so its working myocardium is atrial-CM, not ventricular)
# This isolates pacemaker identity within the node microenvironment at both stages and
# is computed within one dissection/batch. CURATED CONSTRAINT.
KANEMARU_NODE_REGIONS = ("SAN", "AVN")


def kanemaru_node_working_mask(a, node=None):
    """Boolean mask: adult working cardiomyocytes dissected from a node region. Working CM
    = cell_state starting 'vCM' or 'aCM' (a working, non-pacemaker CM cluster). `node` =
    'SAN' or 'AVN' for one region, or None to pool BOTH node regions (the node-dissection
    working reference that mirrors the fetal pooled working_CM). Requires Kanemaru obs."""
    import numpy as np
    cs = a.obs["cell_state"].astype(str)
    reg = a.obs["region"].astype(str).values
    is_work = (cs.str.startswith("vCM") | cs.str.startswith("aCM")).to_numpy()
    in_region = (reg == node) if node is not None else np.isin(reg, KANEMARU_NODE_REGIONS)
    return is_work & in_region


# ---------------------------------------------------------------------------
# Voltage-gated Ca-channel panel used for the conduction-system contrast.
# CACNA1C/D = L-type; CACNA1G/H/E = T/R-type. HCN4/SCN5A are pacemaker-orienting
# reference genes (expected nodal up / down). Curated input.
# ---------------------------------------------------------------------------
CHB_CA_PANEL = ("CACNA1C", "CACNA1D", "CACNA1G", "CACNA1H", "CACNA1E")
CHB_REF_GENES = ("HCN4", "SCN5A")
CHB_CA_BUDGET_GENES = ("CACNA1C", "CACNA1D", "CACNA1G", "CACNA1H", "CACNA1E",
                       "CACNA1A", "CACNA1B", "CACNA1I")


def chb_panel():
    return list(dict.fromkeys(CHB_CA_PANEL + CHB_REF_GENES))


# ---------------------------------------------------------------------------
# AUTOANTIBODY-TARGET -> CLINICAL PHENOTYPE -> LIFE-STAGE map. CURATED INPUT from
# the anti-Ro/autoimmune-channelopathy literature (see project notes). Used to test
# the developmental-vulnerability question: WHY is the fetal heart hit by AV block
# from maternal anti-Ro, while the SAME antibodies leave the maternal heart
# asymptomatic or cause a DIFFERENT (non-block) phenotype?
#
# Prediction tested: if differential channel EXPRESSION TIMING explains this, then
#   - FETAL-BLOCK arm (conduction Ca channels) should be FETAL / pacemaker-enriched
#   - ADULT-PHENOTYPE arm (hERG/IKs, Nav1.5) should be ADULT / working-CM-enriched
#
# fields: (channel_protein, clinical_phenotype, phenotype_stage, predicted_enriched, antibody_note)
#   phenotype_stage : fetal | adult | both
#   predicted_enriched : fetal_pacemaker | adult_working | ubiquitous
TARGET_PHENOTYPE = {
    # --- FETAL-BLOCK arm: anti-Ro Ca-channel targets, conduction system ---
    "CACNA1D": ("Cav1.3 L-type", "AV block (CHB)", "fetal", "fetal_pacemaker",
                "anti-Ro target; Cav1.3 E1 dom-I S5-S6 epitope"),
    "CACNA1G": ("Cav3.1 T-type", "AV block (CHB), irreversible", "fetal", "fetal_pacemaker",
                "anti-Ro target; mapped fetal epitope p305; headline target"),
    "CACNA1H": ("Cav3.2 T-type", "AV block (CHB)", "fetal", "fetal_pacemaker",
                "anti-Ro target (Lazzerini2017/Qu2019); Cav3.2 embryonic-predominant T-type"),
    "HCN4":    ("HCN4 funny current", "nodal automaticity", "fetal", "fetal_pacemaker",
                "pacemaker-orienting reference (not a primary antibody target)"),
    "CACNA1C": ("Cav1.2 L-type", "AV block (fetal) + adult arrhythmia", "both", "ubiquitous",
                "anti-Ro target but ubiquitous/diffuse — within-arm control"),
    # --- ADULT-PHENOTYPE arm: non-block adult autoimmune channelopathy ---
    "KCNH2":   ("hERG / Kv11.1 (IKr)", "autoimmune long-QT / TdP", "adult", "adult_working",
                "anti-Ro52 reduces IKr; S5-pore turret epitope; ADULT-only phenotype"),
    "KCNQ1":   ("Kv7.1 (IKs)", "autoimmune SHORT-QT / VT (agonist-like, increases IKs)", "adult", "adult_working",
                "CORRECTED per Lazzerini2017 Table 1: autoimmune anti-KCNQ1 is AGONIST-like "
                "(increases IKs -> QTc SHORTENING), not long-QT. The autoimmune LONG-QT is hERG/KCNH2."),
    "SCN5A":   ("Nav1.5", "conduction defects (anti-Nav1.5)", "adult", "adult_working",
                "SEPARATE antibody entity (NOT anti-Ro; Brugada/Korkmaz rat model) — "
                "channel context only; node is Na-poor at all ages"),
    "KCNA4":   ("Kv1.4 (Ito)", "autoimmune QTc prolongation / TdP", "adult", "adult_working",
                "anti-Kv1.4 in myasthenia gravis (Lazzerini2017 T1); adult QT target"),
    "KCNJ5":   ("Kir3.4 / GIRK4", "autoimmune atrial fibrillation", "adult", "adult_working",
                "anti-Kir3.4 -> AF (Lazzerini2025 Fig 1, tachy branch)"),
    "KCNJ11":  ("Kir6.2 / KATP", "autoimmune VF / cardiac arrest", "adult", "adult_working",
                "anti-Kir6.2 -> VF/CA (Lazzerini2025 Fig 1, tachy branch)"),
    # --- NEW fetal CHB biomarker (added 2026-06-14, Benjamin 2025): transporter, not a channel ---
    "ATP1A1":  ("Na+/K+-ATPase alpha-1 (TRANSPORTER, not a channel)", "AV block (CHB) — predictive biomarker",
                "fetal", "fetal_pacemaker",
                "anti-ATP1A1 maternal biomarker; 100% of CHB-affected, absent unaffected (Benjamin 2025). "
                "Test whether expression is fetal-enriched like the Ca-block arm."),
}
# the arms, by predicted vulnerable stage
FETAL_BLOCK_ARM = ("CACNA1D", "CACNA1G", "CACNA1H", "HCN4", "ATP1A1")
ADULT_PHENOTYPE_ARM = ("KCNH2", "KCNQ1", "SCN5A", "KCNA4", "KCNJ5", "KCNJ11")

# SHARED Fig-1 gene panel — ONE channel set used by BOTH the cell-type panels (1B,C,
# pacemaker-vs-working) and the stage panel (1D, fetal/adult ratio), so the figure is
# internally consistent. All are autoantibody-target channels (any phenotype), ordered
# fetal/brady -> control -> adult/tachy. (The pacemaker reference HCN4, the non-targeted
# negative control KCND3, and the identity marker GJA5/Cx40 were removed from Fig 1 — they are
# not antibody targets and the figure shows the targeted-channel program only.)
FIG1_PANEL_GENES = (
    "CACNA1G", "CACNA1D", "CACNA1H",                  # fetal / bradyarrhythmia (CHB)
    "CACNA1C",                                         # both (control)
    "KCNH2", "KCNQ1", "KCNA4", "KCNJ5", "KCNJ11", "SCN5A",   # adult / tachyarrhythmia
)
FIG1_IDENTITY_GENES = ()          # (GJA5 removed from Fig 1)
FIG1_RATIO_EXCLUDE = ()           # (HCN4 removed from Fig 1; nothing left to exclude)

# ---------------------------------------------------------------------------
# CANONICAL CHANNEL LABELS — single source of truth for figure text, so every
# panel labels a channel identically. Convention: protein/channel name first,
# gene symbol in parentheses (physiology/EP convention; gene kept for unambiguity,
# e.g. Cav3.1 = CACNA1G). hERG leads for KCNH2 (autoimmune-LQT literature term).
# Ca-handling / identity proteins use their unambiguous protein name (no gene).
# ---------------------------------------------------------------------------
CHANNEL_LABEL = {
    "CACNA1C": "Cav1.2 (CACNA1C)",
    "CACNA1D": "Cav1.3 (CACNA1D)",
    "CACNA1G": "Cav3.1 (CACNA1G)",
    "CACNA1H": "Cav3.2 (CACNA1H)",
    "CACNA1I": "Cav3.3 (CACNA1I)",
    "KCNH2":   "hERG (KCNH2)",
    "KCNQ1":   "Kv7.1 (KCNQ1)",
    "KCNA4":   "Kv1.4 (KCNA4)",
    "KCNA5":   "Kv1.5 (KCNA5)",
    "KCND3":   "Kv4.3 (KCND3)",
    "KCNJ2":   "Kir2.1 (KCNJ2)",
    "KCNJ3":   "Kir3.1 (KCNJ3)",
    "KCNJ5":   "Kir3.4 (KCNJ5)",
    "KCNJ11":  "Kir6.2 (KCNJ11)",
    "SCN5A":   "Nav1.5 (SCN5A)",
    "HCN4":    "HCN4",
    "HCN1":    "HCN1",
    "SLC8A1":  "NCX1",
    "ATP2A2":  "SERCA2a",
    "RYR2":    "RyR2",
    "PLN":     "PLN",
    "CASQ2":   "CASQ2",
    "ATP2B4":  "PMCA4",
    "CALM1":   "CaM1",
    "ATP1A1":  "Na/K-ATPase α1",
    "GJA5":    "Cx40 (GJA5)",
    "GJC1":    "Cx45 (GJC1)",
    "TBX3":    "TBX3",
    "SHOX2":   "SHOX2",
    "ISL1":    "ISL1",
}


def channel_label(gene, short=False):
    """Canonical display label for a gene. short=True returns the protein name only
    (before the parenthesis) for compact panels."""
    lab = CHANNEL_LABEL.get(gene, gene)
    if short and " (" in lab:
        return lab.split(" (")[0]
    return lab


# Ionic current carried by each voltage-gated channel — used only where the current
# context is informative (the Fig-4 detailed heatmap). Channels with no canonical single
# current (handling proteins, transporters, identity markers) are omitted.
CHANNEL_CURRENT = {
    "CACNA1C": "L", "CACNA1D": "L", "CACNA1G": "T", "CACNA1H": "T", "CACNA1I": "T",
    "KCNH2": "IKr", "KCNQ1": "IKs", "KCNA5": "IKur", "KCND3": "Ito",
    "KCNJ2": "IK1", "KCNJ3": "IKACh", "KCNJ5": "IKACh", "KCNJ11": "IKATP",
    "HCN4": "If", "HCN1": "If",
}


def channel_label_current(gene):
    """Fig-4 label: 'protein / current (gene)' where a current is defined, else the
    canonical label. e.g. CACNA1G -> 'Cav3.1 / T (CACNA1G)'; RYR2 -> 'RyR2'."""
    cur = CHANNEL_CURRENT.get(gene)
    if cur is None:
        return CHANNEL_LABEL.get(gene, gene)
    base = CHANNEL_LABEL.get(gene, gene)
    if " (" in base:
        protein, paren = base.split(" (", 1)
        return f"{protein} / {cur} ({paren}"
    return f"{base} / {cur}"

# ---------------------------------------------------------------------------
# FIGURE COLOURING — by the antibody's CLINICAL PHENOTYPE, not by the expression
# conclusion. Colouring channels by their expected arm would make
# the figure look pre-sorted toward the intended result. Instead the visual encoding
# is the EXTERNALLY-KNOWN disease phenotype of each antibody (taken from the clinical
# `TARGET_PHENOTYPE` fields), so the reader can see whether expression falls out along
# those clinical lines rather than being coloured to look as though it does.
#
#   COLOUR  = phenotype STAGE  (the life stage the antibody causes disease in:
#             fetal vs adult vs both/control) — from `phenotype_stage`.
#   HATCH   = phenotype TYPE   (the clinical syndrome: heart-block vs QT-alteration vs
#             other/conduction) — derived from the `phenotype` clinical string.
# Both are properties of the disease, defined in TARGET_PHENOTYPE with citations; the
# expression data are NOT used to assign them.
PHENOTYPE_STAGE_COLOR = {
    "fetal":   "#c44e52",   # red    — antibody causes FETAL disease (CHB)
    "adult":   "#4c72b0",   # blue   — antibody causes ADULT disease (LQT/arrhythmia)
    "both":    "#dd8452",   # orange — both stages (Cav1.2 control)
    "none":    "#8c8c8c",   # grey   — not an autoantibody target (control / non-targeted)
}
# Phenotype TYPE = the arrhythmic phenotype of the antibody, from Lazzerini & Boutjdir 2025
# Fig 1: BRADYarrhythmia (AVB / sinus brady / SA block) vs TACHYarrhythmia (LQTS/SQTS/BrS/
# AF/VF-CA). A channel acting in BOTH branches (Cav1.2: AV block in node cells + VF/QT in
# ventricular CMs; Nav1.5: SA/AV block + Brugada) is 'both' — the same channel gives brady
# or tachy depending on epitope location and cell type. 'none' = not an antibody target.
PHENOTYPE_TYPE_HATCH = {
    "brady": "",        # solid   — bradyarrhythmia (AV block / sinus brady / SA block)
    "tachy": "////",    # hatched — tachyarrhythmia (LQTS / SQTS / BrS / AF / VF-CA)
    "both":  "xx",      # cross   — both (cell-type / epitope-dependent, e.g. Cav1.2, Nav1.5)
    "none":  "",        # no hatch — a non-targeted gene has no antibody, hence no arrhythmic
                        # phenotype; its grey STAGE colour ("not targeted") is the marker, and
                        # it is not given an arrhythmia-type hatch (it is not in that legend).
}

# Non-targeted channel shown as the negative control (Fig 1 stage panel): a channel that is
# ROBUSTLY CARDIAC-EXPRESSED yet NOT in the cited autoimmune-channelopathy classification — the
# "expressed but not attacked" control that excludes the trivial "only-expressed channels get
# targeted" explanation. KCND3/Kv4.3 (Ito, the main transient-outward repolarizing K current)
# is well expressed in cardiomyocytes at both stages yet is absent from the target table.
# NOTE: Cav2.3/CACNA1E (the previous control) was REPLACED because it is at the snRNA detection
# floor in cardiac cells (~0.01-0.2 linear, vs targets 1-60), so it is effectively non-cardiac-
# expressed and cannot serve as an "expressed-yet-not-attacked" control. There is no cardiac-
# expressed non-targeted Cav (the cardiac voltage-gated Ca channels are essentially all anti-Ro
# targets), so the negative control is necessarily cross-family (a non-targeted K channel).
NON_TARGETED_CONTROL = ("KCND3",)


# --- phenotype encoding derived from the CITED classification table -------------------
# channelopathy_classification.json (built by build_channelopathy_classification.py from
# Lazzerini/Qu/Strandberg/Cao/Benjamin etc.) is the authoritative source. The
# figure encoding (stage colour, type hatch, anti-Ro/La bold) is derived from it rather than from a
# parallel hand-coded scheme, so there is one cited source of truth. Loaded lazily.
_CHANNELOPATHY = None


def _load_channelopathy():
    global _CHANNELOPATHY
    if _CHANNELOPATHY is None:
        import json
        p = OUT_TABLES / "channelopathy_classification.json"
        rows = json.loads(p.read_text())["classification"] if p.exists() else []
        by_gene = {}
        for r in rows:
            by_gene.setdefault(r["gene"], []).append(r)
        _CHANNELOPATHY = by_gene
    return _CHANNELOPATHY


def phenotype_stage(gene):
    """Clinical life-stage from the cited classification: fetal | adult | both | none.
    A gene with both a fetal and an adult classification row (e.g. CACNA1C) -> 'both'."""
    rows = _load_channelopathy().get(gene, [])
    if not rows:
        return "none"
    stages = {("fetal" if "fetal" in r.get("life_stage", "").lower() else "adult") for r in rows}
    if stages == {"fetal", "adult"}:
        return "both"
    return "fetal" if "fetal" in stages else "adult"


def phenotype_type(gene):
    """Arrhythmic phenotype TYPE from the cited classification's `arrhythmia_class` field
    (Lazzerini 2025 Fig 1): brady | tachy | both | none. A gene with both a brady and a tachy
    classification row -> 'both' (same channel, different cell type / epitope)."""
    rows = _load_channelopathy().get(gene, [])
    classes = {r.get("arrhythmia_class") for r in rows if r.get("arrhythmia_class")}
    if not classes:
        return "none"
    if classes == {"brady", "tachy"}:
        return "both"
    return classes.pop()


# CONGENITAL-HEART-BLOCK-related antibody targets — shown by BOLDING the gene-name label
# (a typographic mark, NOT a colour/shape data channel). A gene is CHB-related if the cited
# classification gives it an "autoimmune congenital heart block" row (any antibody — anti-Ro,
# anti-La, or the anti-ATP1A1 biomarker), i.e. it is a target of a CHB-causing antibody.
def is_chb_target(gene):
    return any("congenital heart block" in r.get("autoimmune_class", "").lower()
               for r in _load_channelopathy().get(gene, []))


def phenotype_stage_color(gene):
    return PHENOTYPE_STAGE_COLOR.get(phenotype_stage(gene), "#8c8c8c")


def phenotype_type_hatch(gene):
    return PHENOTYPE_TYPE_HATCH.get(phenotype_type(gene), "..")


def target_phenotype(gene):
    return TARGET_PHENOTYPE.get(gene)


# ion family / role classification (subset needed by the panel scripts)
ION_CLASS = {
    "CACNA1C": ("Ca", "Cav", "pore"), "CACNA1D": ("Ca", "Cav", "pore"),
    "CACNA1G": ("Ca", "Cav", "pore"), "CACNA1H": ("Ca", "Cav", "pore"),
    "CACNA1I": ("Ca", "Cav", "pore"), "CACNA1E": ("Ca", "Cav", "pore"),
    "CACNA1A": ("Ca", "Cav", "pore"), "CACNA1B": ("Ca", "Cav", "pore"),
    "HCN4": ("cation", "HCN", "pore"),
    "SCN5A": ("Na", "Nav", "pore"),
    "KCND3": ("K", "Kv", "pore"),   # Kv4.3 (Ito) — non-targeted negative control
}


def ion_type(g):   return ION_CLASS.get(g, ("none", "other", "other"))[0]
def ion_family(g): return ION_CLASS.get(g, ("none", "other", "other"))[1]
def ion_role(g):   return ION_CLASS.get(g, ("none", "other", "other"))[2]
