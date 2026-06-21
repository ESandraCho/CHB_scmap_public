"""
config_channel_classes.py — coarse FUNCTIONAL-CLASS grouping of cardiac ion
channels / Ca-handling machinery, for the within-cell channel-class COMPOSITION
analysis (how the ion-channel budget shifts across developmental stage, read
separately within each cell type).

CURATED INPUT, not a result. The class assignments reuse the two existing curated
configs and add nothing new:
  - calcium localisation (surface vs SR) from config_ca_accessibility.CA_MACHINERY
  - K / Na membership from config_channel_subtypes (ion field)
  - neuronal / skeletal isoforms (config_channel_subtypes.NON_CARDIAC_SUBTYPE_GENES)
    are EXCLUDED — they sit at the snRNA detection floor and only add variance to a
    composition fraction.

Four classes (the comparison axis):
  surface_Ca   antibody-accessible sarcolemmal Ca entry (Cav, NCX, ORAI, TRP)
  SR_Ca        intracellular SR Ca store machinery (RyR, SERCA2, TRIC-A) — the
               cardiac CICR store specifically (ER-signalling IP3R/STIM and the
               skeletal SERCA1 are NOT part of the contractile SR store and are
               excluded to keep the class a clean "SR Ca store")
  K            potassium channels (repolarisation: Kv, Kir)
  Na           cardiac sodium channels (Nav1.5 + auxiliary)

The composition is reported as WITHIN-GROUP FRACTIONS of the combined four-class
pool (summing to 1.0 per group); as an internal ratio it cancels a uniform per-
cell / per-dataset scaling and is robust to library-size offset. NOTE: a *cross-
stage* (fetal-vs-adult) composition is still cross-dataset (fetal Lim vs adult
Kanemaru) and subject to abundance-dependent capture bias; it is therefore
direction-robust, magnitude-cautious (see Discussion / Methods).
"""
from __future__ import annotations

import config_ca_accessibility as acc
import config_channel_subtypes as sub

# --- SR Ca-store members (curated subset of the INTRACELLULAR Ca machinery) ---
# The contractile SR store = RyR release + SERCA2 reuptake + TRIC-A SR counter-ion,
# plus the SERCA2 regulator phospholamban (PLN) and the SR luminal Ca buffer
# calsequestrin-2 (CASQ2) — both physical members of the cardiac SR Ca-store complex
# implicated in AV-block Ca handling (CIRCEP.115.003432). ER-signalling IP3R (ITPR1/2/3)
# and the SOCE ER sensor (STIM2) act in ER Ca signalling, not the contractile SR store;
# skeletal SERCA1 (ATP2A1) and ER TRIC-B (TMEM38B) are non-cardiac-SR — all excluded so
# SR_Ca is a clean store class.
_SR_CA = ("RYR2", "RYR3", "ATP2A2", "TMEM38A", "PLN", "CASQ2")

# --- HCN funny-current channels (pacemaker cation conductance) ---
# HCN4 (dominant nodal If) + HCN1 (fetal/nodal). Central to nodal automaticity and AV
# conduction; a distinct conductance class, not Ca/K/Na. Neuronal-minor HCN3 omitted.
_HCN = ("HCN4", "HCN1", "HCN2")

# --- Cation transporters / pumps that are not voltage-gated channels ---
# Na/K-ATPase alpha-1 (ATP1A1) — the maternal anti-ATP1A1 CHB biomarker (Benjamin 2025);
# included as an AVB-relevant surface transporter. (PMCA4/ATP2B4, a sarcolemmal Ca
# extrusion pump, is classed with surface_Ca via config_ca_accessibility.)
_TRANSPORTER = ("ATP1A1",)


# sarcolemmal Ca extrusion pump PMCA4 (ATP2B4) is a surface Ca-handling member but is
# not in the config_ca_accessibility Cav/NCX/store list; added here for the budget.
_EXTRA_SURFACE_CA = ("ATP2B4",)


def _surface_ca():
    base = tuple(g for g in acc.genes_by_access("SURFACE")
                 if g not in sub.NON_CARDIAC_SUBTYPE_GENES)
    return base + tuple(g for g in _EXTRA_SURFACE_CA if g not in base)


def _ion_class(ion):
    return tuple(g for g in sub.genes_for_ion(ion)
                 if g not in sub.NON_CARDIAC_SUBTYPE_GENES
                 and g not in _HCN)        # HCN carries cation, but is its own class here


# class -> ordered tuple of genes. Built from the curated configs at import time.
CHANNEL_CLASSES = {
    "surface_Ca":  _surface_ca(),
    "SR_Ca":       _SR_CA,
    "K":           _ion_class("K"),
    "Na":          _ion_class("Na"),
    "HCN":         _HCN,
    "transporter": _TRANSPORTER,
}

# display order + colours (consistent with the Fig 2 accessibility palette:
# surface = warm/accessible, SR = muted/buried, K/Na/HCN/transporter = neutral families)
CLASS_ORDER = ("surface_Ca", "SR_Ca", "K", "Na", "HCN", "transporter")
CLASS_LABELS = {
    "surface_Ca": "surface Ca",
    "SR_Ca": "SR Ca",
    "K": "K",
    "Na": "Na",
    "HCN": "HCN (If)",
    "transporter": "Na/K-ATPase",
}
CLASS_COLORS = {
    "surface_Ca": "#dd8452",   # orange — antibody-accessible
    "SR_Ca": "#6f9f6f",        # green  — intracellular store
    "K": "#7a86c4",            # blue-purple
    "Na": "#b0b0b0",           # grey
    "HCN": "#c44e52",          # red    — pacemaker funny current
    "transporter": "#937860",  # brown  — pump/transporter
}

# every gene used (for loading)
ALL_CLASS_GENES = tuple(dict.fromkeys(g for c in CLASS_ORDER for g in CHANNEL_CLASSES[c]))


def class_of(gene):
    for c, genes in CHANNEL_CLASSES.items():
        if gene in genes:
            return c
    return None
