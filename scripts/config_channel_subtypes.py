"""
config_channel_subtypes.py — functional SUBTYPE taxonomy of cardiac ion channels,
for the developmental subtype-ratio analysis.

CURATED INPUT, not a result. Sources: standard cardiac electrophysiology (Nerbonne
& Kass 2005 Physiol Rev; Amin/Wilde 2010; Bartos/Grandi/Ripplinger 2015 Compr
Physiol). This taxonomy resolves the FUNCTIONAL SUBTYPE within each ion
family (e.g. L-type vs T-type Ca; rapid vs slow delayed-rectifier K), to
measure how the subtype RATIO within each family shifts across developmental stage.

Each gene -> (ion, subtype, current, role_note). `ion` = conducted ion; `subtype`
= the functional channel subtype (the tracked ratio); `current` = the
classical cardiac current it carries; `role_note` = one-line functional role.
"""

# ---------------------------------------------------------------------------
# CALCIUM — the CHB-relevant family. L-type (high-voltage, bulk/contraction) vs
# T-type (low-voltage, transient, early diastolic depolarization / pacemaker).
# ---------------------------------------------------------------------------
CA = {
    "CACNA1C": ("Ca", "L-type", "ICaL", "Cav1.2 — bulk L-current, AP upstroke + EC coupling; ubiquitous"),
    "CACNA1D": ("Ca", "L-type", "ICaL", "Cav1.3 — low-threshold L-current, diastolic depol; SA/AV-node pacemaking"),
    "CACNA1G": ("Ca", "T-type", "ICaT", "Cav3.1 — T-current, early diastolic depol; AV-node enriched (CHB epitope p305)"),
    "CACNA1H": ("Ca", "T-type", "ICaT", "Cav3.2 — T-current, fetal/embryonic-prominent"),
    "CACNA1I": ("Ca", "T-type", "ICaT", "Cav3.3 — T-current, mostly neuronal"),
    "CACNA1S": ("Ca", "L-type", "ICaL", "Cav1.1 — skeletal L-type (control/expected-off in heart)"),
    "CACNA1A": ("Ca", "PQ-type", "ICaP/Q", "Cav2.1 — P/Q-type, neuronal"),
    "CACNA1B": ("Ca", "N-type", "ICaN", "Cav2.2 — N-type, neuronal"),
    "CACNA1E": ("Ca", "R-type", "ICaR", "Cav2.3 — R-type"),
}

# ---------------------------------------------------------------------------
# POTASSIUM — repolarizing currents. The ADULT autoimmune target hERG/KCNH2 lives
# here (rapid delayed-rectifier IKr). Track delayed-rectifier (Kv) vs inward-
# rectifier (Kir) vs transient-outward (Ito) subtypes.
# ---------------------------------------------------------------------------
K = {
    "KCNH2": ("K", "Kv_rapid", "IKr", "Kv11.1/hERG — rapid delayed rectifier (ADULT autoimmune LQT target)"),
    "KCNQ1": ("K", "Kv_slow", "IKs", "Kv7.1 — slow delayed rectifier"),
    "KCNA4": ("K", "Kv_Ito", "Ito", "Kv1.4 — transient outward (fast inactivating)"),
    "KCNA5": ("K", "Kv_ultrarapid", "IKur", "Kv1.5 — ultrarapid delayed rectifier (atrial)"),
    "KCND3": ("K", "Kv_Ito", "Ito", "Kv4.3 — transient outward, main Ito"),
    "KCNJ2": ("K", "Kir", "IK1", "Kir2.1 — inward rectifier, resting potential (ADULT-high)"),
    "KCNJ4": ("K", "Kir", "IK1", "Kir2.3 — inward rectifier"),
    "KCNJ12": ("K", "Kir", "IK1", "Kir2.2 — inward rectifier"),
    "KCNJ11": ("K", "Kir_ATP", "IKATP", "Kir6.2 — K-ATP pore"),
    "KCNJ8": ("K", "Kir_ATP", "IKATP", "Kir6.1 — K-ATP pore"),
    "KCNJ5": ("K", "Kir_GIRK", "IKACh", "Kir3.4 — G-protein-gated (nodal, vagal)"),
    "KCNJ3": ("K", "Kir_GIRK", "IKACh", "Kir3.1 — G-protein-gated (nodal, vagal)"),
}

# ---------------------------------------------------------------------------
# SODIUM — depolarization. Cardiac Nav1.5 (SCN5A) vs neuronal/other.
# ---------------------------------------------------------------------------
NA = {
    "SCN5A": ("Na", "Nav_cardiac", "INa", "Nav1.5 — cardiac fast Na, AP upstroke (adult working-high)"),
    "SCN10A": ("Na", "Nav_neuronal", "INa_late", "Nav1.8 — late Na / conduction modulation"),
    "SCN3B": ("Na", "Nav_beta", "INa_aux", "Navβ3 — auxiliary subunit"),
}

# ---------------------------------------------------------------------------
# HCN — funny current (pacemaker cation). HCN4 nodal vs HCN1/2 distribution.
# ---------------------------------------------------------------------------
HCN = {
    "HCN4": ("cation", "HCN", "If", "HCN4 — dominant cardiac funny current, nodal pacemaker"),
    "HCN1": ("cation", "HCN", "If", "HCN1 — funny current, fetal/nodal"),
    "HCN2": ("cation", "HCN", "If", "HCN2 — funny current"),
    "HCN3": ("cation", "HCN", "If", "HCN3 — funny current (minor)"),
}

SUBTYPES = {**CA, **K, **NA, **HCN}

# Channels that are NOT cardiomyocyte-functional — neuronal-type or skeletal isoforms that
# express at trace levels in heart and sit at the snRNA detection floor. Excluded from the
# cross-platform subtype-RATIO concordance (geneset_capture_qc.py), because a within-family
# ratio that pools these in is dominated by whether a platform happens to capture a near-zero
# transcript, not by cardiac biology. CURATED INPUT (Nerbonne & Kass 2005; Catterall 2011).
NON_CARDIAC_SUBTYPE_GENES = frozenset({
    "CACNA1A",   # Cav2.1 P/Q-type — neuronal
    "CACNA1B",   # Cav2.2 N-type — neuronal
    "CACNA1E",   # Cav2.3 R-type — mostly neuronal
    "CACNA1I",   # Cav3.3 — mostly neuronal
    "CACNA1S",   # Cav1.1 — skeletal
    "SCN10A",    # Nav1.8 — neuronal/late-Na
})
CARDIAC_SUBTYPE_GENES = frozenset(SUBTYPES) - NON_CARDIAC_SUBTYPE_GENES


def is_cardiac_subtype(g):
    return g in CARDIAC_SUBTYPE_GENES


def ion(g):     return SUBTYPES.get(g, (None,)*4)[0]
def subtype(g): return SUBTYPES.get(g, (None,)*4)[1]
def current(g): return SUBTYPES.get(g, (None,)*4)[2]
def role(g):    return SUBTYPES.get(g, (None,)*4)[3]


def genes_for_ion(i):
    return [g for g, v in SUBTYPES.items() if v[0] == i]


def subtypes_for_ion(i):
    out = []
    for g, v in SUBTYPES.items():
        if v[0] == i and v[1] not in out:
            out.append(v[1])
    return out


# Expected developmental directions for each subtype ratio, stated explicitly
# (so the analysis reports confirm/refute, not just numbers):
#   Ca:  T-type/L-type ratio HIGHER in fetal, falls toward adult (T-type switch-off)
#   K:   Kir (IK1) and Kv_rapid (IKr) RISE toward adult (mature repolarization)
#   Na:  Nav_cardiac RISES toward adult (mature fast conduction)
#   HCN: If HIGHER in fetal (immature automaticity), retained only in nodal cells
EXPECTED_SHIFTS = {
    "Ca_Ttype_over_Ltype": "higher fetal -> lower adult",
    "K_Kir_fraction":       "lower fetal -> higher adult",
    "K_IKr_fraction":       "lower fetal -> higher adult (adult autoimmune-LQT target)",
    "Na_cardiac_fraction":  "lower fetal -> higher adult",
    "HCN_fraction":         "higher fetal -> lower adult (retained in nodal)",
}


