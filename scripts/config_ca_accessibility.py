"""
config_ca_accessibility.py — subcellular accessibility of the calcium-handling
machinery, for the autoantibody-target-accessibility analysis.

CURATED INPUT (UniProt subcellular-location field + established channel biology),
not a result. The point: a circulating autoantibody can only reach CELL-SURFACE
(plasma membrane / sarcolemma) channels. Intracellular Ca machinery (SR/ER:
RyR2, IP3R, SERCA) carries the bulk of a cardiomyocyte's Ca flux but is
ANTIBODY-INACCESSIBLE. The antibody-accessible Ca target set is therefore the
surface (plasma-membrane) channels.

Localizations cached at data/ca_localization_cache/subcellular.json (UniProt API).
The voltage-gated Cav family is hard-classified SURFACE: UniProt returns a vague
"Membrane" for some isoforms (CACNA1D/I/B/E), but all Cav1.x/2.x/3.x are
plasma-membrane channels by established electrophysiology, which takes precedence
over the under-specified API field.

access: SURFACE (antibody-accessible) | INTRACELLULAR (SR/ER, inaccessible) | MIXED
"""

# gene -> (access, compartment, ca_role)
CA_MACHINERY = {
    # --- voltage-gated Cav (SURFACE — hard-classified) ---
    "CACNA1C": ("SURFACE", "sarcolemma", "L-type Ca entry (bulk)"),
    "CACNA1D": ("SURFACE", "sarcolemma", "L-type Ca entry (nodal pacemaking)"),
    "CACNA1G": ("SURFACE", "sarcolemma", "T-type Ca entry (CHB epitope target)"),
    "CACNA1H": ("SURFACE", "sarcolemma", "T-type Ca entry"),
    "CACNA1I": ("SURFACE", "sarcolemma", "T-type Ca entry"),
    "CACNA1A": ("SURFACE", "sarcolemma", "P/Q-type (neuronal)"),
    "CACNA1B": ("SURFACE", "sarcolemma", "N-type (neuronal)"),
    "CACNA1E": ("SURFACE", "sarcolemma", "R-type"),
    "CACNA1S": ("SURFACE", "sarcolemma", "L-type (skeletal)"),
    # --- Na-Ca exchanger (SURFACE) ---
    "SLC8A1": ("SURFACE", "sarcolemma", "Na/Ca exchange (fetal Ca entry, immature SR)"),
    # --- store-operated Ca entry (SURFACE, with ER sensor) ---
    "ORAI1": ("SURFACE", "plasma_membrane", "store-operated Ca entry"),
    "ORAI2": ("SURFACE", "plasma_membrane", "store-operated Ca entry"),
    "ORAI3": ("SURFACE", "plasma_membrane", "store-operated Ca entry"),
    "STIM1": ("MIXED", "PM+ER", "SOCE ER-Ca sensor"),
    "STIM2": ("INTRACELLULAR", "ER", "SOCE ER-Ca sensor"),
    # --- TRP cation (mostly SURFACE) ---
    "TRPM4": ("MIXED", "PM+ER", "non-selective cation (conduction)"),
    "TRPM7": ("SURFACE", "plasma_membrane", "non-selective cation"),
    "TRPC1": ("SURFACE", "plasma_membrane", "non-selective cation"),
    "TRPC3": ("SURFACE", "plasma_membrane", "non-selective cation"),
    "TRPC6": ("SURFACE", "plasma_membrane", "non-selective cation"),
    # --- ryanodine receptors (INTRACELLULAR, SR — the bulk Ca release) ---
    "RYR2": ("INTRACELLULAR", "SR", "SR Ca release (dominant cardiac Ca flux)"),
    "RYR3": ("INTRACELLULAR", "SR", "SR Ca release"),
    # --- IP3 receptors (INTRACELLULAR, ER) ---
    "ITPR1": ("INTRACELLULAR", "ER", "ER Ca release (signaling)"),
    "ITPR2": ("INTRACELLULAR", "ER", "ER Ca release (signaling)"),
    "ITPR3": ("INTRACELLULAR", "ER", "ER Ca release (signaling)"),
    # --- SERCA pumps (INTRACELLULAR, SR) ---
    "ATP2A2": ("INTRACELLULAR", "SR", "SR Ca reuptake (SERCA2, matures postnatally)"),
    "ATP2A1": ("INTRACELLULAR", "SR", "SR Ca reuptake (SERCA1, skeletal)"),
    # --- SR counter-ion K channel (INTRACELLULAR; couples to RyR2) ---
    "TMEM38A": ("INTRACELLULAR", "SR", "TRIC-A SR K counter-current (binds RyR2)"),
    "TMEM38B": ("INTRACELLULAR", "ER", "TRIC-B ER K counter-current"),
}


def access(g):      return CA_MACHINERY.get(g, (None,)*3)[0]
def compartment(g): return CA_MACHINERY.get(g, (None,)*3)[1]
def ca_role(g):     return CA_MACHINERY.get(g, (None,)*3)[2]


def genes_by_access(a):
    return [g for g, v in CA_MACHINERY.items() if v[0] == a]
