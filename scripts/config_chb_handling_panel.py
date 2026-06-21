"""
config_chb_handling_panel.py — gene panel for the node-resolved fetal-vs-adult
ion-channel + Ca-handling expression comparison (SAN / AVN / vCM).

CURATED INPUT, not a result. The panel extends the channel set used elsewhere in
this analysis with the cardiomyocyte/conduction-cell Ca-handling apparatus that the
CHB literature implicates in AV-block pathogenesis beyond the surface Ca channels
themselves.

Primary source for the Ca-handling additions:
  Strandberg / Lazzerini, "... AV block ..." Circ Arrhythm Electrophysiol
  10.1161/CIRCEP.115.003432 — anti-Ro/SSA dysregulation of L-/T-type Ca channels
  AND downstream Ca handling (Na+/Ca2+ exchanger, SERCA2a, RyR2) in fetal AV block.
Channel classification sources are as in config_channel_subtypes.py (Nerbonne &
Kass 2005; Catterall 2011; Amin/Wilde 2010).

Each gene -> (symbol_label, family, role_note). `family` groups the table.
"""
from __future__ import annotations

# (gene, label, family, role)
PANEL = [
    # --- voltage-gated Ca channels (the surface autoantibody targets) ---
    ("CACNA1C", "Cav1.2 (L)",   "Ca channel", "L-type, bulk ICaL / EC coupling; ubiquitous (anti-Ro target)"),
    ("CACNA1D", "Cav1.3 (L)",   "Ca channel", "L-type low-threshold, nodal diastolic depol (anti-Ro target)"),
    ("CACNA1G", "Cav3.1 (T)",   "Ca channel", "T-type, AV-node enriched; CHB epitope p305 (anti-Ro target)"),
    ("CACNA1H", "Cav3.2 (T)",   "Ca channel", "T-type, fetal/embryonic-prominent (anti-Ro target)"),
    ("CACNA1I", "Cav3.3 (T)",   "Ca channel", "T-type, mostly neuronal (reference)"),

    # --- Ca-handling apparatus (CIRCEP.115.003432 downstream effectors) ---
    ("SLC8A1",  "NCX1",         "Ca handling", "Na+/Ca2+ exchanger; nodal pacemaking + diastolic Ca extrusion"),
    ("ATP2A2",  "SERCA2a",      "Ca handling", "SR Ca2+-ATPase 2a; SR reuptake"),
    ("RYR2",    "RyR2",         "Ca handling", "ryanodine receptor 2; SR Ca2+ release / Ca clock"),
    ("PLN",     "PLN",          "Ca handling", "phospholamban; SERCA2a regulator"),
    ("CASQ2",   "CASQ2",        "Ca handling", "calsequestrin-2; SR Ca2+ buffer"),
    ("ATP2B4",  "PMCA4",        "Ca handling", "plasma-membrane Ca2+-ATPase 4; Ca extrusion"),
    ("CALM1",   "CaM1",         "Ca handling", "calmodulin (Ca sensor; regulates Cav/RyR)"),

    # --- K channels (repolarization; adult autoimmune-LQT arm) ---
    ("KCNH2",   "Kv11.1 (IKr)", "K channel",  "rapid delayed rectifier; adult autoimmune-LQT target (hERG)"),
    ("KCNQ1",   "Kv7.1 (IKs)",  "K channel",  "slow delayed rectifier IKs"),
    ("KCNA5",   "Kv1.5 (IKur)", "K channel",  "ultrarapid delayed rectifier (atrial)"),
    ("KCND3",   "Kv4.3 (Ito)",  "K channel",  "transient outward Ito (non-targeted control)"),
    ("KCNJ2",   "Kir2.1 (IK1)", "K channel",  "inward rectifier; resting potential (adult-high)"),
    ("KCNJ3",   "Kir3.1 (IKACh)","K channel", "G-protein-gated; nodal/vagal"),
    ("KCNJ5",   "Kir3.4 (IKACh)","K channel", "G-protein-gated; nodal/vagal (anti-Kir3.4 -> AF)"),

    # --- Na channel + transporter ---
    ("SCN5A",   "Nav1.5",       "Na channel", "cardiac fast Na; AP upstroke (working-high)"),
    ("ATP1A1",  "Na/K-ATPase a1","Transporter","Na/K-ATPase alpha-1; anti-ATP1A1 CHB biomarker (Benjamin 2025)"),

    # --- pacemaker funny current ---
    ("HCN4",    "HCN4 (If)",    "HCN",        "dominant funny current; nodal pacemaker"),
    ("HCN1",    "HCN1 (If)",    "HCN",        "funny current; fetal/nodal"),

    # --- conduction-cell identity markers (orientation, not targets) ---
    ("GJA5",    "Cx40 (GJA5)",  "Identity",   "gap junction; fast-conduction / pacemaker identity"),
    ("GJC1",    "Cx45 (GJC1)",  "Identity",   "gap junction; nodal identity"),
    ("TBX3",    "TBX3",         "Identity",   "conduction-system transcription factor"),
    ("SHOX2",   "SHOX2",        "Identity",   "SAN pacemaker transcription factor"),
    ("ISL1",    "ISL1",         "Identity",   "SAN/AVN pacemaker transcription factor"),
]

PANEL_GENES = [g for g, *_ in PANEL]
LABEL = {g: lab for g, lab, *_ in PANEL}
FAMILY = {g: fam for g, _, fam, _ in PANEL}
ROLE = {g: r for g, *_, r in PANEL}
FAMILY_ORDER = ["Ca channel", "Ca handling", "K channel", "Na channel",
                "Transporter", "HCN", "Identity"]

# ---------------------------------------------------------------------------
# Cell-group definitions, organised into STAGE BLOCKS. Within a block, groups share
# a developmental stage AND a single sequencing platform, so the question "which cell
# type expresses this gene highest" is answered on like-for-like data:
#
#   FETAL block  (snRNA, nuclei): Lim 2024 SAN (GSE279630) + Protze AVN (GSE297072).
#                Cell-type gradient across SAN/AVN x pacemaker/working.
#   ADULT block  (Multiome-RNA, nuclei): Kanemaru SAN_P / AVN_P / vCM. vCM is
#                RESTRICTED to Multiome nuclei so it matches the nodal cells' platform
#                (SAN_P/AVN_P are 100% Multiome nuclei; unrestricted vCM is mostly
#                3'-snRNA -> different chemistry).
#   SIM ladder   (snRNA, whole-CM, no nodal resolution): the ONLY across-dev-stage
#                comparison (fetal -> child -> adult CM).
#
# matcher: ("isin", {labels}) exact-match, ("prefix", "str") startswith, or a synthetic
# "_..." obs_col resolved in the script.
# ---------------------------------------------------------------------------
GROUPS = [
    # --- FETAL block (snRNA) ---
    ("fetal_SAN_pace", "Fetal SAN pacemaker", "lim2024_fetal_san_qc.h5ad",
     "cell_type", ("isin", {"pacemaker_CM"})),
    ("fetal_SAN_work", "Fetal SAN working",   "lim2024_fetal_san_qc.h5ad",
     "cell_type", ("isin", {"working_CM"})),
    ("fetal_AVN_pace", "Fetal AVN pacemaker", "lim2024_fetal_avn_qc.h5ad",
     "cell_type", ("isin", {"pacemaker_CM"})),
    ("fetal_AVN_work", "Fetal AVN working",   "lim2024_fetal_avn_qc.h5ad",
     "cell_type", ("isin", {"working_CM"})),
    # --- ADULT block (Kanemaru, Multiome-RNA nuclei only) ---
    # NODE-MATCHED to the fetal side: each node's pacemaker is contrasted with the working
    # CMs dissected from the SAME node region (SAN-region aCM; AVN-region vCM/aCM), mirroring
    # the fetal node-adjacent working_CM. This makes pacemaker-vs-working a within-node
    # contrast at both stages, symmetric fetal<->adult.
    ("adult_SAN_pace", "Adult SAN pacemaker", "kanemaru_conduction_qc.h5ad",
     "_kan_multiome_SAN_P",   ("isin", {"__TRUE__"})),
    ("adult_SAN_work", "Adult SAN working",   "kanemaru_conduction_qc.h5ad",
     "_kan_nodework_SAN",     ("isin", {"__TRUE__"})),
    ("adult_AVN_pace", "Adult AVN pacemaker", "kanemaru_conduction_qc.h5ad",
     "_kan_multiome_AVN_P",   ("isin", {"__TRUE__"})),
    ("adult_AVN_work", "Adult AVN working",   "kanemaru_conduction_qc.h5ad",
     "_kan_nodework_AVN",     ("isin", {"__TRUE__"})),
    # separate BULK ventricular reference (true ventricle LV/RV/SP/AX, away from the node) —
    # lets the reader compare node working CM against distal ventricular working myocardium.
    ("adult_vCM",      "Adult vCM (bulk vent.)", "kanemaru_conduction_qc.h5ad",
     "_kan_bulk_vCM",         ("isin", {"__TRUE__"})),
    # --- SIM developmental ladder (snRNA whole-CM; the only cross-stage block) ---
    ("sim_fetal_CM", "Sim fetal CM",  "sim2021_qc.h5ad", "_sim_cm_fetal", ("isin", {"__TRUE__"})),
    ("sim_child_CM", "Sim child CM",  "sim2021_qc.h5ad", "_sim_cm_child", ("isin", {"__TRUE__"})),
    ("sim_adult_CM", "Sim adult CM",  "sim2021_qc.h5ad", "_sim_cm_adult", ("isin", {"__TRUE__"})),
]

# Stage blocks for the "highest-expressing cell type" gradient. The Sim ladder is the
# only cross-DEV-stage comparison; the fetal and adult blocks are within-stage,
# cross-CELL-TYPE.
STAGE_BLOCKS = {
    "fetal (snRNA)":         ["fetal_SAN_pace", "fetal_SAN_work", "fetal_AVN_pace", "fetal_AVN_work"],
    "adult (Multiome nuc.)": ["adult_SAN_pace", "adult_SAN_work", "adult_AVN_pace",
                              "adult_AVN_work", "adult_vCM"],
    "Sim ladder (snRNA)":    ["sim_fetal_CM", "sim_child_CM", "sim_adult_CM"],
}

# Z-SCORE SUB-BATCHES for the figure colour encoding. Per-gene colour z-scores must be
# computed WITHIN a single sample/batch so a colour contrast never confounds biology with
# batch. The fetal SAN and AVN are separate GEO samples (GSE279630 vs GSE297072) processed
# as separate batches (likely different days), so they are normalised separately even though
# they share a study group. The adult (single Kanemaru/Multiome object) and Sim (single
# object) blocks are each one batch. Each sub-batch is a list of groups normalised together.
Z_SUBGROUPS = {
    "fetal (snRNA)":         [["fetal_SAN_pace", "fetal_SAN_work"],     # SAN sample
                              ["fetal_AVN_pace", "fetal_AVN_work"]],    # AVN sample (separate batch)
    # adult: all five groups are one Kanemaru/Multiome object (single platform/atlas), so
    # normalised together — the cross-node colour contrast within one atlas is not a batch
    # boundary the way the two separately-deposited fetal GEO accessions are.
    "adult (Multiome nuc.)": [["adult_SAN_pace", "adult_SAN_work", "adult_AVN_pace",
                               "adult_AVN_work", "adult_vCM"]],
    "Sim ladder (snRNA)":    [["sim_fetal_CM", "sim_child_CM", "sim_adult_CM"]],
}


# ---------------------------------------------------------------------------
# DATASET PROVENANCE for the methods table. CURATED INPUT, transcribed from
# DATA_SOURCES.md (the authoritative provenance doc) — accession, stage,
# chemistry, citation. The per-group cell COUNTS are NOT stored here; they are
# computed live from the objects by node_datasets_methods_table.py so the methods
# table cannot drift from the data actually used. Each entry keys on the h5ad file.
#
# `platform_used`: the chemistry of the cells ACTUALLY used by this analysis (which
# may be narrower than the object as a whole — e.g. Kanemaru's nodal cells are
# Multiome-only single nuclei even though the full object is mixed).
# ---------------------------------------------------------------------------
DATASET_PROVENANCE = {
    # Fetal SAN and AVN are two GEO accessions from the same study group, deposited
    # separately: SAN = GSE279630; AVN = GSE297072 (only the dissected fetal-tissue sample
    # GSM8983395 is used, not the hPSC-derived samples in that series).
    "lim2024_fetal_san_qc.h5ad": dict(
        name="Lim/Protze (fetal SAN)", accession="GEO GSE279630 / GSM8577299",
        tissue_stage="human fetal sinoatrial node, GW19",
        chemistry="10x snRNA-seq", platform_used="10x snRNA (single-nucleus)",
        citation="Lim, … Protze; GSE279630 (fetal SAN, CD34-marker study)",
        donors="single fetal donor"),
    "lim2024_fetal_avn_qc.h5ad": dict(
        name="Lim/Protze (fetal AVN)",
        accession="GEO GSE297072 / GSM8983395 (dissected fetal-tissue sample only; NOT the hPSC-derived samples)",
        tissue_stage="human fetal atrioventricular node",
        chemistry="10x snRNA-seq", platform_used="10x snRNA (single-nucleus)",
        citation="Lohbihler, Lim, … Protze; GSE297072 (same group as SAN; deposited 2025)",
        donors="single fetal donor"),
    "kanemaru_conduction_qc.h5ad": dict(
        name="Kanemaru 2023 (adult conduction)",
        accession="Heart Cell Atlas (heartcellatlas.org), conduction object",
        tissue_stage="adult human SAN/AVN + working ventricle",
        chemistry="object as a whole: 10x snRNA + Multiome-RNA + scRNA",
        platform_used=("Multiome-RNA single nuclei ONLY: SAN_P/AVN_P pacemaker states "
                       "are 100% Multiome (SAN/AVN regions sequenced only by Multiome); "
                       "vCM RESTRICTED to Multiome AND to true ventricular regions "
                       "(LV/RV/SP/AX), excluding AVN-dissected working CMs that are also "
                       "labelled vCM"),
        citation="Kanemaru et al. 2023 Nature", donors="multiple adult donors (20-75 yr)"),
    "sim2021_qc.h5ad": dict(
        name="Sim 2021 (dev ladder)", accession="GEO GSE156703",
        tissue_stage="human heart fetal / child / adult CM (within-study dev axis)",
        chemistry="10x snRNA", platform_used="10x snRNA (single-nucleus); whole-CM, no nodal split",
        citation="Sim et al. 2021", donors="3 fetal + 3 child + 3 adult"),
    # Cross-platform corroboration datasets (raw STRT counts read directly, no _qc.h5ad)
    # and the multi-donor first-trimester replication atlas. Provenance transcribed from
    # DATA_SOURCES.md (authoritative).
    "cui2019_gse106118": dict(
        name="Cui 2019 (fetal STRT)", accession="GEO GSE106118",
        tissue_stage="human fetal heart CM, gestational wk 5-25",
        chemistry="STRT / non-10x (Smart-seq-family)", platform_used="STRT (non-droplet)",
        citation="Cui et al. 2019", donors="multiple fetal samples"),
    "wang2020_gse109816": dict(
        name="Wang 2020 (adult STRT)", accession="GEO GSE109816",
        tissue_stage="human adult heart CM (LA + LV)",
        chemistry="STRT / non-10x (Smart-seq-family)", platform_used="STRT (non-droplet)",
        citation="Wang et al. 2020", donors="multiple adult samples"),
    "lazar2025": dict(
        name="Lázár 2025 (HDCA first-trimester)",
        accession="Mendeley Data 10.17632/fhtb99mdzd.1 (open processed); raw under EGA EGAS50000001029",
        tissue_stage="human developing heart, first trimester (pcw 5.5-14)",
        chemistry="10x snRNA", platform_used="10x snRNA (single-nucleus); whole-heart atlas",
        citation="Lázár et al. 2025 (HDCA)", donors="13 first-trimester donors / 21 samples"),
}

# Analysis -> datasets used. One row per analysis/figure for the reviewer-facing
# supplementary "which analysis used which datasets" table. Keys are DATASET_PROVENANCE
# object keys; the figure/role mapping is editorial (transcribed from REPRODUCIBILITY.md
# and the manuscript), not computed.
ANALYSIS_DATASET_MAP = [
    dict(analysis="Pacemaker-vs-working channel contrast (fetal & adult SAN/AVN)",
         figures="Fig 1B,D,E; Fig S3",
         datasets=["lim2024_fetal_san_qc.h5ad", "lim2024_fetal_avn_qc.h5ad",
                   "kanemaru_conduction_qc.h5ad"],
         cell_types="nodal pacemaker_CM vs node-region working_CM, per node (SAN, AVN)",
         method="two-sided Mann-Whitney U + BH-FDR on log1p means (Cav-pool fractions one-sided)",
         comparison_type="within-dataset cell-type contrast (statistics)"),
    dict(analysis="Constitutive-vs-fetal-specific interaction",
         figures="in-text",
         datasets=["lim2024_fetal_san_qc.h5ad", "lim2024_fetal_avn_qc.h5ad",
                   "kanemaru_conduction_qc.h5ad"],
         cell_types="fetal+adult pacemaker & working CM (4 groups per node)",
         method="OLS expression ~ C(celltype)*C(stage); interaction term",
         comparison_type="celltype x stage interaction (batch-robust ratio-of-ratios)"),
    dict(analysis="Fetal/adult target ratio",
         figures="Fig 1C",
         datasets=["sim2021_qc.h5ad"],
         cell_types="whole cardiomyocytes, fetal vs adult",
         method="within-study ratio of mean-per-cell expression",
         comparison_type="within-study cross-stage ratio"),
    dict(analysis="Node-matched fetal-vs-adult contrasts",
         figures="Fig 1; Fig 2; Fig 4",
         datasets=["lim2024_fetal_san_qc.h5ad", "lim2024_fetal_avn_qc.h5ad",
                   "kanemaru_conduction_qc.h5ad"],
         cell_types="node pacemaker & node-region working CM, fetal vs adult, per node",
         method="linear-scale group means, grouped by life stage",
         comparison_type="cross-dataset (direction-only)"),
    dict(analysis="Ca accessibility & channel-class budget",
         figures="Fig 2A,B,C",
         datasets=["lim2024_fetal_san_qc.h5ad", "lim2024_fetal_avn_qc.h5ad",
                   "kanemaru_conduction_qc.h5ad", "sim2021_qc.h5ad"],
         cell_types="node pacemaker/working CM (+ Sim whole-CM fetal/child/adult)",
         method="within-cell surface fraction, NCX:RYR2, six-class budget (internal ratios)",
         comparison_type="within-cell internal ratios (+ Sim within-study trajectory)"),
    dict(analysis="Developmental hand-off (Cav dominance, RYR2, inward:outward)",
         figures="Fig 3",
         datasets=["sim2021_qc.h5ad"],
         cell_types="whole cardiomyocytes, fetal -> child -> adult (3 donors/stage)",
         method="within-family dominance, per-stage means, inward:outward composition ratio",
         comparison_type="within-study cross-stage axis (3 donors/stage)"),
    dict(analysis="First-trimester donor-level replication",
         figures="Fig S2C",
         datasets=["lazar2025"],
         cell_types="whole cardiomyocytes, per donor (13 donors)",
         method="per-donor RYR2 / Ca T(T+L) vs gestational age, Spearman",
         comparison_type="donor-level Spearman (13 donors)"),
    dict(analysis="Cross-platform corroboration & direction agreement",
         figures="Fig S2A,B",
         datasets=["sim2021_qc.h5ad", "cui2019_gse106118", "wang2020_gse109816"],
         cell_types="whole/working cardiomyocytes, stage-matched (fetal Cui, adult Wang)",
         method="Spearman on log10 expression; fetal/adult log2-ratio direction agreement",
         comparison_type="cross-platform corroboration (direction)"),
    dict(analysis="Platform QC (complexity, depth-matching, gene-set/channel capture)",
         figures="Fig S1",
         datasets=["lim2024_fetal_san_qc.h5ad", "lim2024_fetal_avn_qc.h5ad",
                   "kanemaru_conduction_qc.h5ad", "sim2021_qc.h5ad",
                   "cui2019_gse106118", "wang2020_gse109816"],
         cell_types="working CM (stage-matched) for capture; pacemaker/working for depth-match",
         method="median genes/cell; depth-band-matched FC; GO gene-set & per-channel capture",
         comparison_type="QC / robustness checks"),
]
