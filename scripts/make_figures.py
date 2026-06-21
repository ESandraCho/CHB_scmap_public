#!/usr/bin/env python3
"""
make_figures.py — render the CHB paper figures from the committed result tables.
Reproducible: reads only results/tables/*.{csv,json}; matplotlib only.

Developmental panels are on the SMART-SEQ pipeline (Sim2021 within-study, corroborated by
non-10x Cui2019/Wang2020 STRT) via fetal_vs_adult_vulnerability.json and
ca_subtype_dominance_dev.json. Dissected-nodal panels (Lim/Kanemaru) are unchanged.

seed: 42 (no randomness; layout deterministic)

Figures built here (F4 is built by make_node_handling_figure.py):
  F1  nodal channel program: within-fetal pacemaker-vs-working SAN/AVN (A,B)
      + Sim within-study fetal/adult target ratios (C) + AVN cell-type contrasts (D,E)
  F2  accessibility: surface-accessible Ca fraction (A) + NCX:RYR2 mode (B)
      + six-class channel-class budget (C) — internal per-cell ratios/compositions
  F3  developmental hand-off of Ca reliance (Sim within-study, whole-CM): surface Cav subtypes
      (A) + RYR2 rise (B) + inward:outward composition balance (C)
  S1  platform QC: per-cell complexity + nodal depth-matching + gene-set capture
      + per-gene pattern + per-channel capture. Reads platform_qc.json + geneset_capture_qc.json
  S2  dataset reliability: multi-platform Sim/Cui/Wang agreement + Lázár per-donor + cell counts
  S3  SAN pacemaker-vs-working (the SAN equivalents of F1D,E)

Output: results/figures/*.png (300 dpi) + *.svg (vector, editable <text>)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg
import config_ca_accessibility as acc
import config_channel_classes as cc

T = cfg.OUT_TABLES
FIG = cfg.PUBLIC / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# --- palette ---
FETAL = "#c44e52"     # red
ADULT = "#4c72b0"     # blue
WORK = "#8c8c8c"      # grey
ACCENT = "#dd8452"    # orange (T-type / Cav3.1)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight", "axes.titlesize": 10,
                     "axes.titleweight": "bold",
                     # keep SVG text as editable <text> elements (not glyph paths),
                     # so labels stay selectable/editable in Illustrator/Inkscape
                     "svg.fonttype": "none"})


def _save(fig, stem):
    """Write a figure as both PNG (300 dpi raster) and SVG (vector, editable text)."""
    fig.savefig(FIG / f"{stem}.png")
    fig.savefig(FIG / f"{stem}.svg")


def read(name):
    return pd.read_csv(T / f"{name}.csv", comment="#")


def read_json(name):
    """Load a committed JSON result table (Smart-seq pipeline outputs)."""
    import json
    return json.loads((T / f"{name}.json").read_text())


def star(p_or_fdr):
    return "*" if (pd.notna(p_or_fdr) and p_or_fdr < 0.05) else ""


# =================================================================
# Panel drawers (take an Axes). Merged main figures compose these so
# each panel's plotting logic lives in exactly one place.
# =================================================================

def _panel_schematic(ax, letter):
    """Organizing schematic: which autoimmune-channelopathy PHENOTYPE targets which CELL
    TYPE and which channels. Two clinical syndromes hit two cell types because each targets
    channels enriched in that cell type. Colour = phenotype stage (fetal/adult/both/none),
    hatch = arrhythmia type (brady/tachy/both); bold = CHB target. Encoding derived from the
    cited channelopathy_classification.json (Lazzerini & Boutjdir 2025 Fig 1)."""
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.text(5, 9.6, "Autoimmune channelopathy: each antibody phenotype targets the cell type "
                    "that expresses its channels", ha="center", fontsize=8, fontweight="bold")

    def box(x, y, w, h, fc, ec="black"):
        ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=ec, lw=1.0, zorder=1))

    fcol = cfg.PHENOTYPE_STAGE_COLOR
    # left: bradyarrhythmia -> conduction / pacemaker cell
    box(0.3, 5.6, 4.5, 2.7, "#f6e4e4")
    ax.text(2.55, 7.95, "BRADYARRHYTHMIA (AV block)", ha="center", fontsize=7.5, fontweight="bold", color="#a02020")
    ax.text(2.55, 7.4, "fetal — congenital heart block", ha="center", fontsize=6, color=fcol["fetal"])
    ax.text(2.55, 6.85, "→ conduction / pacemaker cell", ha="center", fontsize=6.5, style="italic")
    ax.text(2.55, 6.25, "Cav3.1 · Cav3.2 · Cav1.3", ha="center", fontsize=6.5)
    # right: tachyarrhythmia -> working CM
    box(5.2, 5.6, 4.5, 2.7, "#e4eaf4")
    ax.text(7.45, 7.95, "TACHYARRHYTHMIA (LQTS / VF)", ha="center", fontsize=7.5, fontweight="bold", color="#2050a0")
    ax.text(7.45, 7.4, "adult — LQTS / SQTS / AF / VF-CA", ha="center", fontsize=6, color=fcol["adult"])
    ax.text(7.45, 6.85, "→ working cardiomyocyte", ha="center", fontsize=6.5, style="italic")
    ax.text(7.45, 6.25, "hERG · Kv7.1 · Kv1.4 · Kir3.4 · Kir6.2", ha="center", fontsize=6.5)
    # middle-bottom: 'both' channels + non-targeted
    box(0.3, 3.4, 4.5, 1.6, "#faf0e6")
    ax.text(2.55, 4.6, "Cav1.2 (CACNA1C)  &  Nav1.5 (SCN5A) — BOTH", ha="center", fontsize=6.4, fontweight="bold")
    ax.text(2.55, 4.0, "brady in node cells + tachy in ventricular CMs", ha="center", fontsize=5.6, color=fcol["both"])
    ax.text(2.55, 3.6, "(epitope- & cell-type-dependent)", ha="center", fontsize=5.4, style="italic", color="#777")

    # encoding key (used by all Fig-1 panels) — fixed positions to avoid overlap
    ax.text(0.3, 2.6, "Encoding (all panels):", fontsize=6.5, fontweight="bold")
    # colour = stage (left half, fixed slots). "not targeted"/none omitted — Fig 1 shows
    # antibody-target channels only (the non-targeted KCND3 was removed).
    for (lab, col), sx in zip([("fetal", fcol["fetal"]), ("adult", fcol["adult"]),
                               ("both", fcol["both"])],
                              [0.3, 1.55, 2.8]):
        ax.add_patch(plt.Rectangle((sx, 1.75), 0.35, 0.35, fc=col, ec="black", lw=0.4))
        ax.text(sx + 0.45, 1.93, lab, fontsize=5.4, va="center")
    ax.text(0.3, 1.25, "colour = phenotype stage", fontsize=5.4, color="#444")
    # hatch = arrhythmia type (right half) — only real phenotypes; non-targeted genes carry no
    # hatch (no antibody -> no phenotype), flagged by the grey "not targeted" colour instead.
    for (lab, h), hx in zip([("brady", cfg.PHENOTYPE_TYPE_HATCH["brady"]),
                             ("tachy", cfg.PHENOTYPE_TYPE_HATCH["tachy"]),
                             ("both", cfg.PHENOTYPE_TYPE_HATCH["both"])],
                            [5.6, 6.9, 8.2]):
        ax.add_patch(plt.Rectangle((hx, 1.75), 0.35, 0.35, fc="white", ec="black", lw=0.5, hatch=h))
        ax.text(hx + 0.42, 1.93, lab, fontsize=5.2, va="center")
    ax.text(5.6, 1.25, "hatch = arrhythmia type (brady/tachy/both)  ·  bold = CHB target", fontsize=5.3, color="#444")
    ax.set_title(f"{letter} · Antibody phenotype → cell type → channels", fontsize=8.5, loc="left")


def _panel_pace_work(ax, node, letter, stage="fetal"):
    """Within-dataset pacemaker-vs-working log2FC for one node, at one stage (fetal or
    adult). BARS encode enrichment DIRECTION (pacemaker-enriched vs depleted — the panel's own
    signal). The phenotype scheme is on the GENE LABELS: label colour = phenotype stage, bold =
    CHB target. fetal -> fetal_pacemaker_vs_working; adult -> adult_pacemaker_vs_working."""
    d = read(f"{stage}_pacemaker_vs_working")
    order = list(reversed(cfg.FIG1_PANEL_GENES)) + list(cfg.FIG1_IDENTITY_GENES)
    d = d[d.gene.isin(order)].copy()
    d["gene"] = pd.Categorical(d["gene"], order); d = d.sort_values("gene")
    sub = d[d.node == node]
    y = np.arange(len(sub))
    vals = sub["log2FC_pace_vs_work"].values
    ax.barh(y, vals, color=["#b9856a" if v > 0 else "#6a89b9" for v in vals],
            edgecolor="black", linewidth=0.4)
    ax.set_yticks(y); ax.set_yticklabels([cfg.channel_label(g) for g in sub["gene"].astype(str)])
    for tick, g in zip(ax.get_yticklabels(), sub["gene"].astype(str)):
        tick.set_color(cfg.phenotype_stage_color(g))
        tick.set_fontweight("bold" if cfg.is_chb_target(g) else "normal")
    ax.axvline(0, color="black", lw=0.8)
    has_fdr = "fdr" in sub.columns
    for yi, v in enumerate(vals):
        s = star(sub["fdr"].iloc[yi]) if has_fdr else ""
        ax.text(v + (0.12 if v >= 0 else -0.12), yi, f"{v:+.1f}{s}",
                va="center", ha="left" if v >= 0 else "right", fontsize=7)
    ax.set_title(f"{letter} · {stage.capitalize()} {node}: pacemaker vs working CM", fontsize=9, pad=6)
    ax.set_xlabel("log2FC (pacemaker-enriched → / ← depleted)")
    ax.set_xlim(-2.8, 6.0 if stage == "adult" else 4.4)


# Cav-subtype composition panel (Fig 1B): the 4 cardiac Cav subtypes as fractions of their sum,
# per node group (stage x cell type). Reads node_channel_handling_expression.csv.
_CAV4 = [("CACNA1C", "Cav1.2 (CACNA1C)", ACCENT),
         ("CACNA1D", "Cav1.3 (CACNA1D)", FETAL),
         ("CACNA1G", "Cav3.1 (CACNA1G)", "#7a4fa0"),
         ("CACNA1H", "Cav3.2 (CACNA1H)", "#3aa37a")]
# node groups, blocked by stage (fetal | adult); 4 cell types each.
_CAV_GROUPS = [
    ("fetal", [("fetal_SAN_pace", "SAN\npace"), ("fetal_SAN_work", "SAN\nwork"),
               ("fetal_AVN_pace", "AVN\npace"), ("fetal_AVN_work", "AVN\nwork")]),
    ("adult", [("adult_SAN_pace", "SAN\npace"), ("adult_SAN_work", "SAN\nwork"),
               ("adult_AVN_pace", "AVN\npace"), ("adult_AVN_work", "AVN\nwork")]),
]


def _panel_cav_subtype_stack(ax, letter):
    """Within-cell COMPOSITION of the four cardiac Cav subtypes (Cav1.2/1.3/3.1/3.2) as
    fractions of their summed expression, per node group (stage x cell type). Shows how the
    Ca-channel program's subtype balance shifts across cell type and stage — e.g. the
    L-type-bulk Cav1.2 share vs the conduction specialists Cav1.3/Cav3.1. Stacked bars,
    blocked fetal | adult; node-resolved (Lim/Protze fetal, Kanemaru adult). Each bar is an
    internal within-cell composition (depth-robust); cross-stage is cross-dataset (direction)."""
    d = read("node_channel_handling_expression").set_index(["group", "gene"])["mean_tp10k"]
    # significance of pacemaker>working for each Cav subtype's COMPOSITION fraction (one-sided
    # MWU, BH-FDR; cav_subtype_pace_work_stats.py). Keyed (pair, gene); applies to pacemaker bars.
    try:
        st = read("cav_subtype_pace_work_stats")
        sig = {(r["pair"], r["gene"]): bool(r["sig_fdr_0.05"]) for _, r in st.iterrows()}
    except Exception:
        sig = {}
    # group key (e.g. fetal_SAN_pace) -> stats pair label (e.g. "fetal SAN"); only pacemaker bars
    def pair_of(group_key):
        parts = group_key.split("_")            # [stage, node, pace|work]
        return (f"{parts[0]} {parts[1]}", parts[2] == "pace")

    pos, ticks, labels, block_first = 0.0, [], [], {}
    first = True
    for stage, groups in _CAV_GROUPS:
        for g, lab in groups:
            pair_label, is_pace = pair_of(g)
            vals = np.array([float(d.get((g, gene), 0.0)) for gene, *_ in _CAV4])
            tot = vals.sum()
            fr = vals / tot if tot > 0 else vals
            bottom = 0.0
            for (gene, glab, col), v in zip(_CAV4, fr):
                ax.bar(pos, v, 0.8, bottom=bottom, color=col, edgecolor="white", lw=0.4,
                       label=glab if first else None)
                if v > 0.10:
                    ax.text(pos, bottom + v / 2, f"{v:.0%}", ha="center", va="center",
                            fontsize=5.6, color="white", fontweight="bold")
                # asterisk on the Cav1.3/Cav3.1 segment of a pacemaker bar if its fraction is
                # significantly > the matched working bar (FDR<0.05, pace>work).
                if is_pace and gene in ("CACNA1D", "CACNA1G") and sig.get((pair_label, gene)):
                    ax.text(pos, bottom + v / 2, "*", ha="center", va="center",
                            fontsize=11, color="white", fontweight="bold")
                bottom += v
            first = False
            ticks.append(pos); labels.append(lab)
            block_first.setdefault(stage, []).append(pos)
            pos += 1.0
        pos += 0.6  # gap between fetal | adult blocks
    ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylim(0, 1.0); ax.set_ylabel("fraction of Cav-subtype pool")
    for stage, ps in block_first.items():
        ax.text(sum(ps) / len(ps), -0.20, stage, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8, fontweight="bold",
                color=FETAL if stage == "fetal" else ADULT)
    ax.legend(fontsize=6.2, frameon=False, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, 1.005), handlelength=1.0, columnspacing=0.9)
    ax.set_title(f"{letter} · Cav-subtype composition by cell type & stage", fontsize=8.5, pad=16)
    ax.text(0.0, -0.30, "* Cav1.3/Cav3.1 fraction higher in pacemaker vs matched working CM "
            "(FDR<0.05, 1-sided MWU)", transform=ax.transAxes, fontsize=5.4, color="#555",
            ha="left", va="top")


def _panel_fetal_ratio(ax, letter):
    """Sim2021 within-study fetal/adult CM expression ratio. Bars are coloured by the
    antibody's CLINICAL phenotype STAGE (fetal/adult/both/none — not the expression conclusion) and
    hatched by arrhythmic phenotype TYPE (brady / tachy / both / none, per Lazzerini 2025
    Fig 1); CHB-target genes have BOLD labels. Shows the antibody-TARGET channels only. Target
    ratios from the vulnerability JSON, else the Ca dominance JSON (same Sim within-study
    source, verified identical for shared genes)."""
    v = read_json("fetal_vs_adult_vulnerability")["targets"]
    dom = read_json("ca_subtype_dominance_dev")["per_gene"]
    eps = 1e-3
    # SHARED Fig-1 channel set (same as 1B,C). Each ratio comes from the vulnerability JSON,
    # else the Ca dominance JSON.
    genes, ratios = [], []
    for g in cfg.FIG1_PANEL_GENES:
        if g in cfg.FIG1_RATIO_EXCLUDE:
            continue
        if g in v and v[g].get("sim_fetal_over_adult") is not None:
            genes.append(g); ratios.append(v[g]["sim_fetal_over_adult"])
        elif g in dom:
            f, a = dom[g]["fetal"]["mean_tp10k"], dom[g]["adult"]["mean_tp10k"]
            genes.append(g); ratios.append(round((f + eps) / (a + eps), 3))
    x = np.arange(len(genes))
    bars = ax.bar(x, ratios, color=[cfg.phenotype_stage_color(g) for g in genes],
                  edgecolor="black", lw=0.5)
    for b, g in zip(bars, genes):                       # hatch = clinical phenotype TYPE
        b.set_hatch(cfg.phenotype_type_hatch(g))
    ax.axhline(1, color="black", lw=0.8, ls="--")
    for xi, r in zip(x, ratios):
        ax.text(xi, r + 0.12, f"{r:.1f}", ha="center", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels([cfg.channel_label(g, short=True) for g in genes],
                       rotation=30, ha="right", fontsize=7)
    for tick, g in zip(ax.get_xticklabels(), genes):    # BOLD label = CHB target
        tick.set_fontweight("bold" if cfg.is_chb_target(g) else "normal")
    ax.set_ylabel("Sim fetal / adult CM (within-study)")
    ax.set_title(f"{letter} · Stage-ratio by antibody phenotype (bold = CHB target)", fontsize=8)
    # legend: colour = stage, hatch = type. ("not targeted"/none is omitted — the only
    # non-targeted gene, KCND3, is no longer shown in Fig 1, so all bars are antibody targets.)
    stage_h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in
               (cfg.PHENOTYPE_STAGE_COLOR["fetal"], cfg.PHENOTYPE_STAGE_COLOR["adult"],
                cfg.PHENOTYPE_STAGE_COLOR["both"])]
    leg1 = ax.legend(stage_h, ["fetal disease", "adult disease", "both"],
                     fontsize=5.5, frameon=False, loc="upper right", title="phenotype stage (colour)",
                     title_fontsize=5.5)
    ax.add_artist(leg1)
    # arrhythmia-type hatch legend — only real phenotypes; non-targeted genes have no
    # phenotype (no antibody) and are flagged by the grey "not targeted" colour, not a hatch.
    type_h = [plt.Rectangle((0, 0), 1, 1, fc="white", ec="black", hatch=h) for h in
              (cfg.PHENOTYPE_TYPE_HATCH["brady"], cfg.PHENOTYPE_TYPE_HATCH["tachy"],
               cfg.PHENOTYPE_TYPE_HATCH["both"])]
    ax.legend(type_h, ["brady", "tachy", "both"],
              fontsize=5.5, frameon=False, loc="upper center", title="arrhythmia type (hatch)",
              title_fontsize=5.5)


# The "pacemaker-enriched at both stages" (constitutive) finding is read directly from the
# fetal-AVN (D) and adult-AVN (E) pacemaker-vs-working panels. The formal cell-type x stage
# interaction backing it is in nodal_specificity_interaction.py (cited in Results/Methods text).


# Bars are ordered and BLOCKED BY DATASET so the load-bearing comparisons are the
# WITHIN-dataset, cross-cell-type ones (clean of platform / abundance-dependent
# capture bias): all-Lim fetal block (SAN, AVN, node wCM) | all-Kanemaru adult block
# (SAN, AVN, node wCM). The cross-stage (fetal-vs-adult) gap straddles the divider and
# is corroborating, not load-bearing — it carries the cross-dataset capture caveat
# (see limitations). The rightmost bar of each block is the NODE-REGION WORKING CM
# (working myocytes dissected within the SAN/AVN block — NOT bulk ventricle; the fetal
# datasets are dissected node tissue, and the adult bar is node-matched to it). It keeps
# the WORK grey fill with a stage-coloured edge.
# NODE-RESOLVED bars, matching F_node_channel_handling: SAN/AVN x pacemaker/working per
# stage, plus an adult bulk-ventricular reference. Bars are separated by LIFE-STAGE
# (fetal | adult) — NOT by dataset name (the old "Lim | Kanemaru" framing was stale now
# that the groups are platform-consistent and node-matched). Pacemaker bars take the stage
# colour; working bars are grey with a stage-coloured edge; the adult bulk-vCM reference is
# a hatched grey bar (no fetal counterpart — no fetal ventricle exists in these data).
_ACC_FETAL = ["fetal_SAN_pace", "fetal_SAN_work", "fetal_AVN_pace", "fetal_AVN_work"]
_ACC_ADULT = ["adult_SAN_pace", "adult_SAN_work", "adult_AVN_pace", "adult_AVN_work",
              "adult_vCM_bulk"]
_ACC_GROUPS = _ACC_FETAL + _ACC_ADULT
_ACC_LABELS = ["SAN\npace", "SAN\nwork", "AVN\npace", "AVN\nwork",
               "SAN\npace", "SAN\nwork", "AVN\npace", "AVN\nwork", "vCM\n(bulk)"]
_ACC_COLS = [FETAL, WORK, FETAL, WORK,  ADULT, WORK, ADULT, WORK, WORK]
_ACC_EDGES = ["black", FETAL, "black", FETAL,  "black", ADULT, "black", ADULT, ADULT]
_ACC_EDGEW = [0.4, 1.4, 0.4, 1.4,  0.4, 1.4, 0.4, 1.4, 1.4]
_ACC_HATCH = [None, None, None, None,  None, None, None, None, "////"]  # bulk vCM hatched
_ACC_DIVIDER = 3.5  # x-position of the fetal | adult stage boundary (after the 4 fetal bars)


# SAN | AVN sub-dividers, within each stage block. The fetal SAN (Lim, GSE279630) and fetal
# AVN (Protze, GSE297072) are SEPARATE samples/dissections (different timelines), so the
# SAN-vs-AVN gap is a cross-sample boundary even within the fetal block; the adult SAN and AVN
# are likewise distinct node dissections. Drawn fainter than the fetal|adult stage divider.
_ACC_SAN_AVN_DIVIDERS = (1.5, 5.5)   # fetal SAN|AVN (after pos 1); adult SAN|AVN (after pos 5)


def _acc_stage_divider(ax):
    """Vertical dividers: a bold fetal | adult life-stage boundary, plus fainter SAN | AVN
    sub-boundaries within each stage (separate samples / dissection timelines). Comparisons
    WITHIN a node and stage are platform-consistent; gaps across SAN/AVN or across stage are
    cross-sample/cross-dataset (capture caveat — see limitations)."""
    for xd in _ACC_SAN_AVN_DIVIDERS:                       # SAN | AVN sub-dividers (faint)
        ax.axvline(xd, color="#bbb", lw=0.7, ls=(0, (2, 2)), zorder=0)
    ax.axvline(_ACC_DIVIDER, color="#888", lw=0.9, ls=(0, (4, 3)), zorder=0)
    yt = ax.get_ylim()[1]
    # stage labels placed clear of the SAN|AVN sub-dividers (at x=1.5 / 5.5)
    ax.text(0.5, yt * 0.995, "fetal", fontsize=6.8, color=FETAL,
            ha="center", va="top", style="italic", fontweight="bold")
    ax.text(6.5, yt * 0.995, "adult", fontsize=6.8, color=ADULT,
            ha="center", va="top", style="italic", fontweight="bold")


def _acc_bars(ax, vals):
    x = np.arange(len(_ACC_GROUPS))
    for xi, v, c, e, w, h in zip(x, vals, _ACC_COLS, _ACC_EDGES, _ACC_EDGEW, _ACC_HATCH):
        ax.bar(xi, v, color=c, edgecolor=e, linewidth=w, hatch=h)
    ax.set_xticks(x); ax.set_xticklabels(_ACC_LABELS, fontsize=6.5)
    return x


def _acc_percell_violins(ax, table, value_col, scale=1.0, clip=None):
    """Overlay the per-cell distribution (within-sample cell dispersion, NOT donor-level
    error) as a thin violin per group, on the same x positions as the bars. Values can be
    scaled (e.g. fraction->%) and clipped to the panel range so a heavy tail doesn't blow up
    the axis (clipping is cosmetic; the bar is the full-data summary)."""
    d = read(table)
    for xi, g in enumerate(_ACC_GROUPS):
        v = d[d.group == g][value_col].to_numpy(float) * scale
        if clip is not None:
            v = v[(v >= clip[0]) & (v <= clip[1])]
        if v.size < 5:
            continue
        parts = ax.violinplot([v], positions=[xi], widths=0.7, showextrema=False,
                              showmedians=False)
        for b in parts["bodies"]:
            b.set_facecolor("#333333"); b.set_alpha(0.16); b.set_edgecolor("none")
            b.set_zorder(2)


def _panel_surface_frac(ax, letter):
    d = read("ca_accessibility")
    surf = []
    for g in _ACC_GROUPS:
        sub = d[d.group == g]; tot = sub["expr"].sum()
        surf.append(100 * sub[sub.access == "SURFACE"]["expr"].sum() / tot)
    x = _acc_bars(ax, surf)
    _acc_percell_violins(ax, "ca_accessibility_percell", "surface_frac", scale=100.0)
    # 50% reference: above = majority of the Ca machinery is antibody-accessible (surface);
    # below = majority is hidden on the intracellular SR. Fetal groups sit above, adult below.
    ax.axhline(50, color="black", lw=0.6, ls="--", zorder=1)
    for xi, s in zip(x, surf):
        ax.text(xi, s + 1.5, f"{s:.0f}%", ha="center", fontsize=6.8)
    ax.set_ylabel("surface-accessible % of Ca machinery")
    ax.set_title(f"{letter} · Antibody-accessible Ca fraction", fontsize=9); ax.set_ylim(0, 100)
    _acc_stage_divider(ax)


def _panel_ncx_mode(ax, letter):
    prox = read("sr_activity_proxies").set_index("group")
    ncx = [prox.loc[g, "NCX_over_RYR2"] for g in _ACC_GROUPS]
    x = _acc_bars(ax, ncx)
    ymax = max(ncx) * 1.5
    # per-cell NCX:RYR2 is right-skewed; clip the violin to the panel range (cosmetic — the
    # bar remains the full-data per-cell-mean ratio).
    _acc_percell_violins(ax, "sr_activity_percell", "ncx_over_ryr2", clip=(0, ymax))
    ax.axhline(1, color="black", lw=0.6, ls="--")
    for xi, v in zip(x, ncx):
        ax.text(xi, v + 0.06, f"{v:.2f}", ha="center", fontsize=6.8)
    ax.set_ylabel("NCX : RYR2  (surface ↔ SR mode)")
    ax.set_title(f"{letter} · Ca-handling mode (~5× flip)", fontsize=9)
    ax.set_ylim(0, ymax)
    ax.text(1.5, max(ncx) * 1.16, "surface-Ca mode", fontsize=6.6, color=FETAL, ha="center")
    ax.text(6.0, 0.72, "SR mode", fontsize=6.6, color=ADULT, ha="center")
    _acc_stage_divider(ax)


# Cell-type blocks for the channel-class composition panel. Each block is
# (block_label, [(stage_label, group_key), ...]). The NODE blocks (SAN/AVN x pace/work)
# are fetal->adult PAIRS, node-matched (node-region working CM, NOT bulk ventricle),
# matching F_node_channel_handling. The final SIM block is the whole-CM developmental
# TRAJECTORY (fetal->child->adult) from the single-platform Sim2021 dataset — added so the
# budget shift is shown on a clean within-study dev axis alongside the node data.
_CB_BLOCKS = [
    ("SAN pace", [("fetal", "fetal_SAN_pace"), ("adult", "adult_SAN_pace")]),
    ("SAN work", [("fetal", "fetal_SAN_work"), ("adult", "adult_SAN_work")]),
    ("AVN pace", [("fetal", "fetal_AVN_pace"), ("adult", "adult_AVN_pace")]),
    ("AVN work", [("fetal", "fetal_AVN_work"), ("adult", "adult_AVN_work")]),
    ("Sim whole-CM", [("fetal", "sim_fetal_CM"), ("child", "sim_child_CM"),
                      ("adult", "sim_adult_CM")]),
]


def _panel_class_budget(ax, letter, horizontal=False):
    """Within-cell ion-channel CLASS composition (surface Ca / SR Ca / K / Na / HCN /
    Na-K-ATPase) as a fraction of the class pool. Node blocks are fetal->adult pairs
    (SAN/AVN x pacemaker/working, node-matched); the final block is the Sim whole-CM
    fetal->child->adult developmental trajectory (single within-study platform). Each bar is
    an internal within-group composition (depth-robust); the cross-stage shift within the
    node blocks is cross-dataset (direction-only), while the Sim block is within-study.
    horizontal=True draws horizontal stacked bars so the thin K/Na/HCN/transporter slices
    have room to read."""
    d = read("channel_class_budget")
    frac = d.pivot(index="group", columns="cls", values="frac")
    pos, ticks, labels = 0.0, [], []
    block_first, first_bar = {}, True
    for ct, stages in _CB_BLOCKS:
        for stage, g in stages:
            start = 0.0
            for c in cc.CLASS_ORDER:
                v = float(frac.loc[g, c])
                lab = cc.CLASS_LABELS[c] if first_bar else None
                if horizontal:
                    ax.barh(pos, v, 0.8, left=start, color=cc.CLASS_COLORS[c],
                            edgecolor="white", linewidth=0.5, label=lab)
                else:
                    ax.bar(pos, v, 0.8, bottom=start, color=cc.CLASS_COLORS[c],
                           edgecolor="white", linewidth=0.5, label=lab)
                if c in ("surface_Ca", "SR_Ca") and v > 0.08:
                    cx, cy = (start + v / 2, pos) if horizontal else (pos, start + v / 2)
                    ax.text(cx, cy, f"{v:.0%}", ha="center", va="center",
                            fontsize=6.2, color="white", fontweight="bold")
                start += v
            first_bar = False
            ticks.append(pos); labels.append(stage)
            block_first.setdefault(ct, []).append(pos)
            pos += 1.0
        pos += 0.6  # gap between blocks

    if horizontal:
        ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=7.0)
        ax.invert_yaxis()                       # first bar at top
        ax.set_xlim(0, 1.0); ax.set_xlabel("fraction of channel-class pool")
        for ct, ps in block_first.items():      # block label centred on its stage rows
            ax.text(-0.115, sum(ps) / len(ps), ct, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=7.2, fontweight="bold", rotation=0)
        ax.legend(fontsize=6.6, frameon=False, ncol=len(cc.CLASS_ORDER), loc="lower center",
                  bbox_to_anchor=(0.5, 1.01), handlelength=1.0, columnspacing=0.9)
    else:
        ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=7.0)
        ax.set_ylim(0, 1.0); ax.set_ylabel("fraction of channel-class pool")
        for ct, ps in block_first.items():
            ax.text(sum(ps) / len(ps), -0.16, ct, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8, fontweight="bold")
        ax.legend(fontsize=6.2, frameon=False, ncol=3, loc="lower center",
                  bbox_to_anchor=(0.5, 1.005), handlelength=1.0, columnspacing=0.9)
    ax.set_title(f"{letter} · Channel-class budget: node (fetal→adult) + Sim dev trajectory",
                 fontsize=9, pad=18)


# Sim 2021 within-study developmental stages (Fig 3 dev-trajectory panels).
_DOM_STAGES = ["fetal", "child", "adult"]


def _panel_sr_traj(ax, letter):
    """RYR2 (SR Ca-release) rises postnatally — Sim2021 within-study, Smart-seq pipeline."""
    dom = read_json("ca_subtype_dominance_dev")
    pg = dom["per_gene"]["RYR2"]; byd = dom.get("per_gene_by_donor", {}).get("RYR2", {})
    vals = [pg[s]["mean_tp10k"] for s in _DOM_STAGES]
    ax.plot(range(3), vals, "-o", color=ADULT, lw=2, zorder=3)
    for i, s in enumerate(_DOM_STAGES):            # donor-level dispersion (3 Sim donors/stage)
        _donor_points(ax, i, byd.get(s, {}), color=ADULT)
    for i, v in enumerate(vals):
        ax.text(i, v + 8, f"{v:.0f}", ha="center", fontsize=7.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(_DOM_STAGES)
    ax.set_ylabel("RYR2 (SR release) mean TP10K")
    allv = vals + [v for s in _DOM_STAGES for v in byd.get(s, {}).values()]
    ax.set_ylim(0, max(allv) * 1.2)
    ax.set_title(f"{letter} · SR (RYR2) matures postnatally", fontsize=9)


def _donor_points(ax, xpos, donor_dict, color="#222", jitter=0.0):
    """Overlay per-donor values at x=xpos as small points + a vertical SD/range whisker
    (donor-level dispersion — Sim has 3 donors/stage, genuine biological replication, NOT
    pseudoreplicated cells). donor_dict maps donor_id -> value."""
    vals = [v for v in donor_dict.values() if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if not vals:
        return
    vals = np.array(vals, float)
    # range whisker (min..max of the donors) — honest for n=3, no normality assumption
    ax.plot([xpos, xpos], [vals.min(), vals.max()], color=color, lw=0.8, zorder=4)
    xs = xpos + (np.linspace(-1, 1, len(vals)) * 0.06 if len(vals) > 1 else 0)
    ax.scatter(xs, vals, s=7, color=color, edgecolors="white", linewidths=0.3, zorder=5)


def _panel_cav_dominance(ax, letter):
    """Cav1.2/CACNA1C is the dominant L-type, constantly expressed across stages
    (~99% of L-type, ~97-99% of CMs), while the specialized anti-Ro conduction
    targets (Cav1.3, Cav3.1) are minor and decline. Sim2021 within-study; Smart-seq."""
    dom = read_json("ca_subtype_dominance_dev")
    pg = dom["per_gene"]; byd = dom.get("per_gene_by_donor", {})
    genes = ["CACNA1C", "CACNA1D", "CACNA1G"]
    labels = {g: cfg.channel_label(g, short=True) for g in genes}
    cols = {"CACNA1C": ACCENT, "CACNA1D": FETAL, "CACNA1G": "#7a4fa0"}
    x = np.arange(len(_DOM_STAGES)); w = 0.26
    for i, g in enumerate(genes):
        vals = [pg[g][s]["mean_tp10k"] for s in _DOM_STAGES]
        bx = x + (i - 1) * w
        ax.bar(bx, vals, w, label=labels[g], color=cols[g], edgecolor="black", lw=0.4)
        for xpos, s in zip(bx, _DOM_STAGES):       # donor-level dispersion (3 Sim donors/stage)
            _donor_points(ax, xpos, byd.get(g, {}).get(s, {}))
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(_DOM_STAGES)
    ax.set_ylabel("mean TP10K (log scale)")
    ax.set_title(f"{letter} · Cav1.2 dominant & constant; specialists decline (whole-CM)", fontsize=8.5)
    ax.legend(fontsize=6.5, frameon=False, title="L/T-type", title_fontsize=6.5, loc="lower left")
    # annotate Cav1.2 %-expressing constancy (the 'workhorse does not switch' point)
    pe = [pg["CACNA1C"][s]["pct_expr"] for s in _DOM_STAGES]
    ax.text(0.98, 0.97, f"Cav1.2 in {pe[0]:.0f}/{pe[1]:.0f}/{pe[2]:.0f}% of CMs\n(fetal/child/adult)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6, color=ACCENT, style="italic")


def _panel_inward_outward(ax, letter):
    """The CONSEQUENCE of the hand-off: the cell's preferred transcriptional BALANCE between the
    reactivatable inward channels (bulk L-type ICaL Cav1.2 + late-Na SCN5A) and its outward
    repolarizing-K reserve, across fetal->child->adult. The fetal cell is provisioned toward a
    larger PROPORTION of L-Ca/Na channels that could be potentiated under a hERG block, against a
    thinner K reserve; the balance falls with maturation. This is a channel-composition ratio
    (provision/preference), NOT a current/flux magnitude and NOT a measurement of afterdepolarizations."""
    krd = read_json("k_repolarization_dev")
    ps = krd["per_stage"]; byd = krd.get("inward_over_outward_by_donor", {})
    vals = [ps[s]["EAD_substrate"]["inward_over_outward_Cav12only"] for s in _DOM_STAGES]
    cols = [FETAL, "#9a6a6a", ADULT]
    ax.bar(range(3), vals, color=cols, edgecolor="black", lw=0.4)
    for i, s in enumerate(_DOM_STAGES):            # donor-level dispersion (3 Sim donors/stage)
        _donor_points(ax, i, byd.get(s, {}))
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.03, f"{v:.1f}", ha="center", fontsize=7.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(_DOM_STAGES)
    ax.set_ylabel("inward(L-Ca+Na) : outward(K)\nchannel-transcript balance")
    allv = vals + [v for s in _DOM_STAGES for v in byd.get(s, {}).values()]
    ax.set_ylim(0, max(allv) * 1.2)
    ax.set_title(f"{letter} · Reactivatable-inward : repolarizing balance\n"
                 "(fetal provisioned toward potentiatable L-Ca/Na)", fontsize=8)


def fig1():
    """The nodal channel program, AVN-focused (the CHB block site). A (full-width schematic):
    which antibody arrhythmic phenotype targets which cell type / channels (cited encoding).
    B: composition of the four cardiac Cav subtypes (Cav1.2/1.3/3.1/3.2) by node cell type and
    stage. STAGE axis — C: Sim within-study fetal/adult ratio per channel (the Ca-block targets are
    strongly fetal-elevated while the K/QT targets sit near 1 — the developmental enrichment is
    Ca-specific). CELL-TYPE axis — D fetal AVN, E adult AVN pacemaker-vs-working (same
    channel set / order): the anti-Ro Ca targets are pacemaker-enriched in BOTH D and E (the
    constitutive finding; the formal cell-type × stage interaction is in Methods +
    nodal_specificity_interaction). SAN equivalents are in the supplement. One consistent
    phenotype encoding (colour=stage, hatch=arrhythmia type, bold=CHB target) throughout."""
    fig = plt.figure(figsize=(11, 10.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.72, 1.0, 1.0], hspace=0.6, wspace=0.28)
    # row 0: organizing schematic (full width)
    _panel_schematic(fig.add_subplot(gs[0, :]), "A")
    # row 1: B = Cav-subtype composition + C = STAGE axis (Sim fetal/adult ratio), single row
    _panel_cav_subtype_stack(fig.add_subplot(gs[1, 0]), "B")
    axC = fig.add_subplot(gs[1, 1]); _panel_fetal_ratio(axC, "C")
    axC.set_ylabel("Sim fetal / adult CM (STAGE)", fontsize=7.5)
    # row 2: CELL-TYPE axis — pacemaker vs working in the AVN, fetal (D) and adult (E)
    axD = fig.add_subplot(gs[2, 0]); axE = fig.add_subplot(gs[2, 1])
    _panel_pace_work(axD, "AVN", "D", stage="fetal")
    _panel_pace_work(axE, "AVN", "E", stage="adult")
    axD.set_ylabel("← CELL-TYPE axis (AVN, fetal & adult) →", fontsize=7, color="#666", style="italic")
    fig.suptitle("F1 · The nodal channel program (AVN), by antibody phenotype, cell type, and "
                 "developmental stage",
                 fontsize=10.5, fontweight="bold", y=0.995)
    _save(fig, "F1_fetal_nodal_channel_program")
    plt.close(fig)


# ---------------------------------------------------------------- F2 accessibility
def fig2():
    """Why the FETUS, not the adult. Given that the anti-Ro Ca targets are pacemaker-enriched
    at both ages (the constitutive finding is now shown in F1B,C), the fetal selectivity is one
    of accessibility: the fetal node exposes far more of its Ca machinery at the antibody-
    reachable surface (A), running a surface-Ca-entry mode (B) while the adult hides the same
    channels behind a mature SR. Both panels are INTERNAL per-cell ratios (surface/total; NCX/
    RYR2), so they cancel uniform depth/library-size differences and are robust across the two
    datasets; bars are blocked by dataset (Lim fetal | Kanemaru adult) so the load-bearing
    comparisons are the WITHIN-dataset, cross-cell-type contrasts. The absolute SR-module
    ensemble (an absolute cross-dataset sum, not depth-robust) is moved to the supplement
    (S6); its content is already carried here by RYR2 inside the NCX:RYR2 mode ratio (B).
    (C) the within-cell ion-channel CLASS budget (surface Ca / SR Ca / K / Na / HCN /
    Na-K-ATPase) shifts fetal->adult in BOTH the AVN and the working myocardium, read
    separately: the surface-Ca
    share falls and the SR-Ca share rises with maturation (an internal composition, depth-
    robust; the fetal-vs-adult shift is cross-dataset, direction-only)."""
    # Two rows: A,B side-by-side on top; C as a full-width HORIZONTAL stacked-bar panel
    # on the bottom (the 6-class budget reads better horizontally — thin K/Na/HCN/
    # transporter slices get length, and the 6-entry legend has room).
    fig = plt.figure(figsize=(9.6, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.0], hspace=0.5, wspace=0.30)
    _panel_surface_frac(fig.add_subplot(gs[0, 0]), "A")
    _panel_ncx_mode(fig.add_subplot(gs[0, 1]), "B")
    _panel_class_budget(fig.add_subplot(gs[1, :]), "C", horizontal=True)
    fig.suptitle("F2 · The Ca targets are nodal at both ages (F1) — fetal vulnerability comes "
                 "from accessibility, not nodal-specificity",
                 fontsize=8.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, "F2_accessibility")
    plt.close(fig)


# ---------------------------------------------------------------- F3 developmental Ca-reliance hand-off (Smart-seq)
def fig3():
    """F3: developmental axis (Sim 2021 within-study, whole-CM, fetal->child->adult).
    (A) surface Cav subtypes — Cav1.2 constant, Cav1.3/Cav3.1 decline; (B) RYR2 (SR) rises;
    (C) inward:outward composition ratio (reactivatable inward L-type ICaL + late-Na over the
    repolarizing-K reserve) falls with maturation. Panel C is a channel-composition ratio, not a
    current/flux magnitude. Whole-CM only (the Sim atlas has no dissected node); the node-specific
    evidence is Fig 1."""
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
    _panel_cav_dominance(axes[0], "A")
    _panel_sr_traj(axes[1], "B")
    _panel_inward_outward(axes[2], "C")
    fig.suptitle("F3 · Developmental hand-off of calcium reliance: surface specialists → SR, and a "
                 "fetal arrhythmia-substrate balance (whole-CM, Sim within-study)",
                 fontsize=9.5, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, "F3_developmental_backdrop")
    plt.close(fig)


# ---------------------------------------------------------------- S7 developmental replication / detail
def fig_san():
    """Supplementary: SAN pacemaker-vs-working — the sinoatrial-node equivalents of the AVN
    panels in F1 (B,C). The main figure focuses on the AVN because it is the CHB block site;
    the SAN contrasts are shown here for completeness (same encoding: bars = enrichment
    direction, labels = phenotype stage colour + CHB bold)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    _panel_pace_work(axes[0], "SAN", "A", stage="fetal")
    _panel_pace_work(axes[1], "SAN", "B", stage="adult")
    fig.suptitle("S3 · Sinoatrial node (SAN) pacemaker vs working CM — fetal (A) and adult (B); "
                 "AVN equivalents are in Fig 1B,C",
                 fontsize=9, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "S3_SAN_pacemaker_vs_working")
    plt.close(fig)


# ---------------------------------------------------------------- S3 dataset reliability / consistency
def fig_reliability():
    """Reliability of the Smart-seq developmental results: (A) non-10x Cui/Wang STRT
    data corroborate the 10x Sim within-study magnitudes; (B) the fetal/adult direction
    agrees between the within-study 10x ratio and the cross-platform Cui/Wang ratio;
    (C) the developmental trend holds across 13 independent Lázár donors; (D) the cell
    counts behind every contrast. All on the current Smart-seq pipeline."""
    v = read_json("fetal_vs_adult_vulnerability")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))

    # --- A: multi-platform agreement (10x Sim vs non-10x Cui fetal / Wang adult) ---
    ax = axes[0, 0]
    tg = v["targets"]
    genes = [g for g in ["CACNA1C", "CACNA1D", "CACNA1G", "HCN4", "KCNH2", "KCNQ1", "SCN5A"]
             if g in tg]
    # fetal: Sim fetal vs Cui late_fetal; adult: Sim adult vs Wang adult. Plot log10(x+eps).
    eps = 1e-3
    sim_f = [np.log10(tg[g]["sim"].get("fetal", np.nan) + eps) for g in genes]
    cui_f = [np.log10(tg[g]["cui"].get("late_fetal", np.nan) + eps) for g in genes]
    sim_a = [np.log10(tg[g]["sim"].get("adult", np.nan) + eps) for g in genes]
    wang_a = [np.log10(tg[g].get("wang_adult", np.nan) + eps) for g in genes]
    ax.scatter(sim_f, cui_f, color=FETAL, s=34, edgecolor="black", lw=0.4, label="fetal: Sim(10x) vs Cui(STRT)", zorder=3)
    ax.scatter(sim_a, wang_a, color=ADULT, s=34, edgecolor="black", lw=0.4, label="adult: Sim(10x) vs Wang(STRT)", zorder=3)
    lims = [min(sim_f + cui_f + sim_a + wang_a) - 0.3, max(sim_f + cui_f + sim_a + wang_a) + 0.3]
    ax.plot(lims, lims, color="#999", lw=0.8, ls="--")
    # Spearman of the platform agreement (direction), pooled
    from scipy.stats import spearmanr
    xs = sim_f + sim_a; ys = cui_f + wang_a
    ok = [(a, b) for a, b in zip(xs, ys) if np.isfinite(a) and np.isfinite(b)]
    rho, p = spearmanr([a for a, _ in ok], [b for _, b in ok]) if len(ok) > 2 else (np.nan, np.nan)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("10x Sim  log10 mean TP10K"); ax.set_ylabel("non-10x (Cui/Wang) log10 mean")
    ax.set_title(f"A · Non-10x STRT corroborates 10x (ρ={rho:.2f}, p={p:.1e})", fontsize=8.5)
    ax.legend(fontsize=6, frameon=False, loc="upper left")

    # --- B: cross-platform direction agreement (Smart-seq-native robustness): the
    # within-study 10x fetal/adult ratio vs the cross-platform (Cui-fetal / Wang-adult)
    # ratio. Same fetal-elevated vs adult-flat verdict on both = the result is not a
    # single-platform artifact. ---
    ax = axes[0, 1]
    genes_b = [g for g in ["CACNA1G", "CACNA1D", "HCN4", "CACNA1C", "KCNH2", "KCNQ1", "SCN5A"]
               if g in tg]
    sim_r, xplat_r = [], []
    for g in genes_b:
        sim_r.append(tg[g].get("sim_fetal_over_adult", np.nan))
        cf = tg[g]["cui"].get("late_fetal", np.nan); wa = tg[g].get("wang_adult", np.nan)
        xplat_r.append((cf + eps) / (wa + eps) if np.isfinite(cf) and np.isfinite(wa) else np.nan)
    x = np.arange(len(genes_b)); w = 0.4
    ax.bar(x - w/2, np.log2(sim_r), w, label="Sim 10x (within-study)", color=ADULT, edgecolor="black", lw=0.4)
    ax.bar(x + w/2, np.log2(xplat_r), w, label="Cui/Wang STRT (cross-platform)", color=WORK, edgecolor="black", lw=0.4)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(genes_b, rotation=30, ha="right")
    ax.set_ylabel("log2 fetal / adult ratio")
    ax.set_title("B · Fetal/adult direction agrees across platforms", fontsize=8.5)
    ax.legend(fontsize=6, frameon=False)

    # --- C: per-donor spread (Lázár n=13), cell-count-sized points ---
    ax = axes[1, 0]
    lz = read("lazar_dev_replication")
    sizes = 8 + 30 * (lz["n_cm"] / lz["n_cm"].max())
    ax.scatter(lz["age_pcw"], lz["RYR2"], s=sizes, color=ADULT, edgecolor="black", lw=0.4, alpha=0.85, zorder=3)
    z = np.polyfit(lz["age_pcw"], lz["RYR2"], 1)
    xs = np.linspace(lz["age_pcw"].min(), lz["age_pcw"].max(), 20)
    ax.plot(xs, np.polyval(z, xs), color=ACCENT, lw=1.5)
    ax.set_xlabel("gestational age (pcw)"); ax.set_ylabel("RYR2 (linear, whole-CM)")
    ax.set_title("C · Trend holds across 13 donors (ρ=+0.63, p=0.022)\npoint size = cells/donor", fontsize=8)

    # --- D: cell counts behind every contrast ---
    ax = axes[1, 1]
    nc = dict(v["n_cells"])
    # add dissected-nodal pacemaker counts from the within-fetal table
    fpw = read("fetal_pacemaker_vs_working")
    for node in ("SAN", "AVN"):
        sub = fpw[fpw.node == node]
        if len(sub):
            nc[f"lim_{node}_pace"] = int(sub["n_pacemaker"].iloc[0])
    order = ["sim_fetal", "sim_child", "sim_adult", "cui_early_fetal", "cui_late_fetal",
             "wang_adult", "lim_SAN_pace", "lim_AVN_pace"]
    order = [k for k in order if k in nc]
    vals = [nc[k] for k in order]
    barcols = [FETAL if "fetal" in k or "pace" in k else ADULT if "adult" in k else WORK for k in order]
    yy = np.arange(len(order))
    ax.barh(yy, vals, color=barcols, edgecolor="black", lw=0.4)
    for yi, val in zip(yy, vals):
        ax.text(val, yi, f" {val:,}", va="center", fontsize=6.5)
    ax.set_yticks(yy); ax.set_yticklabels([k.replace("_", " ") for k in order], fontsize=7)
    ax.set_xscale("log"); ax.invert_yaxis()
    ax.set_xlabel("cells per group (log scale)")
    ax.set_title("D · Cells behind each contrast (nodal = n=1 donor)", fontsize=8.5)

    fig.suptitle("S2 · Dataset reliability: the Smart-seq developmental results are "
                 "platform-corroborated, direction-consistent, and donor-replicated",
                 fontsize=10, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _save(fig, "S2_dataset_reliability")
    plt.close(fig)


# ---------------------------------------------------------------- S4 platform-suitability QC
def _panel_complexity(ax, letter):
    q = read_json("platform_qc")
    comp = q["A_per_cell_complexity"]
    items = [(k, v) for k, v in comp.items() if "median_genes_per_cell" in v]
    labels = [k.replace("_", "\n") for k, _ in items]
    vals = [v["median_genes_per_cell"] for _, v in items]
    plats = [v.get("platform", "") for _, v in items]
    # colour: droplet 10x (snRNA + Multiome) vs STRT/Smart-seq; Multiome hatched so it is not
    # read as identical to 3'-snRNA (it is the platform carrying the adult NODAL result).
    cols = [ADULT if "10x" in p else ACCENT for p in plats]
    hatches = ["////" if "Multiome" in p else "" for p in plats]
    x = np.arange(len(items))
    for xi, val, c, h in zip(x, vals, cols, hatches):
        ax.bar(xi, val, color=c, edgecolor="black", lw=0.4, hatch=h)
        ax.text(xi, val + 40, f"{val:.0f}", ha="center", fontsize=6.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("median genes per cell")
    ax.set_title(f"{letter} · Smart-seq recovers ≥ droplet 10x complexity", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=ADULT),
               plt.Rectangle((0, 0), 1, 1, color=ADULT, hatch="////"),
               plt.Rectangle((0, 0), 1, 1, color=ACCENT)]
    ax.legend(handles, ["10x snRNA", "10x Multiome\n(adult nodal)", "STRT (Smart-seq)"],
              fontsize=5.5, frameon=False)


def _panel_depth_robust(ax, letter):
    q = read_json("platform_qc")
    nodal = q["B_nodal_depth_robustness"]
    genes = ["CACNA1D", "CACNA1G", "HCN4", "KCNH2", "SCN5A"]
    avn = nodal["fetal_AVN"]
    full = [avn["log2fc_full"].get(g, np.nan) for g in genes]
    matched = [avn["log2fc_depth_matched"].get(g, np.nan) for g in genes]
    x = np.arange(len(genes)); w = 0.38
    ax.bar(x - w/2, full, w, label="full data", color=FETAL, edgecolor="black", lw=0.4)
    ax.bar(x + w/2, matched, w, label="depth-matched", color=WORK, edgecolor="black", lw=0.4)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(genes, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("pacemaker/working log2 FC (AVN)")
    ax.set_title(f"{letter} · Nodal contrast survives depth-matching", fontsize=8)
    ax.legend(fontsize=6, frameon=False)


_GS_SETS = ["electrophysiology", "structural", "housekeeping"]
_GS_LAB = {"electrophysiology": "electro-\nphysiology", "structural": "structural",
           "housekeeping": "house-\nkeeping"}


def _panel_geneset_capture(ax, letter):
    q = read_json("geneset_capture_qc")
    x = np.arange(len(_GS_SETS)); w = 0.2
    bar_specs = []
    for stage, col in [("fetal", FETAL), ("adult", ADULT)]:
        groups = q["stages"][stage]["groups"]
        sm = q["stages"][stage]["smartseq_ref"]
        sn = [g for g in groups if groups[g]["platform"] == "snRNA"]
        sn_pct = [np.mean([groups[g]["genesets"][s]["pct_detected"] for g in sn]) for s in _GS_SETS]
        sm_pct = [groups[sm]["genesets"][s]["pct_detected"] for s in _GS_SETS]
        bar_specs.append((f"{stage} snRNA", sn_pct, col, 0.55))
        bar_specs.append((f"{stage} Smart-seq", sm_pct, col, 1.0))
    for i, (lab, vals, col, alpha) in enumerate(bar_specs):
        ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=col, alpha=alpha, edgecolor="black", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels([_GS_LAB[s] for s in _GS_SETS], fontsize=7)
    ax.set_ylabel("% of gene set detected"); ax.set_ylim(0, 118)
    ax.set_title(f"{letter} · Gene-set capture: snRNA vs Smart-seq", fontsize=8)
    ax.legend(fontsize=5.5, frameon=False, ncol=2, loc="upper left")


def _panel_pattern_corr(ax, letter):
    q = read_json("geneset_capture_qc")
    rows = []
    for stage in ("fetal", "adult"):
        for tag, sets in q["stages"][stage]["snrna_vs_smartseq_pattern_rho"].items():
            for s in _GS_SETS:
                r = sets[s]["rho"]
                if r is not None:
                    rows.append((s, r))
    setvals = {s: [r for (ss, r) in rows if ss == s] for s in _GS_SETS}
    means = [np.mean(setvals[s]) for s in _GS_SETS]
    x = np.arange(len(_GS_SETS))
    ax.bar(x, means, color=[ACCENT, "#6a8cbf", WORK], edgecolor="black", lw=0.4)
    for s_i, s in enumerate(_GS_SETS):
        ax.scatter([s_i] * len(setvals[s]), setvals[s], color="black", s=9, zorder=3)
    for xi, m in zip(x, means):
        ax.text(xi, m + 0.03, f"{m:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([_GS_LAB[s] for s in _GS_SETS], fontsize=7)
    ax.set_ylabel("snRNA vs Smart-seq ρ (per-gene)"); ax.set_ylim(0, 1)
    ax.set_title(f"{letter} · Expression pattern preserved", fontsize=8)


import config_channel_subtypes as _cst


def _panel_channel_capture(ax, letter):
    """Per-channel ABSOLUTE detection % across platforms for every cardiac-expressed channel
    subtype, STAGE-MATCHED: fetal working CM (snRNA Lim/Protze/Sim vs Smart-seq Cui) and adult
    working CM (snRNA Sim/Kanemaru vs Smart-seq Wang). Each channel speaks for itself — unlike a
    within-family ratio, this is not distorted by the abundant L-type denominator. snRNA detects
    most cardiac channels comparably to or better than Smart-seq; for CACNA1G this holds at the
    fetal stage (the basis of the study's enrichment claim), while at the adult stage CACNA1G is
    near its detection floor on every platform and Smart-seq (Wang) detects it marginally higher."""
    q = read_json("geneset_capture_qc")

    def stage_vals(stage):
        cc = q["stages"][stage]["cardiac_channel_capture"]
        sm = q["stages"][stage]["smartseq_ref"]
        sn_groups = [g for g in cc if g != sm]
        return cc, sm, sn_groups

    fcc, fsm, fsn = stage_vals("fetal")
    acc_, asm, asn = stage_vals("adult")
    # order channels by ion family then fetal snRNA detection
    chans = sorted(_cst.CARDIAC_SUBTYPE_GENES,
                   key=lambda g: (_cst.ion(g), -np.mean([fcc[t][g]["detect_pct"] for t in fsn])))
    f_sn = [np.mean([fcc[t][g]["detect_pct"] for t in fsn]) for g in chans]
    f_ss = [fcc[fsm][g]["detect_pct"] for g in chans]
    a_sn = [np.mean([acc_[t][g]["detect_pct"] for t in asn]) for g in chans]
    a_ss = [acc_[asm][g]["detect_pct"] for g in chans]

    y = np.arange(len(chans)); h = 0.18; gap = 0.07   # gap separates the fetal pair from the adult pair
    ax.barh(y - 1.5*h - gap, f_sn, h, label="fetal snRNA", color=FETAL, edgecolor="black", lw=0.25)
    ax.barh(y - 0.5*h - gap, f_ss, h, label="fetal Smart-seq", color=FETAL, alpha=0.45, edgecolor="black", lw=0.25)
    ax.barh(y + 0.5*h + gap, a_sn, h, label="adult snRNA", color=ADULT, edgecolor="black", lw=0.25)
    ax.barh(y + 1.5*h + gap, a_ss, h, label="adult Smart-seq", color=ADULT, alpha=0.45, edgecolor="black", lw=0.25)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g} ({_cst.subtype(g).replace('Kv_','').replace('Kir_','Kir-')})" for g in chans],
                       fontsize=4.6)
    ax.invert_yaxis()
    ax.set_xlabel("% cells detecting channel (working CM)")
    ax.set_title(f"{letter} · Per-channel capture, all cardiac subtypes (stage-matched)", fontsize=7.5)
    ax.legend(fontsize=5, frameon=False, loc="lower right", ncol=2)
    for yi, g in zip(y, chans):
        if g in ("CACNA1G", "CACNA1D", "CACNA1C"):
            ax.text(-1, yi, "◄", fontsize=5.5, color="#a02020", va="center", ha="right")


def _panel_nodal_enrichment(ax, letter):
    """Within-snRNA nodal pacemaker-vs-working detection % of the F1D channels in the dissected
    fetal AVN (the CHB lesion site): the conduction-Ca targets are well-detected AND
    pacemaker-enriched, so the headline contrasts are within-platform, not a capture floor."""
    q = read_json("geneset_capture_qc")
    avn = q["nodal_pacemaker_vs_working"]["fetal_AVN"]
    genes = ["CACNA1D", "CACNA1G", "HCN4", "CACNA1C", "SCN5A"]
    pace = [avn["pacemaker_CM"]["channels"][g]["detect_pct"] for g in genes]
    work = [avn["working_CM"]["channels"][g]["detect_pct"] for g in genes]
    y = np.arange(len(genes)); h = 0.38
    ax.barh(y - h/2, pace, h, label="pacemaker CM", color=FETAL, edgecolor="black", lw=0.4)
    ax.barh(y + h/2, work, h, label="working CM", color=WORK, edgecolor="black", lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{g}\n({_cst.subtype(g).replace('Nav_','').replace('Kv_','')})" for g in genes],
                       fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("% cells detecting (fetal AVN, snRNA)")
    ax.set_title(f"{letter} · Targets detected & pacemaker-enriched (within-platform)", fontsize=8)
    ax.legend(fontsize=6, frameon=False, loc="lower right")


def fig_platform():
    """S1 — foundation method-support figure (6 panels) vs the 'snRNA/10x under-captures
    cardiomyocyte genes' critique. (A) per-cell complexity by platform; (B) the nodal
    within-donor contrast survives depth-matching; (C) gene-set capture (EP/structural/HK)
    snRNA vs Smart-seq in matched working CMs; (D) EP expression pattern preserved across
    platforms; (E) per-channel detection for ALL cardiac subtypes, snRNA vs Smart-seq — each
    channel on its own, not a ratio; (F) the F1D channels are detected and pacemaker-enriched
    within the dissected node (within-platform). Reads platform_qc.json + geneset_capture_qc.json."""
    fig = plt.figure(figsize=(13, 8.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=0.40, wspace=0.42)
    _panel_complexity(fig.add_subplot(gs[0, 0]), "A")
    _panel_depth_robust(fig.add_subplot(gs[0, 1]), "B")
    _panel_geneset_capture(fig.add_subplot(gs[0, 2]), "C")
    _panel_pattern_corr(fig.add_subplot(gs[1, 0]), "D")
    _panel_channel_capture(fig.add_subplot(gs[1, 1]), "E")
    _panel_nodal_enrichment(fig.add_subplot(gs[1, 2]), "F")
    fig.suptitle("S1 · Method-support foundation: snRNA captures the cardiac channels the study "
                 "uses comparably to Smart-seq; the headline target contrasts are within-platform "
                 "(bounded checks — no Smart-seq conduction-tissue data exists)",
                 fontsize=9, fontweight="bold", y=0.995)
    _save(fig, "S1_platform_qc")
    plt.close(fig)


def main():
    # F1-F3 + S1-S3 here; F4 (node-channel dot-plot) is built by make_node_handling_figure.py.
    fig1(); fig2(); fig3()
    fig_san(); fig_platform(); fig_reliability()
    pngs = sorted(p.name for p in FIG.glob("*.png"))
    print("wrote", len(pngs), "figures (PNG+SVG) to", FIG)
    for p in pngs:
        print("  ", p)


if __name__ == "__main__":
    main()
