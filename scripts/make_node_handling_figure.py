#!/usr/bin/env python3
"""
make_node_handling_figure.py — figure for the node-resolved fetal-vs-adult channel +
Ca-handling expression analysis (reads node_channel_handling_expression.csv).

seed: 42

Three panels, sharing the gene rows (grouped by family):
  A  FETAL block (snRNA): SAN/AVN x pacemaker/working cell-type gradient
  B  ADULT block (Multiome nuclei): SAN_P / AVN_P / vCM gradient
  C  SIM ladder (snRNA whole-CM): fetal -> child -> adult dev gradient (the only
     cross-dev-stage comparison)
Dot SIZE = detection rate (fraction of cells expressing); dot COLOUR = log2 fold-change of
that group's mean vs the gene's SUB-BATCH mean (per gene, so each gene's cell-type/stage
pattern is visible regardless of absolute level). The reference mean is computed within each
SUB-BATCH (fetal SAN and fetal AVN separately — different samples/timelines; adult Kanemaru as
one; Sim as one), so a colour contrast is comparable WITHIN a block only and never across the
cross-dataset/cross-platform boundaries between blocks.

Output: results/figures/F4_node_channel_handling.{png,svg}  (this dot-plot IS Figure 4)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg
import config_chb_handling_panel as panel

FIG = cfg.PUBLIC / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight",
                     "axes.titlesize": 9, "axes.titleweight": "bold", "svg.fonttype": "none"})

CMAP = "RdBu_r"      # diverging: high z = red (high for that gene), low = blue
GROUP_LABELS = {     # short x-tick labels per group
    "fetal_SAN_pace": "SAN\npace", "fetal_SAN_work": "SAN\nwork",
    "fetal_AVN_pace": "AVN\npace", "fetal_AVN_work": "AVN\nwork",
    "adult_SAN_pace": "SAN\npace", "adult_SAN_work": "SAN\nwork",
    "adult_AVN_pace": "AVN\npace", "adult_AVN_work": "AVN\nwork", "adult_vCM": "vCM\n(bulk)",
    "sim_fetal_CM": "fetal", "sim_child_CM": "child", "sim_adult_CM": "adult",
}
BLOCK_TITLE = {
    "fetal (snRNA)":         "A · Fetal node\n(Lim SAN | Protze AVN, snRNA)",
    "adult (Multiome nuc.)": "B · Adult node\n(Kanemaru, Multiome)",
    "Sim ladder (snRNA)":    "C · Sim dev ladder\n(whole-CM, snRNA)",
}


def _dot_size(detect):
    """detection fraction (0-1) -> marker area."""
    return 12 + detect * 240


def main():
    df = pd.read_csv(cfg.OUT_TABLES / "node_channel_handling_expression.csv", comment="#")
    # gene row order = panel order, grouped by family; keep only genes present
    genes = [g for g in panel.PANEL_GENES if g in set(df.gene)]
    fam_of = {g: panel.FAMILY[g] for g in genes}
    # canonical labels with the ionic current annotation (protein / current (gene)),
    # consistent with the F1/F3 base label but keeping the current context useful in
    # this detailed heatmap; falls back to the panel label for any unmapped gene.
    ylabels = [cfg.channel_label_current(g) if g in cfg.CHANNEL_LABEL else panel.LABEL[g]
               for g in genes]
    y_of = {g: i for i, g in enumerate(genes[::-1])}   # top gene at top

    blocks = list(panel.STAGE_BLOCKS.items())
    width_ratios = [len(gs) for _, gs in blocks]
    fig, axes = plt.subplots(1, len(blocks), figsize=(10.5, 0.30 * len(genes) + 1.8),
                             gridspec_kw={"width_ratios": width_ratios, "wspace": 0.12},
                             sharey=True)

    mean_by = df.set_index(["group", "gene"])["mean_tp10k"].to_dict()
    det_by = df.set_index(["group", "gene"])["detect_frac"].to_dict()

    present = set(df.group)
    for ax, (block, groups) in zip(axes, blocks):
        groups = [g for g in groups if g in present]
        # COLOUR = per-gene log2 fold-change of each group's mean vs the geometric mean of
        # that gene ACROSS the sub-batch's groups. This is computed WITHIN each sample/batch
        # sub-group (fetal SAN and AVN are separate batches -> normalised separately) so a
        # colour contrast is never a batch effect. Unlike a z-score, log2FC preserves the
        # MAGNITUDE of the contrast even when a sub-batch has only 2 groups (a z-score of 2
        # values is always +/-1, which made the old panel look binary). zmat: (grp,gene)->log2FC.
        subbatches = [[g for g in sb if g in present]
                      for sb in panel.Z_SUBGROUPS[block]]
        zmat = {}                                  # (grp, gene) -> log2FC vs sub-batch geomean
        for sb in subbatches:
            for g in genes:
                vals = np.array([mean_by.get((grp, g), np.nan) for grp in sb], float)
                lv = np.log2(vals + cfg.EPS)
                ref = np.nanmean(lv)               # log2 geometric mean across the sub-batch
                for grp, l in zip(sb, lv):
                    zmat[(grp, g)] = l - ref
        zall = np.array([zmat[(grp, g)] for grp in groups for g in genes], float)
        vmax = np.nanpercentile(np.abs(zall), 98) or 1.0

        for xi, grp in enumerate(groups):
            for g in genes:
                d = det_by.get((grp, g), 0.0)
                ax.scatter(xi, y_of[g], s=_dot_size(d), c=[zmat[(grp, g)]],
                           cmap=CMAP, vmin=-vmax, vmax=vmax,
                           edgecolors="black", linewidths=0.3, zorder=3)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([GROUP_LABELS.get(grp, grp) for grp in groups], fontsize=7)
        ax.set_xlim(-0.6, len(groups) - 0.4)
        ax.set_ylim(-0.6, len(genes) - 0.4)
        # vertical divider at each sample/batch boundary (e.g. SAN | AVN in the fetal panel,
        # where SAN=Lim and AVN=Protze are separate samples normalised separately). Colours are
        # comparable only WITHIN a sub-batch, so the divider marks where cross-reading stops.
        boundary = 0
        for sb in subbatches[:-1]:
            boundary += len(sb)
            ax.axvline(boundary - 0.5, color="#999", lw=0.8, ls=(0, (2, 2)), zorder=1)
        ax.set_title(BLOCK_TITLE.get(block, block), fontsize=8.5)
        ax.tick_params(length=0)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_visible(False)

    # y labels + family separators on the leftmost axis
    ax0 = axes[0]
    ax0.set_yticks(range(len(genes)))
    ax0.set_yticklabels(ylabels[::-1], fontsize=7)
    # family group bars (light separators between families)
    fams_top_down = [fam_of[g] for g in genes]
    for i in range(1, len(genes)):
        if fams_top_down[i] != fams_top_down[i - 1]:
            yline = len(genes) - i - 0.5
            for ax in axes:
                ax.axhline(yline, color="0.85", lw=0.6, zorder=0)
    # family labels at right margin of last axis
    axL = axes[-1]
    seen = set()
    for g in genes:
        fam = fam_of[g]
        if fam in seen:
            continue
        seen.add(fam)
        idxs = [y_of[gg] for gg in genes if fam_of[gg] == fam]
        axL.text(1.02, np.mean(idxs), fam, transform=axL.get_yaxis_transform(),
                 fontsize=6.5, rotation=270, va="center", ha="left", color="0.4")

    # colourbar + size legend
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-1, 1))
    cb = fig.colorbar(sm, ax=axes, fraction=0.018, pad=0.06)
    cb.set_label("log2 fold-change vs gene's\nsub-batch mean (red=high).\nComparable WITHIN a block only\n(datasets differ across blocks)", fontsize=6.0)
    cb.ax.tick_params(labelsize=6)
    size_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="0.5",
                           markeredgecolor="black", markeredgewidth=0.3,
                           markersize=np.sqrt(_dot_size(d)), label=f"{int(d*100)}%")
                    for d in (0.1, 0.5, 0.9)]
    axes[0].legend(handles=size_handles, title="% cells\nexpressing", fontsize=6,
                   title_fontsize=6, loc="upper left", bbox_to_anchor=(0.0, -0.04),
                   ncol=3, frameon=False, handletextpad=0.1, columnspacing=0.8)

    fig.suptitle("Cardiac channels & Ca-handling: cell-type and developmental gradients\n"
                 "(anti-Ro Ca targets + CIRCEP.115.003432 handling apparatus)",
                 fontsize=9.5, y=1.02)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"F4_node_channel_handling.{ext}")
    print("wrote", FIG / "F4_node_channel_handling.png", "(+ svg)")
    # retire the pre-rename output (this dot-plot is now Figure 4)
    for ext in ("png", "svg"):
        old = FIG / f"F_node_channel_handling.{ext}"
        if old.exists():
            old.unlink()


if __name__ == "__main__":
    main()
