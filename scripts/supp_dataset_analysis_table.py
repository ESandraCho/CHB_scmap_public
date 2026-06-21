#!/usr/bin/env python3
"""
supp_dataset_analysis_table.py — reviewer-facing supplementary tables that make explicit
which datasets were used and which analysis used which dataset(s). Several datasets are
combined across the study; this table prevents that from being confusing.

seed: 42

Emits:
  (1) DATASET provenance table — one row per dataset (accession, stage, chemistry,
      platform used, donors, citation), from config_chb_handling_panel.DATASET_PROVENANCE
      (transcribed from DATA_SOURCES.md, the authoritative source).
  (2) ANALYSIS -> DATASETS map — one row per analysis/figure listing the datasets it used
      and the comparison type, from config_chb_handling_panel.ANALYSIS_DATASET_MAP.

No data values are hard-coded here; provenance + mapping are curated config, formatted
for the manuscript supplement.

Outputs:
  results/tables/supp_datasets_provenance.csv
  results/tables/supp_analysis_dataset_map.csv
  results/reports/supp_dataset_analysis_tables.md   (both, formatted for the supplement)
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg
import config_chb_handling_panel as panel

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("supp_tbl")


def df_to_md(df: pd.DataFrame) -> str:
    """Minimal GitHub-flavoured markdown table (avoids the optional `tabulate` dep)."""
    cols = list(df.columns)
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = [
        "| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |"
        for row in df.itertuples(index=False)
    ]
    return "\n".join([head, sep, *body])


def dataset_rows() -> pd.DataFrame:
    rows = []
    for key, p in panel.DATASET_PROVENANCE.items():
        rows.append(dict(
            dataset=p["name"],
            object_or_dir=key,
            accession=p["accession"],
            tissue_stage=p["tissue_stage"],
            chemistry=p["chemistry"],
            platform_used=p["platform_used"],
            donors=p["donors"],
            citation=p["citation"],
        ))
    return pd.DataFrame(rows)


def analysis_rows() -> pd.DataFrame:
    name = {k: v["name"] for k, v in panel.DATASET_PROVENANCE.items()}
    rows = []
    for a in panel.ANALYSIS_DATASET_MAP:
        missing = [d for d in a["datasets"] if d not in name]
        if missing:
            raise KeyError(f"ANALYSIS_DATASET_MAP references unknown datasets: {missing}")
        rows.append(dict(
            analysis=a["analysis"],
            datasets="; ".join(name[d] for d in a["datasets"]),
            cell_types=a["cell_types"],
            method=a["method"],
            comparison_type=a["comparison_type"],
            figures=a["figures"],
        ))
    return pd.DataFrame(rows)


def main():
    ds = dataset_rows()
    am = analysis_rows()

    cfg.OUT_TABLES.mkdir(parents=True, exist_ok=True)
    p_ds = cfg.OUT_TABLES / "supp_datasets_provenance.csv"
    p_am = cfg.OUT_TABLES / "supp_analysis_dataset_map.csv"
    with open(p_ds, "w") as f:
        f.write(cfg.header("supplementary dataset provenance (one row per dataset); "
                           "from config_chb_handling_panel.DATASET_PROVENANCE") + "\n")
        ds.to_csv(f, index=False)
    with open(p_am, "w") as f:
        f.write(cfg.header("supplementary analysis->datasets map (one row per analysis/"
                           "figure); from config_chb_handling_panel.ANALYSIS_DATASET_MAP") + "\n")
        am.to_csv(f, index=False)

    cfg.OUT_REPORTS.mkdir(parents=True, exist_ok=True)
    rp = cfg.OUT_REPORTS / "supp_dataset_analysis_tables.md"
    with open(rp, "w") as f:
        f.write(cfg.header("supplementary dataset + analysis-mapping tables") + "\n\n")
        f.write("## Supplementary Table 1. Single-cell/-nucleus datasets used\n\n")
        f.write(df_to_md(ds) + "\n\n")
        f.write("## Supplementary Table — analysis-to-dataset map "
                "(which analysis used which datasets)\n\n")
        f.write(df_to_md(am) + "\n")

    logger.info("wrote %s (%d datasets), %s (%d analyses), %s",
                p_ds.name, len(ds), p_am.name, len(am), rp.name)


if __name__ == "__main__":
    main()
