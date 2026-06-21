# CHB single-cell channel-expression analysis — reproducibility repository

Public code and computed results for the manuscript *"Developmental and cell-type
channel-expression programs as a candidate substrate for the fetal selectivity of
autoimmune congenital heart block."*

Single-cell/-nucleus RNA-seq analysis of cardiac ion-channel expression across cell
type (nodal pacemaker vs working cardiomyocyte) and developmental stage (fetal → adult),
focused on the calcium-channel program relevant to anti-Ro/SSA congenital heart block.

All results are reproducible from public data by the scripts here, in a pinned
environment. **No values are hard-coded** — every number is computed from the data at
run time. No data are bundled; datasets are fetched from their public accessions.

---

## What is in this repository

```
scripts/        analysis + config + fetch scripts (29 .py)
results/
  tables/       computed result tables (CSV/JSON/TSV) + CHB_supplementary_tables.xlsx
  figures/      rendered figures (PNG + SVG); final/ = submission JPG+SVG
requirements.txt
README.md
REPRODUCIBILITY.md   full data → script → table → figure trace (start here to verify a number)
DATA_SOURCES.md      dataset accessions + processing provenance
```

The committed `results/tables/` and `results/figures/` are the exact outputs of the
scripts; you can inspect them directly, or regenerate them from the data (below).

## Quick start: regenerate the figures from the committed tables

This needs **no data download** — it re-renders the figures from the committed result tables:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/make_figures.py              # F1–F3, S1–S3 (PNG + SVG)
python scripts/make_node_handling_figure.py # F4
```

## Full reproduction: from raw data

```bash
# 1. environment
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. fetch data + annotations (one-time, network)
python scripts/fetch_data.py                # fetal-node 10x triplets (Lim SAN/AVN) from GEO
python scripts/fetch_annotation_caches.py   # UniProt subcellular-localization cache
python scripts/fetch_genesets.py            # EBI QuickGO gene sets
#   The adult (Kanemaru), developmental (Sim 2021), Cui/Wang STRT and Lázár objects are
#   obtained from their accessions and QC-processed to the *_qc.h5ad the analysis reads —
#   see DATA_SOURCES.md for accessions and the QC/normalization that produces them.

# 3. classification table that drives the Fig-1 encoding
python scripts/build_channelopathy_classification.py

# 4. analyses (each reads the h5ad/raw objects + configs, writes results/tables/)
python scripts/fetal_pacemaker_vs_working.py
python scripts/adult_pacemaker_vs_working.py
python scripts/nodal_specificity_interaction.py
python scripts/cav_subtype_pace_work_stats.py     # after node_channel_handling_expression.py
python scripts/node_channel_handling_expression.py
python scripts/ca_accessibility_analysis.py
python scripts/sr_activity_proxies.py
python scripts/channel_class_budget.py
python scripts/ca_subtype_dominance_dev.py
python scripts/k_repolarization_dev.py
python scripts/fetal_vs_adult_vulnerability.py
python scripts/lazar_dev_replication.py
python scripts/channel_subtype_conduction.py
python scripts/geneset_capture_qc.py
python scripts/platform_qc.py

# 5. methods / supplementary tables
python scripts/node_datasets_methods_table.py
python scripts/supp_dataset_analysis_table.py
python scripts/build_supplementary_tables.py      # -> CHB_supplementary_tables.xlsx

# 6. figures
python scripts/make_figures.py
python scripts/make_node_handling_figure.py
```

`REPRODUCIBILITY.md` gives the full **figure → script → table** map and the run-order with
dependency notes. Configuration (paths, QC thresholds, gene panels, channel taxonomy) lives in
the `config_*.py` modules; the analysis scripts import them and hard-code nothing.

## Notes for reviewers

- **Determinism:** seed 42 throughout; layout is deterministic. Re-running a script reproduces
  its committed table byte-for-byte (apart from the generation-date stamp in the header comment).
- **Statistics scope:** only the within-fetal pacemaker-vs-working contrasts (per-gene
  Mann–Whitney + BH-FDR) and the Lázár donor-level replication (Spearman) carry inferential
  statistics; all cross-dataset fetal-vs-adult contrasts are direction-only (single fetal donor
  per node). This is stated in the manuscript and in each figure legend.
- **Interpretive `.md` report summaries** are written to `results/reports/`. They are derived from
  the same tables and are not required to reproduce any figure or number.

## License

Released under the MIT License (see `LICENSE`). When using this code or the derived
results, please cite the associated manuscript.
