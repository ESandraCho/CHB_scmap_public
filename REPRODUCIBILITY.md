# Reproducibility map — CHB developmental-expression paper

How to regenerate every figure and number in the manuscript, tracing each
**data file → script → output table → figure/manuscript** link. All analysis lives
in `scripts/`; curated inputs are in the `config_*.py` files; committed outputs
are in `results/tables/` and `results/figures/`. Seed 42 throughout.

---

## 1. Environment

- Python with: `anndata`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `leidenalg`, `matplotlib`
  (no `scanpy`, by design — numba/numpy conflict). Conda env name used during development:
  `structural_epitope`.
- No network needed to regenerate figures from the committed tables. Network **is** needed
  for the one-time data/annotation fetches (§2, §6).

## 2. Input data files

Single-cell/-nucleus objects are large and **not** committed; they are fetched from their
public accessions and QC-processed into the `*_qc.h5ad` objects the analysis reads.
Provenance, accessions and the processing that produces each object: `DATA_SOURCES.md`.

| Object (processed) | Accession | Role |
|---|---|---|
| `lim2024_fetal_san_qc.h5ad` | GEO GSE279630 / GSM8577299 (Lim/Protze) | fetal SAN (dissected) |
| `lim2024_fetal_avn_qc.h5ad` | GEO GSE297072 / GSM8983395 (Lim/Protze, same group) | fetal AVN (dissected) |
| `kanemaru_conduction_qc.h5ad` | Heart Cell Atlas (Kanemaru 2023) | adult SAN/AVN + working CM (Multiome nuclei used) |
| `sim2021_qc.h5ad` | GEO GSE156703 (Sim 2021) | fetal/child/adult whole-CM dev axis |
| `cui2019_gse106118/` (raw STRT) | GEO GSE106118 (Cui 2019) | fetal Smart-seq corroboration |
| `wang2020_gse109816/` (raw STRT) | GEO GSE109816 (Wang 2020) | adult Smart-seq corroboration |
| `lazar2025/*.h5` (TenX*) | Lázár 2025 (Mendeley fhtb99mdzd.1) | 13-donor first-trimester replication |

Fetch helpers (one-time, network): `fetch_data.py` (fetal-node 10x triplets),
`fetch_annotation_caches.py` (UniProt topology/localization caches), `fetch_genesets.py`
(EBI QuickGO gene sets). Curated inputs (no data): `config_chb.py`, `config_chb_handling_panel.py`,
`config_channel_subtypes.py`, `config_channel_classes.py`, `config_ca_accessibility.py`.

## 3. Run order

```
# one-time (network):
python fetch_data.py
python fetch_annotation_caches.py
python fetch_genesets.py
python build_channelopathy_classification.py     # -> channelopathy_classification.json (Fig 1 encoding)

# analysis (any order; each reads h5ad/raw + configs, writes results/tables + results/reports):
python fetal_pacemaker_vs_working.py
python adult_pacemaker_vs_working.py
python nodal_specificity_interaction.py
python ca_accessibility_analysis.py
python sr_activity_proxies.py
python channel_class_budget.py
python ca_subtype_dominance_dev.py
python k_repolarization_dev.py
python fetal_vs_adult_vulnerability.py
python lazar_dev_replication.py
python cav_subtype_pace_work_stats.py            # depends on node_channel_handling_expression.csv
python node_channel_handling_expression.py
python node_datasets_methods_table.py
python supp_dataset_analysis_table.py             # Supp Tables 1-2 as CSV (datasets + analysis map)
python build_supplementary_tables.py             # Supp Tables 1-3 -> CHB_supplementary_tables.xlsx (Excel)
python channel_subtype_conduction.py
python geneset_capture_qc.py                      # network (QuickGO) on first run
python platform_qc.py

# figures (read the tables above):
python make_figures.py                            # F1, F2, F3, S1, S2, S3
python make_node_handling_figure.py               # F4
```
Dependency note: `cav_subtype_pace_work_stats.py` reads `node_channel_handling_expression.csv`,
so run `node_channel_handling_expression.py` first. All other analysis scripts read only
data + configs.

## 4. Figure → script → table map

Figures are written to `public/results/figures/` as PNG + SVG.

| Figure | Panel(s) | Source table(s) | Producing script |
|---|---|---|---|
| **F1A** | antibody-phenotype schematic | `channelopathy_classification.json` (encoding) | `build_channelopathy_classification.py` |
| **F1B** | Cav-subtype composition + pace>work asterisks | `node_channel_handling_expression.csv`, `cav_subtype_pace_work_stats.csv` | `node_channel_handling_expression.py`, `cav_subtype_pace_work_stats.py` |
| **F1C** | Sim fetal/adult ratio | `fetal_vs_adult_vulnerability.json`, `ca_subtype_dominance_dev.json` | `fetal_vs_adult_vulnerability.py`, `ca_subtype_dominance_dev.py` |
| **F1D,E** | AVN pacemaker-vs-working (fetal, adult) | `fetal_pacemaker_vs_working.csv`, `adult_pacemaker_vs_working.csv` | `fetal_pacemaker_vs_working.py`, `adult_pacemaker_vs_working.py` |
| **F2A** | surface-accessible Ca fraction (+per-cell violins) | `ca_accessibility.csv`, `ca_accessibility_percell.csv` | `ca_accessibility_analysis.py` |
| **F2B** | NCX:RYR2 mode (+per-cell violins) | `sr_activity_proxies.csv`, `sr_activity_percell.csv` | `sr_activity_proxies.py` |
| **F2C** | 6-class channel budget (+Sim ladder) | `channel_class_budget.csv` | `channel_class_budget.py` |
| **F3A** | Cav subtype dev trajectory | `ca_subtype_dominance_dev.json` | `ca_subtype_dominance_dev.py` |
| **F3B** | RYR2 dev trajectory (+donor points) | `ca_subtype_dominance_dev.json` | `ca_subtype_dominance_dev.py` |
| **F3C** | inward:outward balance (+donor points) | `k_repolarization_dev.json` | `k_repolarization_dev.py` |
| **F4** | node-channel dot-plot (fetal/adult/Sim) | `node_channel_handling_expression.csv` | `make_node_handling_figure.py` |
| **S1** | platform QC (complexity, capture, depth) | `platform_qc.json`, `geneset_capture_qc.json` | `platform_qc.py`, `geneset_capture_qc.py` |
| **S2** | dataset reliability (Sim/Cui/Wang + Lázár) | `fetal_vs_adult_vulnerability.json`, `lazar_dev_replication.csv` | `fetal_vs_adult_vulnerability.py`, `lazar_dev_replication.py` |
| **S3** | SAN pacemaker-vs-working (fetal, adult) | `fetal_pacemaker_vs_working.csv`, `adult_pacemaker_vs_working.csv` | `fetal_pacemaker_vs_working.py`, `adult_pacemaker_vs_working.py` |

## 5. Methods / supplementary tables

`node_datasets_methods_table.py` → `node_datasets_methods.csv` / `_summary.csv` (+`.md`):
the datasets-and-cell-types table for the Methods section, with cell counts computed live
from the objects and provenance from `DATA_SOURCES.md`.

`supp_dataset_analysis_table.py` → `supp_datasets_provenance.csv` (**Supp Table 1**, all datasets)
and `supp_analysis_dataset_map.csv` (**Supp Table 2**, which analysis used which datasets / cell
types / method / figures), from the curated provenance + analysis map in
`config_chb_handling_panel.py`. **Supp Table 3** = `channelopathy_classification.tsv`
(`build_channelopathy_classification.py`).

`build_supplementary_tables.py` → `CHB_supplementary_tables.xlsx`: the single formatted Excel
workbook for the supplement, one sheet per table — **S1** datasets (all 7; cell counts and object
totals read live from `node_datasets_methods_summary.csv`, `fetal_vs_adult_vulnerability.json`,
`lazar_dev_replication.csv`, and the raw Cui/Wang/Lázár matrices), **S2** analysis-to-dataset/figure
map, **S3** the channel classification (`channelopathy_classification.tsv`). No data values are
hard-coded; only the S2 cross-walk (documentation) is curated in-script.

## 6. Other live outputs (cited in text, not a main figure)

| Output | Script | Used for |
|---|---|---|
| `nodal_specificity_interaction.{csv,md}` | `nodal_specificity_interaction.py` | constitutive-vs-fetal-specific interaction (Methods) |
| `channel_subtype_conduction.{csv,md}` | `channel_subtype_conduction.py` | subtype-ratio supporting analysis |
| `cav_subtype_pace_work_stats.{csv,md}` | `cav_subtype_pace_work_stats.py` | F1B pace>work statistics (Results) |

