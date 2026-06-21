# Data sources

All primary data are public. No data files are bundled in this repository; they are
fetched from their accessions by `scripts/fetch_data.py` and processed into the
`*_qc.h5ad` objects the analysis reads. Annotation caches are fetched by
`scripts/fetch_annotation_caches.py`.

## Single-cell datasets

| Object (PROC_DIR) | Accession | Tissue / stage | Chemistry | Source |
|---|---|---|---|---|
| `lim2024_fetal_san_qc.h5ad` | GEO GSE279630 / GSM8577299 | human fetal sinoatrial node, GW19 | 10x snRNA-seq | Lim … Protze 2024 |
| `lim2024_fetal_avn_qc.h5ad` | GEO GSE297072 / GSM8983395 (dissected fetal-tissue sample only; NOT the hPSC-derived samples in the series) | human fetal atrioventricular node | 10x snRNA-seq | Lohbihler, Lim … Protze (GSE297072) |
| `kanemaru_conduction_qc.h5ad` | Heart Cell Atlas — conduction object (heartcellatlas.org) | adult human SAN/AVN + working CM | 10x snRNA + Multiome-RNA + scRNA; the SAN_P_cell/AVN_P_cell pacemaker states used here are 100% single-nucleus Multiome-RNA | Kanemaru et al. 2023 Nature |
| `sim2021_qc.h5ad` | GEO GSE156703 | human heart fetal/child/adult CM (within-study dev axis) | 10x snRNA | Sim et al. 2021 |
| `cui2019_gse106118/` (raw STRT counts) | GEO GSE106118 | human heart fetal CM, gestational wk 5–25 | STRT / non-10x | Cui et al. 2019 |
| `wang2020_gse109816/` (raw STRT matrix) | GEO GSE109816 | human adult heart CM (LA + LV) | STRT / non-10x | Wang et al. 2020 |
| `lazar2025/*.h5` (TenX*) | Mendeley Data 10.17632/fhtb99mdzd.1 (open processed; raw under EGA EGAS50000001029) | developing human heart, first trimester (pcw 5.5–14) | 10x snRNA | Lázár et al. 2025 |

## How the data are obtained and processed

- `fetch_data.py` downloads the two fetal-node 10x triplets (Lim SAN/AVN) directly from GEO.
- The Cui/Wang STRT matrices are downloaded from their GEO accessions and read directly (gzip),
  TP10K-normalized in-script; no `_qc.h5ad` is built for them.
- The adult (Kanemaru), developmental (Sim 2021) and first-trimester (Lázár) objects are large
  processed objects obtained from their accessions and put through standard QC (min genes,
  max %mt) + TP10K (+ log1p where stated) to produce the `*_qc.h5ad` files the analysis reads.
  QC thresholds are in `scripts/config_chb.py` (no hard-coded per-run values).

The developmental axis uses Sim 2021 as a within-study 10x contrast (batch-clean), corroborated
by the two independent non-droplet STRT datasets (Cui fetal, Wang adult). There is no Smart-seq
dataset of dissected conduction tissue, so the nodal arm rests on single-nucleus data plus
within-donor robustness checks (`platform_qc.py`, Fig S1) rather than cross-platform replication.

## Annotation caches (fetched live from UniProt)

| Cache | Source | Used by |
|---|---|---|
| `ca_localization_cache/subcellular.json` | UniProt subcellular-location field | `ca_accessibility_analysis.py` (surface vs intracellular) |

Gene→accession maps are curated inputs in `scripts/fetch_annotation_caches.py`.

## Gene-set capture QC (Fig S1)

`fetch_genesets.py` fetches the electrophysiology (GO:0005216 + GO:0022857) and structural
(GO:0030017 + GO:0008307) human gene sets live from the EBI QuickGO API and writes
`geneset_cache/genesets.json` (housekeeping = curated Eisenberg & Levanon 2013).
`geneset_capture_qc.py` then compares snRNA vs Smart-seq capture of these sets, and the cardiac
ion-channel subtype ratios, in stage-matched working cardiomyocytes (→ `geneset_capture_qc.json`).
Cardiac-vs-non-cardiac subtype membership is curated in `config_channel_subtypes.py`.

## Curated scientific inputs (in config, not results)

- Channel functional-subtype taxonomy: `scripts/config_channel_subtypes.py`
  (Nerbonne & Kass 2005; Amin et al. 2010; Catterall 2011).
- Ca-machinery subcellular accessibility: `scripts/config_ca_accessibility.py`
  (UniProt + standard channel biology; Cav family classified surface).
- Anti-Ro Ca-channel target panel and channelopathy classification: `scripts/config_chb.py`,
  `scripts/build_channelopathy_classification.py` (Lazzerini et al. 2017; Qu et al. 2019;
  Capecchi et al. 2019; Strandberg et al. 2013; Benjamin et al. 2025; Lazzerini & Boutjdir 2025).
