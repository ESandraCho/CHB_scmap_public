#!/usr/bin/env python3
"""
fetch_data.py — download the primary single-cell datasets from GEO. Reproducible
from scratch; no data is bundled, every file is fetched from its public accession.

Datasets (10x triplet, raw counts):
  GSE279630 / GSM8577299 — human fetal sinoatrial node, snRNA-seq (Lim et al. 2024)
  GSE297072 / GSM8983395 — dissected human fetal atrioventricular node tissue, snRNA-seq (Protze lab, unpublished; GSM8983395 is the dissected-tissue sample, NOT the hPSC-derived cells in the series)

Other datasets used by the analysis (adult Kanemaru conduction object; Sim 2021
developmental atlas; Cui/Wang STRT; Lázár first-trimester) are obtained from their
accessions; see DATA_SOURCES.md for accessions and the QC/normalization that produces
the *_qc.h5ad files under PROC_DIR.

Accessions are curated inputs. Output paths derive from config_chb (no hard-coding).
"""
from __future__ import annotations

import logging
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("fetch_data")

GEO_SUPPL = ("https://ftp.ncbi.nlm.nih.gov/geo/samples/{stem}nnn/{gsm}/suppl/"
             "{gsm}_{fname}")
# (subdir, GSM, file-prefix). 10x triplet: barcodes/features/matrix.
SAMPLES = [
    ("lim2024_fetal_san", "GSM8577299", "fetal"),
    ("lim2024_fetal_avn", "GSM8983395", "fetal"),
]
TRIPLET = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]


def download(gsm, prefix, dest):
    stem = gsm[:-3]
    dest.mkdir(parents=True, exist_ok=True)
    for part in TRIPLET:
        fname = f"{prefix}_{part}"
        out = dest / f"{gsm}_{fname}"
        if out.exists() and out.stat().st_size > 0:
            logger.info("exists %s", out.name); continue
        url = GEO_SUPPL.format(stem=stem, gsm=gsm, fname=fname)
        logger.info("downloading %s", url)
        urllib.request.urlretrieve(url, out)


def main():
    base = cfg.DATA_DIR
    for subdir, gsm, prefix in SAMPLES:
        download(gsm, prefix, base / subdir / gsm)
    logger.info("raw matrices under %s", base)


if __name__ == "__main__":
    main()
