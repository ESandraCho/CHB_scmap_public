#!/usr/bin/env python3
"""
fetch_annotation_caches.py — regenerate the curated annotation caches from primary
databases (UniProt). Run once before the analysis scripts; reproducible from scratch.

Produces:
  <DATA_DIR>/ca_localization_cache/subcellular.json — UniProt subcellular-location
      field per Ca-handling gene; used to classify surface vs intracellular.

Gene->accession maps are curated inputs (config). No hard-coded results — every
value is fetched live from UniProt at run time.
"""
from __future__ import annotations

import json
import logging
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("fetch_caches")

# Curated gene -> UniProt accession (canonical human). Inputs, not results.
_CAV_ACC = {
    "CACNA1C": "Q13936", "CACNA1D": "Q01668", "CACNA1G": "O43497", "CACNA1H": "O95180",
    "CACNA1I": "Q9P0X4", "CACNA1A": "O00555", "CACNA1B": "Q00975", "CACNA1E": "Q15878",
}
CA_MACHINERY_ACC = {
    **_CAV_ACC, "CACNA1S": "Q13698",
    "RYR2": "Q92736", "RYR3": "Q15413",
    "ITPR1": "Q14643", "ITPR2": "Q14571", "ITPR3": "Q14573",
    "ATP2A2": "P16615", "ATP2A1": "O14983",
    "ORAI1": "Q96D31", "ORAI2": "Q96SN7", "ORAI3": "Q9BRQ5",
    "STIM1": "Q13586", "STIM2": "Q9P246",
    "TRPM4": "Q8TD43", "TRPM7": "Q96QT4", "TRPC1": "P48995", "TRPC3": "Q13507",
    "TRPC6": "Q9Y210", "SLC8A1": "P32418", "TMEM38A": "Q9H6F2", "TMEM38B": "Q9NVV0",
}
UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.json?fields={fields}"


def fetch_json(acc, fields, retries=3):
    url = UNIPROT.format(acc=acc, fields=fields)
    for i in range(retries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=30))
        except Exception as e:
            if i == retries - 1:
                logger.warning("fetch failed %s: %s", acc, e)
                return None
            time.sleep(1.5)


def fetch_localization(data_dir):
    out = data_dir / "ca_localization_cache"
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for gene, acc in CA_MACHINERY_ACC.items():
        j = fetch_json(acc, "cc_subcellular_location")
        locs = []
        if j is not None:
            for c in j.get("comments", []):
                if c.get("commentType") == "SUBCELLULAR LOCATION":
                    for s in c.get("subcellularLocations", []):
                        locs.append(s.get("location", {}).get("value", ""))
        results[gene] = {"uniprot": acc, "locations": list(dict.fromkeys(locs))}
        logger.info("localization %s: %s", gene, "; ".join(results[gene]["locations"]) or "none")
    (out / "subcellular.json").write_text(json.dumps(results, indent=1))


def main():
    data_dir = cfg.DATA_DIR
    data_dir.mkdir(parents=True, exist_ok=True)
    fetch_localization(data_dir)
    logger.info("caches written under %s", data_dir)


if __name__ == "__main__":
    main()
