#!/usr/bin/env python3
"""
fetch_genesets.py — regenerate the gene sets used by the platform gene-capture QC
(geneset_capture_qc.py / Fig method-support) from primary annotation sources.

Produces <DATA_DIR>/geneset_cache/genesets.json with three human gene-symbol sets:

  electrophysiology : GO ion-channel + transporter activity (GO:0005216 ion channel
                      activity, GO:0022857 transmembrane transporter activity), human,
                      via the EBI QuickGO annotation API (includes descendant terms).
  structural        : GO sarcomere / contractile-fibre structural constituents
                      (GO:0030017 sarcomere, GO:0008307 structural constituent of
                      muscle), human, QuickGO.
  housekeeping      : curated validated human housekeeping genes (Eisenberg & Levanon
                      2013, Trends Genet — the standard stably-expressed reference set).

The point of the QC is the criticism that snRNA-seq under-captures low-abundance genes;
electrophysiology genes (ion channels) are typically lower-abundance than structural
(sarcomere) genes, so EP-vs-structural-vs-housekeeping capture is the informative contrast.

Gene-set DEFINITIONS (which GO terms, which HK source) are curated inputs; the gene
MEMBERSHIPS are fetched live. Run once before geneset_capture_qc.py.

Env: structural_epitope (urllib only).
"""
from __future__ import annotations

import json
import logging
import time
import urllib.request
from pathlib import Path

import config_chb as cfg

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("fetch_genesets")

# Curated GO-term definitions per set (inputs, not results).
GO_SETS = {
    "electrophysiology": ["GO:0005216", "GO:0022857"],   # ion channel activity; transmembrane transporter activity
    "structural":        ["GO:0030017", "GO:0008307"],   # sarcomere; structural constituent of muscle
}

# Eisenberg & Levanon 2013 housekeeping genes — a compact, widely-used validated subset
# (curated input; the full list is ~3800, this stable core is sufficient as a reference).
HOUSEKEEPING = [
    "ACTB", "GAPDH", "B2M", "HPRT1", "PGK1", "PPIA", "RPL13A", "RPLP0", "TBP",
    "UBC", "YWHAZ", "SDHA", "TFRC", "GUSB", "HMBS", "PGAM1", "RPS18", "ACTG1",
    "EEF1A1", "EIF4A2", "GPI", "PSMB4", "RPL19", "RPL27", "RPS29", "TUBB",
    "VPS29", "CHMP2A", "EMC7", "REEP5", "SNRPD3", "VCP", "RAB7A", "C1orf43",
]

QUICKGO = ("https://www.ebi.ac.uk/QuickGO/services/annotation/search?"
           "goId={go}&taxonId=9606&geneProductType=protein&limit=200&page={page}")


def fetch_go_symbols(go_id, max_pages=60):
    """All human gene symbols annotated to a GO term (+ descendants), paginated."""
    symbols, page = set(), 1
    while page <= max_pages:
        url = QUICKGO.format(go=go_id, page=page)
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        for attempt in range(3):
            try:
                j = json.load(urllib.request.urlopen(req, timeout=30))
                break
            except Exception as e:
                if attempt == 2:
                    logger.warning("fetch failed %s p%d: %s", go_id, page, e)
                    return symbols
                time.sleep(1.5)
        results = j.get("results", [])
        for r in results:
            s = r.get("symbol")
            if s:
                symbols.add(s)
        total = j.get("numberOfHits", 0)
        if page * 200 >= total or not results:
            break
        page += 1
        time.sleep(0.2)
    return symbols


def main():
    out = cfg.DATA_DIR / "geneset_cache"
    out.mkdir(parents=True, exist_ok=True)
    sets = {}
    for name, gos in GO_SETS.items():
        syms = set()
        for go in gos:
            s = fetch_go_symbols(go)
            logger.info("%s %s: %d symbols", name, go, len(s))
            syms |= s
        sets[name] = sorted(syms)
        logger.info("%s TOTAL: %d genes", name, len(syms))
    sets["housekeeping"] = sorted(HOUSEKEEPING)
    logger.info("housekeeping (curated Eisenberg 2013): %d genes", len(sets["housekeeping"]))

    meta = {
        "sources": {
            "electrophysiology": "GO:0005216 + GO:0022857 (QuickGO, human, incl. descendants)",
            "structural": "GO:0030017 + GO:0008307 (QuickGO, human)",
            "housekeeping": "Eisenberg & Levanon 2013 Trends Genet (curated stable core)",
        },
        "genesets": sets,
        "counts": {k: len(v) for k, v in sets.items()},
    }
    (out / "genesets.json").write_text(json.dumps(meta, indent=1))
    logger.info("wrote %s | counts=%s", (out / "genesets.json").name, meta["counts"])


if __name__ == "__main__":
    main()
