#!/usr/bin/env python3
"""
Add OMIM gene IDs to STRchive-loci.json.

Looks up each locus's gene symbol in HGNC and fills omim_gene where it is empty.
Existing values are left unchanged.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import jsbeautifier
import requests

REPO = Path(__file__).resolve().parent.parent
LOCI_PATH = REPO / "data" / "STRchive-loci.json"

HGNC_URL = "https://rest.genenames.org/fetch/symbol/{symbol}"

def write_loci(loci, loci_path):
    """Write loci to path"""
    options = jsbeautifier.default_options()
    options.indent_size = 2
    options.brace_style = "expand"
    with open(loci_path, "w") as fh:
        fh.write(jsbeautifier.beautify(json.dumps(loci, ensure_ascii=False), options))
        fh.write("\n")

def hgnc_omim(symbol):
    """OMIM gene ID(s) for an approved HGNC gene symbol."""
    response = requests.get(
        HGNC_URL.format(symbol=symbol),
        headers={"Accept": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    docs = response.json()["response"]["docs"]
    if not docs:
        return None
    return docs[0].get("omim_id", [])

def with_omim_gene(locus):
    """Return locus with an omim_gene key placed after omim (schema order)."""
    if "omim_gene" in locus:
        return locus
    out = {}
    for key, value in locus.items():
        out[key] = value
        if key == "omim":
            out["omim_gene"] = []
    out.setdefault("omim_gene", [])
    return out

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--loci",
        type=Path,
        default=LOCI_PATH,
        help="path to STRchive-loci.json (default: %(default)s)",
    )
    args = parser.parse_args()
    loci_path = args.loci.resolve()
    if not loci_path.exists():
        sys.exit(f"ERROR: {loci_path} not found")
    with open(loci_path) as fh:
        loci = [with_omim_gene(locus) for locus in json.load(fh)]

    n_filled = 0
    cache = {}
    for locus in loci:
        if locus.get("omim_gene"):
            continue
        symbol = locus.get("gene")
        if symbol not in cache:
            cache[symbol] = hgnc_omim(symbol)
            time.sleep(0.1)  # HGNC asks for <10 requests per second
        omim = cache[symbol]
        if omim is None:
            print(f"WARNING: {locus['id']}: {symbol} is not an approved HGNC symbol", file=sys.stderr)
        elif not omim:
            print(f"WARNING: {locus['id']}: no OMIM gene ID in HGNC for {symbol}", file=sys.stderr)
        else:
            locus["omim_gene"] = omim
            n_filled += 1

    print(f"Loci: {len(loci)}")
    print(f"  omim_gene filled: {n_filled}")
    write_loci(loci, loci_path)
    print(f"Wrote {loci_path}")


if __name__ == "__main__":
    main()
