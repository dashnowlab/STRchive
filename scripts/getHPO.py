#!/usr/bin/env python3
"""
Add HPO terms to STRchive-loci.json.
"""
import argparse
import json
import sys
from pathlib import Path

import jsbeautifier
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
LOCI_PATH = REPO / "data" / "STRchive-loci.json"

HPO_URL = ("https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/disease_to_phenotypic_feature_association.all.tsv.gz")

def write_loci(loci, loci_path):
    """Write loci to path"""
    options = jsbeautifier.default_options()
    options.indent_size = 2
    options.brace_style = "expand"
    with open(loci_path, "w") as fh:
        fh.write(jsbeautifier.beautify(json.dumps(loci, ensure_ascii=False), options))
        fh.write("\n")

def normalise_mondo(m):
    """add 'mondo' to name"""
    if m is None:
        return None
    m = str(m).strip()
    if not m:
        return None
    return m if m.startswith("MONDO:") else f"MONDO:{m}"

def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]

def load_monarch():
    """mondo id to {HPO id: HPO name}."""
    hpo = pd.read_csv(
        HPO_URL,
        sep="\t",
        usecols=["subject", "object", "object_label"],
        dtype=str,
    )
    hpo = hpo.loc[
        hpo["subject"].str.startswith("MONDO:", na=False)
        & hpo["object"].str.startswith("HP:", na=False)
        & hpo["object_label"].notna()
    ]
    mondo_to_hpo = {}
    for mondo, hpo_id, hpo_name in hpo[
        ["subject", "object", "object_label"]
    ].itertuples(index=False):
        mondo_to_hpo.setdefault(mondo, {})[hpo_id] = hpo_name
    return mondo_to_hpo

def annotate(loci, mondo_to_hpo):
    """Set hpo terms in place."""
    n_annotated = n_changed = n_no_mondo = 0

    for locus in loci:
        before = locus.get("hpo_terms")
        terms = {}
        for existing in as_list(before):
            existing = str(existing).strip()
            if existing:
                terms[existing.split()[0]] = existing

        mondos = list(filter(None, map(normalise_mondo, as_list(locus.get("mondo")))))
        if not mondos:
            n_no_mondo += 1
        for mondo in mondos:
            for hpo_id, hpo_name in mondo_to_hpo.get(mondo, {}).items():
                terms[hpo_id] = f"{hpo_id} {hpo_name}"

        after = sorted(terms.values()) or before
        locus["hpo_terms"] = after
        n_annotated += bool(after)
        n_changed += after != before
    return n_annotated, n_changed, n_no_mondo

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
        loci = json.load(fh)

    n_annotated, n_changed, n_no_mondo = annotate(loci, load_monarch())
    print(f"Loci: {len(loci)}")
    print(f"  with HPO terms after update: {n_annotated}")
    print(f"  hpo terms changed:           {n_changed}")
    print(f"  no mondo id:                 {n_no_mondo}")
    write_loci(loci, loci_path)
    print(f"Wrote {loci_path}")


if __name__ == "__main__":
    main()
