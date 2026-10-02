#!/usr/bin/env python3
"""
Add HPO terms to STRchive-loci.json.

Terms come from the HPO annotations (phenotype.hpoa) for each locus's OMIM and
Orphanet ids. Orphanet annotations are limited to terms at or above a minimum
frequency (default: occasional, i.e. >=5% of patients). OMIM annotations rarely
record a frequency, so all are kept. Term names come from the HPO ontology.
"""
import argparse
import json
import sys
from pathlib import Path

import jsbeautifier
import pandas as pd
import requests

REPO = Path(__file__).resolve().parent.parent
LOCI_PATH = REPO / "data" / "STRchive-loci.json"

HPO_RELEASE = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
HPOA_URL = f"{HPO_RELEASE}/phenotype.hpoa"
HPO_ONTOLOGY_URL = f"{HPO_RELEASE}/hp.json"

# HPO frequency terms used by Orphanet annotations, most to least frequent
ORPHANET_FREQUENCIES = {
    "obligate": "HP:0040280",       # 100%
    "very-frequent": "HP:0040281",  # 80-99%
    "frequent": "HP:0040282",       # 30-79%
    "occasional": "HP:0040283",     # 5-29%
    "very-rare": "HP:0040284",      # 1-4%
}

def write_loci(loci, loci_path):
    """Write loci to path"""
    options = jsbeautifier.default_options()
    options.indent_size = 2
    options.brace_style = "expand"
    with open(loci_path, "w") as fh:
        fh.write(jsbeautifier.beautify(json.dumps(loci, ensure_ascii=False), options))
        fh.write("\n")

def normalise_id(value, prefix):
    """add prefix e.g. 'OMIM:' to id"""
    if value is None:
        return None
    value = str(value).strip()
    if not value:
        return None
    return value if value.startswith(prefix) else f"{prefix}{value}"

def as_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]

def load_hpo_names():
    """HPO id to HPO name."""
    graph = requests.get(HPO_ONTOLOGY_URL, timeout=120).json()["graphs"][0]
    names = {}
    for node in graph["nodes"]:
        if node.get("lbl") and "/HP_" in node["id"]:
            names[node["id"].rsplit("/", 1)[1].replace("_", ":")] = node["lbl"]
    return names

def load_annotations(min_frequency, hpo_names):
    """disease id (OMIM:/ORPHA:) to {HPO id: HPO name}."""
    levels = list(ORPHANET_FREQUENCIES)
    keep = {ORPHANET_FREQUENCIES[f] for f in levels[: levels.index(min_frequency) + 1]}
    hpoa = pd.read_csv(HPOA_URL, sep="\t", comment="#", dtype=str)
    is_omim = hpoa["database_id"].str.startswith("OMIM:", na=False)
    is_orpha = hpoa["database_id"].str.startswith("ORPHA:", na=False)
    hpoa = hpoa.loc[
        hpoa["qualifier"].isna()   # skip "NOT" annotations
        & (hpoa["aspect"] == "P")  # phenotypic abnormalities only
        & (is_omim | (is_orpha & hpoa["frequency"].isin(keep)))
    ]
    disease_to_hpo = {}
    missing = set()
    for disease, hpo_id in hpoa[["database_id", "hpo_id"]].itertuples(index=False):
        if hpo_id not in hpo_names:
            missing.add(hpo_id)
            continue
        disease_to_hpo.setdefault(disease, {})[hpo_id] = hpo_names[hpo_id]
    if missing:
        print(f"WARNING: skipped {len(missing)} HPO terms with no name in the ontology", file=sys.stderr)
    return disease_to_hpo

def annotate(loci, disease_to_hpo, replace=False):
    """Set hpo terms in place."""
    n_annotated = n_changed = n_no_ids = 0

    for locus in loci:
        before = locus.get("hpo_terms")
        terms = {}
        if not replace:
            for existing in as_list(before):
                existing = str(existing).strip()
                if existing:
                    terms[existing.split()[0]] = existing

        diseases = [normalise_id(i, "OMIM:") for i in as_list(locus.get("omim"))]
        diseases += [normalise_id(i, "ORPHA:") for i in as_list(locus.get("orphanet"))]
        diseases = list(filter(None, diseases))
        if not diseases:
            n_no_ids += 1
        for disease in diseases:
            for hpo_id, hpo_name in disease_to_hpo.get(disease, {}).items():
                terms[hpo_id] = f"{hpo_id} {hpo_name}"

        after = sorted(terms.values()) or before
        locus["hpo_terms"] = after
        n_annotated += bool(after)
        n_changed += after != before
    return n_annotated, n_changed, n_no_ids

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
    parser.add_argument(
        "--orphanet-min-frequency",
        choices=list(ORPHANET_FREQUENCIES),
        default="occasional",
        help="least frequent Orphanet annotation to include (default: %(default)s)",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="replace existing hpo terms instead of adding to them",
    )
    args = parser.parse_args()
    loci_path = args.loci.resolve()
    if not loci_path.exists():
        sys.exit(f"ERROR: {loci_path} not found")
    with open(loci_path) as fh:
        loci = json.load(fh)

    disease_to_hpo = load_annotations(args.orphanet_min_frequency, load_hpo_names())
    n_annotated, n_changed, n_no_ids = annotate(loci, disease_to_hpo, args.replace)
    print(f"Loci: {len(loci)}")
    print(f"  with HPO terms after update: {n_annotated}")
    print(f"  hpo terms changed:           {n_changed}")
    print(f"  no omim or orphanet id:      {n_no_ids}")
    write_loci(loci, loci_path)
    print(f"Wrote {loci_path}")


if __name__ == "__main__":
    main()
