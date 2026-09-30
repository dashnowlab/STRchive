
#!/usr/bin/env python3
"""
Retrieve PubMed Literature for monthly lit review
Usage:
  python get-literature.py STRchive-loci.json literature-dir out-citations.json out-loci-literature.json
"""

import csv
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOMART_HOSTS = ["https://useast.ensembl.org", "https://asia.ensembl.org", "https://www.ensembl.org"]
SYNONYM_FILE = "data/literature/gene_synonyms.tsv"  #Backup
SYNONYM_MAX_AGE_DAYS = float(os.environ.get("STRCHIVE_SYNONYM_MAX_AGE_DAYS", 180))

# Locus fields that can have citation markers
REF_FIELDS = ["age_onset", "mechanism_detail", "details", "prevalence_details",
              "year", "disease_description", "detection"]

# Synonyms to exclude
EXCLUDED_SYNONYMS = ["B37", "MHP", "MED", "DM", "DM1", "FA", "GAC", "SPD",
                     "PRP", "A1", "CCD", "PHP", "VCF", "PEM", "MCD", "EMA"]

# Synonyms BioMart misses
EXTRA_SYNONYMS = {"FMR1": "FMR-1", "NUTM2B-AS1": "LOC642361/NUTM2B-AS1", "ATXN10": "SCA10"}

#Repeat terms to inclide 
TERMS_REPEAT = ['"repeat expansion"[Title/Abstract]', '"repeat expansions"[Title/Abstract]',
                '"tandem repeat"[Title/Abstract]', '"tandem repeats"[Title/Abstract]',
                '"repeat sequence"[Title/Abstract]', '"repeat sequences"[Title/Abstract]',
                '"repeat length"[Title/Abstract]', '"repeat lengths"[Title/Abstract]',
                '"expansion"[Title]', '"expansions"[Title]', '"repeats"[Title]']

#Disease terms to inclide 
TERMS_DISEASE = ["disease*[Title/Abstract]", "disorder*[Title/Abstract]", "syndrome*[Title/Abstract]",
                 "patient*[Title/Abstract]", "proband*[Title/Abstract]"]
#types of publications
TERMS_PUBTYPE = ['"journal article"[Publication Type]', '"letter"[Publication Type]',
                 '"Case Reports"[Publication Type]']

TERMS_LANGUAGE = '"English"[Language]'
TERMS_EXCLUDE = '"review"[Publication Type]'

# New locus query
NEW_TERMS_REPEAT = ['"repeat expansion"[Title/Abstract]', '"tandem repeat"[Title/Abstract]']
NEW_TERMS_DISCOVERY = ['"discovered"[Title/Abstract]', '"identified"[Title/Abstract]',
                       '"causative"[Title/Abstract]', '"underlie"[Title/Abstract]',
                       '"basis"[Title/Abstract]']

NEW_TERMS_DISEASE = ['"disease"[Title/Abstract]', '"disorder"[Title/Abstract]',
                     '"syndrome"[Title/Abstract]', '"condition*"[Title/Abstract]']


def log(*args):
    print(*args, file=sys.stderr)


def dedupe(items):
    # unique order
    seen = {}
    for x in items:
        seen.setdefault(x, None)
    return list(seen)


def text_leaves(value):
   # Return all text leaves
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, str):
        yield value
    elif isinstance(value, (int, float)):
        yield str(value)
    elif isinstance(value, list):
        for item in value:
            yield from text_leaves(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from text_leaves(item)
    else:
        log(f"Ignoring unexpected {type(value).__name__} in reference field")



def http_post(url, params, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, data=urllib.parse.urlencode(params).encode(),
                headers={"User-Agent": "STRchive get-literature.py"})
            return urllib.request.urlopen(req, timeout=180).read().decode("utf-8", "replace")
        except (urllib.error.URLError, OSError) as e:
            if attempt == tries - 1:
                raise
            log(f"  request failed ({e}), retrying...")
            time.sleep(2 ** attempt)


def eutils(endpoint, params):
    params = dict(params)
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    if os.environ.get("NCBI_EMAIL"):
        params["email"] = os.environ["NCBI_EMAIL"]
    params["tool"] = "STRchive"
    text = http_post(EUTILS + endpoint, params)
    time.sleep(0.11 if os.environ.get("NCBI_API_KEY") else 0.34)  # rate limit
    return text



# synonyms

def synonym_file_age_days():
    #Age of synonym file 
    try:
        ts = subprocess.run(["git", "log", "-1", "--format=%ct", "--", SYNONYM_FILE],
                            capture_output=True, text=True, timeout=10).stdout.strip()
        if ts:
            return (time.time() - int(ts)) / 86400
    except Exception:
        pass
    return (time.time() - os.path.getmtime(SYNONYM_FILE)) / 86400


def read_synonym_file(genes):
    # Fallback synonyms from the local TSV
    if not os.path.exists(SYNONYM_FILE):
        sys.exit(f"BioMart is unreachable and there is no fallback file")

    age_days = synonym_file_age_days()
    log("Failed to create the mart object. Using previous gene synonyms from file:", SYNONYM_FILE)
    log(f"  last modified {age_days:.0f} days ago") # To better know when it was last updated
    if age_days > SYNONYM_MAX_AGE_DAYS and not os.environ.get("STRCHIVE_ALLOW_STALE_SYNONYMS"):
        sys.exit(f" {SYNONYM_FILE} is {age_days:.0f} days old (limit "
                 f"{SYNONYM_MAX_AGE_DAYS:.0f})")

    with open(SYNONYM_FILE, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = {"hgnc_symbol", "external_synonym"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{SYNONYM_FILE} is missing required column(s): {', '.join(sorted(missing))}")
        rows = [(r["hgnc_symbol"], r["external_synonym"] or "") for r in reader]

    rows.sort(key=lambda r: (r[1], r[0]))
    return dedupe(rows + [(g, "") for g in genes])


def get_synonyms(genes):
    #(hgnc_symbol, external_synonym) pairs 
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
           '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0"'
           ' count="" datasetConfigVersion="0.6">'
           '<Dataset name="hsapiens_gene_ensembl" interface="default">'
           f'<Filter name="hgnc_symbol" value="{",".join(genes)}"/>'
           '<Attribute name="hgnc_symbol"/><Attribute name="external_synonym"/>'
           '</Dataset></Query>')

    requested = set(genes)
    for host in BIOMART_HOSTS: #attempt to retrieve rows
        try: 
            text = http_post(host + "/biomart/martservice", {"query": xml}, tries=1)
            if "Query ERROR" in text:
                raise RuntimeError(text.strip().splitlines()[0])
            rows = [line.split("\t") for line in text.splitlines() if line.strip()]
            rows = [(r[0], r[1] if len(r) > 1 else "") for r in rows if r[0]]
            if not rows:
                raise RuntimeError("empty response")
            if not requested.intersection(sym for sym, _ in rows):
                raise RuntimeError("response contained none of the requested gene symbols")
            log(f"Retrieved {len(rows)} synonym rows from {host}") 
            return rows
        except Exception as e:
            log("Failed to create the mart object for host:", host, ". Trying another.")
            log(f"  ({e})")  # to actually log exceptions 

    return read_synonym_file(genes)


def report_synonym_coverage(rows, genes):
    # Log which genes will be searched under only their primary symbol
    have = {sym for sym, syn in rows if syn and syn != sym}
    absent = [g for g in genes if g not in have and g not in EXTRA_SYNONYMS]
    log(f"Synonyms found for {len(genes) - len(absent)}/{len(genes)} genes.")
    if absent:
        log(f"  Primary symbol only (no synonyms): {' '.join(absent)}")


def build_gene_terms(genes):
    # Return (gene_names, consolidated_strings)
    rows = get_synonyms(genes)
    report_synonym_coverage(rows, genes)

    # Genes with no synonym stand in for themselves
    rows = [(sym, syn or sym) for sym, syn in rows]
    rows.sort(key=lambda r: r[0]) 

    excluded = re.compile("|".join(EXCLUDED_SYNONYMS))
    rows = [r for r in rows if not excluded.search(r[1])]

    dropped = [g for g in genes if g not in {sym for sym, _ in rows} and g not in EXTRA_SYNONYMS]
    if dropped:
        log(f"  WARNING: dropped by EXCLUDED_SYNONYMS, will not be searched: {' '.join(dropped)}")

    rows = [(f'"{sym}"', f'"{syn}"') for sym, syn in rows]
    for sym, syn in EXTRA_SYNONYMS.items():
        if sym in genes:
            rows.append((f'"{sym}"', f'"{syn}"'))

    groups = {}
    for sym, syn in rows:
        groups.setdefault(sym, []).append(syn)

    consolidated = [" OR ".join(dedupe([sym] + syns)) for sym, syns in sorted(groups.items())]
    # Substitute to avoid unrelated PMIDs
    consolidated = [c.replace("BMD", "Becker muscular dystrophy") for c in consolidated]

    return dedupe([sym for sym, _ in rows]), consolidated



# PubMed

def build_query(terms_repeat, terms_middle, terms_disease):
    query = " ".join([
        "(", " OR ".join(terms_repeat), ")",
        "AND (", " OR ".join(terms_middle), ")",
        "AND", TERMS_LANGUAGE,
        "AND (", " OR ".join(terms_disease), ")",
        "AND (", " OR ".join(TERMS_PUBTYPE), ")",
        "NOT", TERMS_EXCLUDE,
    ])
    return re.sub(r"\s{2,}", " ", query)


def fetch_medline(query, out_prefix, gene=None):
    # Search PubMed and write medline records to <out_prefix>_batch_01.txt
    try:
        result = json.loads(eutils("esearch.fcgi", {
            "db": "pubmed", "term": query, "retmax": "10000", "retmode": "json"}))["esearchresult"]
    except Exception as e:
        log("batch pubmed download error: ", e, "")
        sys.exit(1)
    pmids = result.get("idlist", [])
    if gene is not None:
        log("Found", result.get("count", 0), "articles for gene:", gene, "") 
        if not pmids:
            log("Skipping fetch for:", gene, "")
            return ""
    else:
        log("Found", result.get("count", 0), "articles for new loci", "") 
        if not pmids:
            log("Skipping fetch for new loci - No articles found.")
            return ""

    chunks, dropped = [], 0
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
        try:
            chunks.append(eutils("efetch.fcgi", {"db": "pubmed", "id": ",".join(batch),
                                                 "rettype": "medline", "retmode": "text"}))
        except Exception as e:
            dropped += len(batch)
            log("entrez_fetch error:", e, "")
    if not chunks:
        return ""

    out_file = out_prefix + "_batch_01.txt"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(chunks) + "\n")  # Not sure if this is actually needed but it was the only way I could get it to exactly match R output. 
    if gene is not None:
        log(out_file, "")
    log("Full file path:", out_file, "")
    if not os.path.exists(out_file):
        log("Error: File not found -", out_file, "")
    return "\n".join(chunks)



def main(in_json, lit_dir, out_citations, out_loci_citations):
    log("Arguments: ", in_json, lit_dir, out_citations, out_loci_citations, "")

    with open(in_json, encoding="utf-8") as f:
        loci = json.load(f)

    bracketed = re.compile(r"(?<=\[)[^\]]+(?=\])")
    for locus in loci:
        found = []
        for field in REF_FIELDS:
            for text in text_leaves(locus.get(field)):
                found += bracketed.findall(text)
        locus["references"] = ",".join(dedupe(m.replace("; ", ",").strip() for m in found))

    genes = dedupe(locus["gene"] for locus in loci if locus.get("gene"))
    log("Searching for genes:", " ".join(genes), "")

    gene_names, consolidated = build_gene_terms(genes)

    os.makedirs(lit_dir, exist_ok=True)
    pmids_by_gene = {}
    for gene_name in gene_names:
        gene = gene_name.strip('"')
        log("Processing gene:", gene_name, "") 
        matched = [c for c in consolidated if re.search(gene_name, c, re.IGNORECASE)]
        terms_gene = [t + "[Title/Abstract]" for t in " OR ".join(matched).split(" OR ")]
        query = build_query(TERMS_REPEAT, terms_gene, TERMS_DISEASE)
        medline = fetch_medline(query, os.path.join(lit_dir, gene), gene=gene)
        pmids_by_gene[gene] = dedupe(re.findall(r"(?<=PMID- )\d+", medline))

    # add search results without PMIDs already there
    for locus in loci:
        found = ["@pmid:" + p for p in pmids_by_gene.get(locus.get("gene"), [])]
        cited = set(re.split(r",\s*", locus["references"]))
        redundant = [p for p in found if p in cited]
        if redundant:
            log("Removing the following elements from additional_literature: ", ", ".join(redundant), "")
        locus["additional_literature"] = ",".join(p for p in found if p not in cited)

    with open(out_loci_citations, "w", encoding="utf-8") as f:
        json.dump([{k: locus[k] for k in ("id", "additional_literature", "references")}
                   for locus in loci], f, separators=(",", ":"))

    # papers that maybe have new loci
    fetch_medline(build_query(NEW_TERMS_REPEAT, NEW_TERMS_DISCOVERY, NEW_TERMS_DISEASE),
                  os.path.join(lit_dir, "new_loci"))

    # Every citation
    citations = []
    for field in ("references", "additional_literature"):
        for locus in loci:
            citations += [t.strip()[1:] for t in locus[field].split(",")
                          if re.match(r"@\S+", t.strip())]
    with open(out_citations, "w", encoding="utf-8") as f:
        json.dump(dedupe(citations), f, separators=(",", ":"))
    log(f"Wrote all citations to {out_citations}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
