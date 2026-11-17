"""
adapted from: https://github.com/greenelab/lab-website-template/blob/main/_cite/util.py
"""

import subprocess
import json
import jsbeautifier
import argparse
import sys
from citation_utils import cite_with_ncbi, cite_mondo

# These types are fetched from NCBI E-utilities (see citation_utils.py) rather than Manubot
ncbi_types = ["pmid", "genereviews", "omim"]

# These types are fetched from the EBI Ontology Lookup Service (see citation_utils.py)
ols_types = ["mondo"]

# These types are generated locally from the ID, with no network lookup
local_types = ["gnomad"]

non_manubot_types = ncbi_types + ols_types + local_types

# Notes on citations made by this script rather than Manubot. Cached entries of
# non_manubot_types without one of these were made by Manubot and are looked up again.
own_note_prefixes = ("Generated from NCBI", "Generated from the EBI", "Generated locally")

# Special cases not handled by Manubot.
# In future should pull these from identifiers.js or the json spec

# These will be converted to the citation type "url" before being passed to Manubot
database_urls = {
    "malacard": "https://www.malacards.org/card/{id}",
    "stripy": "https://stripy.org/database/{id}",
    "gnomad": "https://gnomad.broadinstitute.org/short-tandem-repeat/{id}?dataset=gnomad_r4",
    "orphanet": "https://www.orpha.net/en/disease/detail/{id}",
    }

# Titles returned when a site serves a bot-check/landing page instead of the record.
# A citation with one of these (or an empty title) is treated as a failed lookup,
# so it is not reused from --cache and will be retried on the next run.
placeholder_title_substrings = [
    "checking your browser",
    "recaptcha",
    "ontology lookup service",
    "vérification de la connexion",
    "just a moment",
    "attention required",
    "access denied",
    ]
placeholder_titles_exact = [
    "gnomad",
    "omim",
    "orphanet",
    "malacards",
    "genereviews",
    "pubmed",
    ]

def bad_title_reason(title):
    """
    Return a reason string if the title is missing or a known placeholder, else None
    """
    t = (title or "").strip()
    if not t:
        return "empty title"
    tl = t.lower()
    if tl in placeholder_titles_exact:
        return f"placeholder title '{t}'"
    for s in placeholder_title_substrings:
        if s in tl:
            return f"placeholder title '{t}'"
    return None

def cite_gnomad(_id):
    """
    Build a citation for a gnomAD short tandem repeat page from its ID.
    The gnomAD site is a JavaScript app, so scraping it only ever returns the title "gnomAD".
    """
    locus = _id.split(":", 1)[1]
    return {
        "id": _id,
        "manubot_success": True,
        "link": database_urls["gnomad"].format(id=locus),
        "title": f"{locus} short tandem repeat",
        "type": "dataset",
        "doi": "",
        "authors": [],
        "publisher": "gnomAD v4",
        "issn": "",
        "date": "",
        "abstract": "",
        "language": "en",
        "note": "Generated locally from the gnomAD STR ID (no fetch)",
    }

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate full citations from ids with Manubot")
    parser.add_argument("input", help="JSON file of ids to generate citations for")
    parser.add_argument("output", help="Output JSON file off full citations")
    parser.add_argument("--cache", help="Previous citations JSON (e.g. STRchive-citations.json). Entries with a valid title are reused instead of being looked up again (see README).", default=None)
    # old name for --cache, kept so existing commands still work
    parser.add_argument("--append", dest="cache", help=argparse.SUPPRESS)
    parser.add_argument("--inloci", help="JSON input is in STRchive loci format. Only items in the references field will be queried. Recommend using with --cache", action='store_true')
    parser.add_argument("--incurations", help="JSON input is in criTRia curations format. Only items in the references field will be queried. Recommend using with --cache", action='store_true')
    parser.add_argument("--refresh", help="Look up all input IDs again. Cached entries are only used where a lookup fails.", action='store_true')
    parser.add_argument("--curations", help="Optional criTRia curations JSON file. References from these curations will be added to the IDs queried.", default=None)
    return parser.parse_args()

def cite_with_manubot(ids, cache_ids=None, refresh_cache=None, keep_cached=True):
    """
    generate full citations from ids with Manubot (or NCBI / locally for some types)
    cache_ids: ids to skip; refresh_cache: previous citations, used if a lookup fails;
    keep_cached: also output cached ids that aren't in ids
    """
    refresh_cache = refresh_cache or {}

    # remove duplicates, preserves order python 3.7+
    ids = list(dict.fromkeys(list(ids) + (list(refresh_cache) if keep_cached else [])))
    cache_ids = list(dict.fromkeys(cache_ids or []))

    # Remove ids that are already in the cache
    if cache_ids:
        ids = [i for i in ids if i not in cache_ids]

    # ids = ids[:10] #XXX for testing only

    # output list of full citation details
    citations = []

    # Group ids by type
    id_types = {}
    for _id in ids:
        id_type = _id.split(":")[0]

        if id_type not in id_types:
            id_types[id_type] = []
        id_types[id_type].append(_id)

    # NCBI types
    ncbi_ids = [i for t in ncbi_types for i in id_types.get(t, [])]
    if ncbi_ids:
        citations += list(cite_with_ncbi(ncbi_ids, bad_title_reason=bad_title_reason).values())

    # MONDO terms
    mondo_ids = id_types.get("mondo", [])
    if mondo_ids:
        citations += list(cite_mondo(mondo_ids).values())

    # gnomAD STR pages: generated locally, no network lookup
    gnomad_ids = id_types.get("gnomad", [])
    if gnomad_ids:
        sys.stderr.write(f"Generating {len(gnomad_ids)} gnomad citations locally\n")
        citations += [cite_gnomad(i) for i in gnomad_ids]

    # process ids by type (useful for testing)
    for id_type, id_list in id_types.items():
        if id_type in non_manubot_types:
            continue
        # if id_type != "gnomad": #XXX for testig only
        #     continue
        sys.stderr.write(f"Processing {len(id_list)} citations of type: {id_type}\n")

        for _id in id_list:
            # new citation
            citation = {}

            # original id
            citation["id"] = _id

            # Manubot success flag
            citation["manubot_success"] = False

            id_type = _id.split(":")[0]

            # Convert some types to URLs using the rules above
            if id_type in database_urls:
                url = database_urls[id_type].format(id=_id[len(id_type)+1:]) # effectively splitting on first colon only
                sys.stderr.write(f"WARNING: Manubot does not support {_id}. Converted to URL: {url}\n")
                _id = "url:" + url
                citation["link"] = url
            
            # run Manubot
            try:
                commands = ["manubot", "cite", _id]
                result = subprocess.run(commands, 
                    capture_output=True,
                    timeout=3)
                output = result.stdout
                # throwing away stderr for now to hide all the DeprecationWarnings

            except Exception as e:
                sys.stderr.write("WARNING: Manubot could not generate citation\n")
                sys.stderr.write(f"{e}\n")
                citation["note"] = f"WARNING: Manubot could not generate citation: {e}"
                citations.append(citation)
                continue
            
            # parse results as json
            try:
                manubot = json.loads(output)[0]
            except Exception as e:
                sys.stderr.write("WARNING: Couldn't parse Manubot response\n")
                sys.stderr.write(f"{e}\n")
                sys.stderr.write(f"{output}\n")
                citation["note"] = f"WARNING: Couldn't parse Manubot response: {e}"
                citations.append(citation)
                continue
            
            # title
            citation["title"] = get_safe(manubot, "title").strip()

            # type
            citation["type"] = get_safe(manubot, "type").strip()

            # doi
            citation["doi"] = get_safe(manubot, "DOI").strip()

            # authors
            citation["authors"] = []
            for author in get_safe(manubot, "author", {}):
                given = get_safe(author, "given").strip()
                family = get_safe(author, "family").strip()
                if given or family:
                    citation["authors"].append([given, family])

            # publisher
            container = get_safe(manubot, "container-title").strip()
            collection = get_safe(manubot, "collection-title").strip()
            publisher = get_safe(manubot, "publisher").strip()
            source = get_safe(manubot, "source").strip()
            citation["publisher"] = container or publisher or collection or source or ""
            citation["issn"] = get_safe(manubot, "ISSN")

            # dates
            # citation["accessed_date"] = get_date(get_safe(manubot, "accessed", {}))
            citation["date"] = get_date(get_safe(manubot, "issued", {}))

            # link
            citation["link"] = get_safe(manubot, "URL").strip()

            # abstract
            citation["abstract"] = get_safe(manubot, "abstract").strip()

            # language
            citation["language"] = get_safe(manubot, "language", "en").strip()

            # note
            citation["note"] = get_safe(manubot, "note").strip()

            # Manubot returned JSON, but it may be a bot-check page or have no title
            reason = bad_title_reason(citation["title"])
            if reason:
                sys.stderr.write(f"WARNING: {_id} returned {reason}; marking as failed\n")
                citation["note"] = f"WARNING: Manubot returned {reason}"
                citation["title"] = ""
                citations.append(citation)
                continue

            # successfully generated citation with manubot
            citation["manubot_success"] = True

            # add citation to list
            citations.append(citation)

    # if a lookup failed (e.g. a site is down), keep the previous citation instead
    n_fallback = 0
    for i, cite in enumerate(citations):
        if not cite["manubot_success"] and cite["id"] in refresh_cache:
            citations[i] = refresh_cache[cite["id"]]
            n_fallback += 1
    if n_fallback:
        sys.stderr.write(f"WARNING: lookup failed for {n_fallback} citations; kept previous versions\n")

    return citations


# get YYYY-MM-DD date string from wonky Manubot split-date format
def get_date(date_parts):
    # extract date part
    def get_part(date_parts, index):
        try:
            return int(date_parts["date-parts"][0][index])
        except Exception:
            return 0

    # date
    year = get_part(date_parts, 0)
    if year:
        # fallbacks for month and day
        month = get_part(date_parts, 1) or 1
        day = get_part(date_parts, 2) or 1
        return f"{year:04}-{month:02}-{day:02}"
    else:
        # if no year, consider date missing data
        return ""


def get_safe(item, path, default=""):
    """
    safely access value in nested lists/dicts
    """

    for part in str(path).split("."):
        try:
            part = int(part)
        except ValueError:
            part = part
        try:
            item = item[part]
        except (KeyError, IndexError, AttributeError, TypeError):
            return default
    return item

def main(args):
    """
    Expecting:

    Files:
    args.input: JSON file of ids to generate citations for
    args.output: Output JSON file of full citations

    Options:
    args.cache
    args.inloci
    """
    # read cached citations JSON
    cache_ids = []
    cache_json = []
    refresh_cache = {}
    if args.cache:
        try:
            with open(args.cache, "r") as file:
                cache_json_all = json.load(file)
                for cite in cache_json_all:
                    status = True
                    try:
                        status = cite["manubot_success"]
                    except KeyError:
                        try:
                            if cite["note"].startswith("WARNING"):
                                status = False
                        except KeyError:
                            status = False
                    # Don't keep entries whose title is empty or a placeholder,
                    # even if they were previously marked successful
                    if status == True and bad_title_reason(cite.get("title")):
                        status = False
                    cite["manubot_success"] = status
                    id_type = cite["id"].split(":")[0]
                    made_by_manubot = id_type in non_manubot_types and not (cite.get("note") or "").startswith(own_note_prefixes)
                    if status == True and (args.refresh or id_type in local_types or made_by_manubot):
                        # regenerate these; keep the old version as a fallback only
                        refresh_cache.setdefault(cite["id"], cite)
                    elif status == True:
                        if cite["id"] not in cache_ids:
                            # only keep the first entry for each id, to avoid duplicates
                            cache_ids.append(cite["id"])
                            cache_json.append(cite)
        except Exception as e:
            sys.stderr.write(f"WARNING: Couldn't read cache JSON\n")
            sys.stderr.write(f"{e}\n")
            cache_json = []
            refresh_cache = {}
    sys.stderr.write(f"Skipping lookup of {len(cache_ids)} existing citations\n")
    sys.stderr.write(f"Re-fetching {len(refresh_cache)} existing citations (previous versions kept as fallback)\n")

    # read input JSON
    with open(args.input, "r") as file:
        if args.inloci:
            data = []
            loci_data = json.load(file)
            for record in loci_data:
                data = data + [x.lstrip('@') for x in record['references']]
                data = data + [x.lstrip('@') for x in record['additional_literature']]
        elif args.incurations:
            data = []
            curations_data = json.load(file)
            for record in curations_data:
                data = data + [x.lstrip('@') for x in record.get('references', [])]
        else:
            data = json.load(file)

    # add references from supplementary curations file
    if args.curations:
        with open(args.curations, "r") as f:
            curations_data = json.load(f)
            for record in curations_data:
                data = data + [x.lstrip('@') for x in (record.get('references') or [])]

    # write output JSON
    with open(args.output, "w") as file:
        new_json = cite_with_manubot(data, cache_ids=cache_ids, refresh_cache=refresh_cache,
                                     keep_cached=not args.refresh)
        options = jsbeautifier.default_options()
        options.indent_size = 2
        options.brace_style="expand"
        file.write(jsbeautifier.beautify(json.dumps(cache_json + new_json), options))

#test = ["doi:10.1101/2023.10.09.23296582", "pmid:11246464"]
#print(json.dumps(cite_with_manubot(test), indent=4))

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    main(args)
