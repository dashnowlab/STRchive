"""
adapted from: https://github.com/greenelab/lab-website-template/blob/main/_cite/util.py
"""

import subprocess
import json
import argparse
import sys

# Special cases not handled by Manubot
database_urls = {
    "malacard": "https://www.malacards.org/card/{id}",
    "genereviews": "https://www.ncbi.nlm.nih.gov/books/{id}",
    "stripy": "https://stripy.org/database/{id}",
    "gnomad": "https://gnomad.broadinstitute.org/short-tandem-repeat/{id}?dataset=gnomad_r4",
    "omim": "https://omim.org/entry/{id}",
    }

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate full citations from ids with Manubot")
    parser.add_argument("input", help="JSON file of ids to generate citations for")
    parser.add_argument("output", help="Output JSON file off full citations")
    parser.add_argument("--append", help="JSON file of previous citations to be skipped and appended to the output", default=None)
    return parser.parse_args()

def cite_with_manubot(ids, append_ids=None):
    """
    generate full citations from ids with Manubot
    """

    # Remove ids that are already in the append list
    if append_ids:
        ids = [i for i in ids if i not in append_ids]

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

    # process ids by type (useful for testing)
    for id_type, id_list in id_types.items():
        # if id_type != "gnomad": #XXX for testig only
        #     continue
        sys.stderr.write(f"Processing {len(id_list)} citations of type: {id_type}\n")

        for _id in id_list:
            # new citation
            citation = {}

            # original id
            citation["id"] = _id

            id_type = _id.split(":")[0]

            # Convert some types to URLs using the rules above
            if id_type in database_urls:
                url = database_urls[id_type].format(id=_id[len(id_type)+1:])
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
                sys.stderr.write(f"{output[0]}\n")
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
            citation["issued_date"] = get_date(get_safe(manubot, "issued", {}))

            # link
            citation["link"] = get_safe(manubot, "URL").strip()

            # abstract
            citation["abstract"] = get_safe(manubot, "abstract").strip()

            # language
            citation["language"] = get_safe(manubot, "language", "en").strip()

            # note
            citation["note"] = get_safe(manubot, "note").strip()

            # add citation to list
            citations.append(citation)

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
        month = get_part(date_parts, 1) or "01"
        day = get_part(date_parts, 2) or "01"
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
    args.input: JSON file of ids to generate citations for
    args.output: Output JSON file of full citations
    """
    # read "append" JSON
    append_ids = []
    if args.append:
        try:
            with open(args.append, "r") as file:
                append_json = json.load(file)
                append_ids = [cite["id"] for cite in append_json]
        except Exception as e:
            sys.stderr.write(f"WARNING: Couldn't read append JSON\n")
            sys.stderr.write(f"{e}\n")
            append_json = []
    sys.stderr.write(f"Skipping lookup of {len(append_ids)} existing citations\n")

    # read input JSON
    with open(args.input, "r") as file:
        data = json.load(file)

    # write output JSON
    with open(args.output, "w") as file:
        new_json = cite_with_manubot(data, append_ids=append_ids)
        json.dump(append_json + new_json, file, indent=4)

#test = ["doi:10.1101/2023.10.09.23296582", "pmid:11246464"]
#print(json.dumps(cite_with_manubot(test), indent=4))

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    main(args)
