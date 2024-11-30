"""
adapted from: https://github.com/greenelab/lab-website-template/blob/main/_cite/util.py
"""

import subprocess
import json
import argparse

# Special cases not handled by Manubot
database_urls = {
    "malacard": "https://www.malacards.org/card/{id}",
    "genereviews": "https://www.ncbi.nlm.nih.gov/books/{id}",
    }

def parse_args():
    """
    Parse command line arguments
    """

    parser = argparse.ArgumentParser(description="Generate full citations from ids with Manubot")
    parser.add_argument("input", help="JSON file of ids to generate citations for")
    parser.add_argument("output", help="Output JSON file off full citations")
    return parser.parse_args()

def cite_with_manubot(ids):
    """
    generate full citations from ids with Manubot
    """

    # output list of full citation details
    citations = []

    for _id in ids:
        # new citation
        citation = {}

        if _id.split(":")[0] in database_urls:
            print(_id)
            id_type = database_urls[_id.split(":")[0]]
            url = database_urls[_id.split(":")[0]].format(id=_id.split(":")[1])
            print(f"WARNING: Manubot does not support {_id}. Converted to URL: {url}")
            _id = "url:" + url

        # original id
        citation["id"] = _id

        # run Manubot
        try:
            commands = ["manubot", "cite", _id, "--log-level=WARNING"]
            output = subprocess.Popen(commands, stdout=subprocess.PIPE).communicate()
        except Exception as e:
            print("WARNING: Manubot could not generate citation")
            print(e)
            citations.append(citation)
            continue

        # parse results as json
        try:
            manubot = json.loads(output[0])[0]
        except Exception as e:
            print("WARNING: Couldn't parse Manubot response")
            print(e)
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
    # read input JSON
    with open(args.input, "r") as file:
        data = json.load(file)

    # write output JSON
    with open(args.output, "w") as file:
        json.dump(cite_with_manubot(data), file, indent=4)

#test = ["doi:10.1101/2023.10.09.23296582", "pmid:11246464"]
#print(json.dumps(cite_with_manubot(test), indent=4))

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()
    main(args)
