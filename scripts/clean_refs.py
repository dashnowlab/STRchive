import json
import re

with open('data/STRchive_literature_PMID.json') as f:
    data = json.load(f)

    newdict = {}
    for gene, genedict in data.items():
        newdict[gene] = {}
        for pmid, pmidlist in genedict.items():
            newdict[gene][pmid] = {}
            for x in pmidlist:
                for key, value in x.items():
                    if key == "PublicationDate":
                        key = "Year"
                    value = '; '.join(value)
                    value = re.sub(' +', ' ', value)
                    newdict[gene][pmid][key] = value

with open('data/STRchive_PMID.json', 'w') as f:
    json.dump(newdict, f, indent=4)
