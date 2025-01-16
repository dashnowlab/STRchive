# Find the last updated date for each locus in the json

import git
import re
import time

repo = git.Repo(".")

ignore_fields = ["references", "additional_literature"]

locus_id = None
locus_date = None
for commit, lines in repo.blame("HEAD", "data/STRchive-loci.json"):
    line_date = time.gmtime(commit.committed_date)
    for line in lines:
        if line.startswith("}"):
            print(locus_id, locus_date.tm_year, locus_date.tm_mon, locus_date.tm_mday)
        elif line.startswith('{'):
            locus_id = None
            locus_date = None
        else:
            if line.strip().startswith('"id":'):
                locus_id = line.split(":")[-1].strip().replace('"', "")
            if not any(ignore_field in line for ignore_field in ignore_fields):
                # Get the latest date
                if locus_date is None:
                    locus_date = line_date
                elif line_date > locus_date:
                    locus_date = line_date
