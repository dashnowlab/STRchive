#!/usr/bin/env python3
"""
Check that renamed loci and criTRia curations keep redirects.

Compares locus IDs and criTRia Locus_IDs against a base git ref (e.g. origin/main).
Any ID that was removed must be listed in a locus's previous_ids, so the site
redirects the old page to the new one. Also reports the redirects this change adds.

Exits with an error if an ID was removed without a redirect, or if previous_ids
clash with current IDs or each other. Warns about criTRia curations whose
Locus_ID doesn't match a STRchive locus.
"""
import argparse
import json
import os
import subprocess
import sys

LOCI_PATH = "data/STRchive-loci.json"
CURATIONS_PATH = "data/criTRia-curations.json"

def load(path, ref=None):
    """Load JSON from the working tree, or from a git ref."""
    if ref is None:
        with open(path) as fh:
            return json.load(fh)
    text = subprocess.run(
        ["git", "show", f"{ref}:{path}"], capture_output=True, text=True, check=True
    ).stdout
    return json.loads(text)

def find_line(path, field, value):
    """1-based line number of '"field": value' in a JSON file, or None."""
    target = f'"{field}": {json.dumps(value, ensure_ascii=False)}'
    with open(path) as fh:
        for number, line in enumerate(fh, 1):
            if target in line:
                return number
    return None

def redirects(loci):
    """old id to new id, from previous_ids."""
    return {old: locus["id"] for locus in loci for old in locus.get("previous_ids") or []}

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--base", default="origin/main", help="git ref to compare against (default: %(default)s)")
    args = parser.parse_args()

    loci = load(LOCI_PATH)
    base_loci = load(LOCI_PATH, args.base)
    ids = {locus["id"] for locus in loci}
    curation_ids = {c["Locus_ID"] for c in load(CURATIONS_PATH)}
    base_ids = {locus["id"] for locus in base_loci}
    base_curation_ids = {c["Locus_ID"] for c in load(CURATIONS_PATH, args.base)}
    new = redirects(loci)
    old = redirects(base_loci)

    errors = []
    for kind, removed in [
        ("Locus", base_ids - ids),
        ("criTRia curation", base_curation_ids - curation_ids),
    ]:
        for removed_id in sorted(removed):
            if removed_id not in new:
                errors.append(
                    f"{kind} ID {removed_id} was removed or renamed. If it was renamed, "
                    f'add "{removed_id}" to previous_ids of the new locus so old links redirect.'
                )
    seen = {}
    for locus in loci:
        for prev in locus.get("previous_ids") or []:
            if prev in ids:
                errors.append(f"{locus['id']}: previous ID {prev} is also a current locus ID")
            if prev in seen:
                errors.append(f"{locus['id']}: previous ID {prev} is also listed for {seen[prev]}")
            seen[prev] = locus["id"]

    warnings = [
        (
            f"criTRia curation {curation_id} has no matching STRchive locus, so its page "
            "can't link to a locus. Its Locus_ID should match a locus id.",
            find_line(CURATIONS_PATH, "Locus_ID", curation_id),
        )
        for curation_id in sorted(curation_ids - ids)
    ]

    added = sorted((o, n) for o, n in new.items() if old.get(o) != n)
    lines = ["## Redirects", ""]
    if added:
        lines += ["New or changed redirects in this PR:", "", "| Old ID | New ID |", "|---|---|"]
        lines += [f"| {o} | {n} |" for o, n in added]
    else:
        lines.append("No new redirects in this PR.")
    if errors:
        lines += ["", "**Problems:**", ""] + [f"- {e}" for e in errors]
    if warnings:
        lines += ["", "**Warnings:**", ""] + [f"- {w}" for w, _ in warnings]
    report = "\n".join(lines)
    print(report)

    # Show the report on the GitHub Actions run summary
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as fh:
            fh.write(report + "\n")
    for e in errors:
        print(f"::error file={LOCI_PATH}::{e}")
    # With a line number, GitHub also shows the warning inline in the PR diff
    for w, line in warnings:
        location = f"file={CURATIONS_PATH}" + (f",line={line}" if line else "")
        print(f"::warning {location}::{w}")
    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
