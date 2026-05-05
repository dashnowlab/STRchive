# Read criTRia curations from a json file and update them applying curation rules.
import sys
import argparse
import doctest
import os
import time
import json
import math
import jsbeautifier
import jsonschema
import pandas as pd
import requests

def parse_args():
    parser = argparse.ArgumentParser(description='Read criTRia curations from a JSON file and update them applying curation rules.')
    parser.add_argument('--json', default='data/criTRia-curations.json', 
                        help='Input JSON file path (default: data/criTRia-curations.json)')
    parser.add_argument('--out', default='data/criTRia-curations.json', 
                        help='Output JSON file path that will be overwritten if it exists (default: data/criTRia-curations.json)')
    parser.add_argument('--schema', default='data/criTRia-curations.schema.json', 
                        help='Optional JSON schema file path for validating the input JSON file (default: data/criTRia-curations.schema.json)')
    return parser.parse_args()

SUPERCATEGORY_LOOKUP = {
    "Genetic Evidence": {"max_score": 12},
    "Experimental Evidence": {"max_score": 6},
}

EVIDENCE_LOOKUP = {
    "Probands": {
        "evidence_category": "Singular Evidence",
        "evidence_supercategory": "Genetic Evidence",
        "evidence_max_score": 6,
        "category_max_score": 6,
    },
    "Allele": {
        "evidence_category": "Collective Evidence",
        "evidence_supercategory": "Genetic Evidence",
        "evidence_max_score": 2,
        "category_max_score": 3,
    },
    "Computational": {
        "evidence_category": "Collective Evidence",
        "evidence_supercategory": "Genetic Evidence",
        "evidence_max_score": 3,
        "category_max_score": 3,
    },
    "Segregation": {
        "evidence_category": "Collective Evidence",
        "evidence_supercategory": "Genetic Evidence",
        "evidence_max_score": 3,
        "category_max_score": 3,
    },
    "Case-control data": {
        "evidence_category": "Statistics",
        "evidence_supercategory": "Genetic Evidence",
        "evidence_max_score": 12,
        "category_max_score": 12,
    },
    "Biochemical function": {
        "evidence_category": "Function",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Protein interaction": {
        "evidence_category": "Function",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Regulatory impact": {
        "evidence_category": "Function",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Patient cells": {
        "evidence_category": "Functional Alteration",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Non-patient cells": {
        "evidence_category": "Functional Alteration",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 1,
        "category_max_score": 2,
    },
    "Non-human model organism": {
        "evidence_category": "Models",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 4,
        "category_max_score": 4,
    },
    "Cell culture": {
        "evidence_category": "Models",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Human treatment": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 4,
        "category_max_score": 4,
    },
    "Rescue in non-human model organism": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 4,
        "category_max_score": 4,
    },
    "Rescue in cell culture": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Rescue in patient cells": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
}

CLASSIFICATION_LOOKUP = {
    "Definitive": {
        "score_range": (12, 18),
        "replication_requirement": "≥2 publications, ≥3 years apart",
        "min_pubs": 2,
        "min_years": 3,
    },
    "Strong": {
        "score_range": (12, 18),
        "replication_requirement": "Compelling, not repeated over time",
        "min_pubs": 1,
        "min_years": 0,
    },
    "Moderate": {
        "score_range": (7, 11),
        "replication_requirement": "Some support, no contradictions",
        "min_pubs": 1,
        "min_years": 0,
    },
    "Limited": {
        "score_range": (0.1, 6),
        "replication_requirement": "Some evidence, not compelling",
        "min_pubs": 1,
        "min_years": 0,
    },
    "Disputed": {
        "score_range": None,
        "replication_requirement": "Contradictory evidence present",
        "min_pubs": 1,
        "min_years": 0,
    },
    "Refuted": {
        "score_range": None,
        "replication_requirement": "Causality ruled out",
        "min_pubs": 1,
        "min_years": 0,
    },
    "No Known Relationship": {
        "score_range": (-1, 0),
        "replication_requirement": "No evidence",
        "min_pubs": 0,
        "min_years": 0,
    },
}

def _normalize_evidence_type(evidence_type):
    if evidence_type is None or pd.isna(evidence_type):
        return ""
    return " ".join(str(evidence_type).strip().lower().split())

NORMALIZED_EVIDENCE_LOOKUP = {
    _normalize_evidence_type(evidence_type): mapping
    for evidence_type, mapping in EVIDENCE_LOOKUP.items()
}

CATEGORY_ORDER = []
for evidence_type in EVIDENCE_LOOKUP:
    category = EVIDENCE_LOOKUP[evidence_type]['evidence_category']
    if category not in CATEGORY_ORDER:
        CATEGORY_ORDER.append(category)

def lookup_evidence(evidence_type):
    return NORMALIZED_EVIDENCE_LOOKUP.get(_normalize_evidence_type(evidence_type), {
        "evidence_category": None,
        "evidence_supercategory": None,
    })

def fetch_pubmed_publication_date(pmid, max_retries=5):
    """Fetch the publication date for a given PubMed ID using the NCBI E-utilities API.
    Retries with exponential backoff on 429 Too Many Requests responses."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    # NCBI allows ~3 requests/second without an API key; add a small base delay.
    time.sleep(0.34)
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 429:
                wait = 2 ** attempt  # 1, 2, 4, 8, 16 seconds
                sys.stderr.write(f"Warning: Rate limited by NCBI for PMID {pmid}, retrying in {wait}s (attempt {attempt + 1}/{max_retries})\n")
                time.sleep(wait)
                continue
            response.raise_for_status()
            data = response.json()
        except requests.exceptions.RequestException as e:
            sys.stderr.write(f"Warning: Unable to fetch publication date for PMID {pmid}: {e}\n")
            return None
        except ValueError as e:
            sys.stderr.write(f"Warning: Unable to decode PubMed response for PMID {pmid}: {e}\n")
            return None

        if 'result' in data and pmid in data['result']:
            pub_date = data['result'][pmid].get('pubdate')
            return pub_date
        return None

    sys.stderr.write(f"Warning: Unable to fetch publication date for PMID {pmid} after {max_retries} retries (rate limited)\n")
    return None

def get_publication_dates(citation):
    """Extract all publication dates from a citation field containing one or more citations.
    Returns a list of ISO date strings, one for each resolvable PMID found."""
    if not isinstance(citation, str):
        return []

    publication_dates = []
    for citation_part in citation.split(';'):
        citation_part = citation_part.strip()
        if not citation_part or not citation_part.lower().startswith('pmid:'):
            continue

        pmid = citation_part.split(':', 1)[1].strip()
        if not pmid:
            continue

        pub_date = fetch_pubmed_publication_date(pmid)
        if pub_date:
            try:
                publication_dates.append(pd.to_datetime(pub_date).date().isoformat())
            except Exception as e:
                sys.stderr.write(f"Warning: Unable to parse publication date '{pub_date}' for PMID {pmid}: {e}\n")

    return publication_dates

def publication_interval(pub_dates):
    """Calculate the interval in years between the earliest and latest publication dates for a list of PubMed IDs."""
    if len(pub_dates) < 2:
        return None
    
    pub_dates = sorted(pd.to_datetime(d) for d in pub_dates if d is not None)

    if len(pub_dates) < 2:
        return None

    earliest_date = pd.to_datetime(pub_dates[0])
    latest_date = pd.to_datetime(pub_dates[-1])
    interval_years = (latest_date - earliest_date).days / 365.25
    return interval_years

def summarize_curations(locus):
    """ 
    Expecting curations to be in the the following format (only relevent fields shown):
    {
        "genetic_evidence_details": [
            {
                "Evidence type": "Probands",
                "Score": 6.0,
                "Citation": "pmid:39068203",
            }
        ],
        "experimental_evidence_details": [
            {
                "Evidence type": "Regulatory Impact",
                "Score": 0.5,
                "Citation": "pmid:39068203",
            }
        ]
    }
    The function should summarize the evidence for each locus and add/update the following values:
    "category_summary":
        {
            "Collective Evidence": 2.5,
            "Function": 1.0,
            "Functional Alteration": 1.5,
            "Singular Evidence": 6.0,
            "Statistics": 0,
            "Models": 0,
            "Rescue": 0
        },
        "supercategory_summary":
        {
            "Experimental Evidence": 2.5,
            "Genetic Evidence": 8.5
        },
        "total_score": 11.0,
        "classification": "Moderate",
        "publication_count": 1,
        "publication_interval_years": null,
    """
    # Add publication year for each citation
    all_pub_dates = []
    for evidence_entry in locus.get('genetic_evidence_details', []) + locus.get('experimental_evidence_details', []):
        citations = evidence_entry.get('Citation', None)
        publication_dates = evidence_entry.get('publication_dates', [])
        if citations is None:
            publication_dates = []
        elif len(publication_dates) < citations.count(';') + 1:
            publication_dates = get_publication_dates(citations)
        evidence_entry['publication_dates'] = publication_dates
        all_pub_dates.extend(publication_dates)

    locus_id = locus.get('Locus_ID', 'Unknown Locus')
    df = pd.DataFrame(locus.get('genetic_evidence_details', []) + locus.get('experimental_evidence_details', []))

    #sys.stdout.write(f"Processing locus '{locus_id}' with {len(df)} evidence entries...\n")
    if len(df) == 0:
        sys.stderr.write(f"Warning: No evidence entries found for locus '{locus_id}'\n")
        return locus

    # Report any evidence types that are not recognized in the EVIDENCE_LOOKUP
    for evidence_type in df['Evidence type'].unique():
        if _normalize_evidence_type(evidence_type) not in NORMALIZED_EVIDENCE_LOOKUP:
            sys.stderr.write(f"Warning: Unrecognized evidence type '{evidence_type}' in locus '{locus_id}'\n")


    # Sum the scores for each evidence category and supercategory
    category_scores = df.groupby('evidence_category')['Score'].sum().to_dict()
    category_summary = {
        category: category_scores.get(category, 0)
        for category in CATEGORY_ORDER
    }

    # Enforce the maximum score for each category        
    for category, score in category_summary.items():
        max_score = max(EVIDENCE_LOOKUP[evidence_type]['category_max_score'] for evidence_type in EVIDENCE_LOOKUP if EVIDENCE_LOOKUP[evidence_type]['evidence_category'] == category)
        if score is not None and max_score is not None:
            if score > max_score:
                category_summary[category] = max_score

    # Sum the scores for each supercategory
    supercategory_summary = df.groupby('evidence_supercategory')['Score'].sum().to_dict()
    # If there are any supercategories that are missing from the summary, add them with a score of 0
    for supercategory in SUPERCATEGORY_LOOKUP:
        if supercategory not in supercategory_summary:
            supercategory_summary[supercategory] = 0
    # Enforce the maximum score for each supercategory
    for supercategory, score in supercategory_summary.items():
        max_score = SUPERCATEGORY_LOOKUP[supercategory]['max_score']
        if score is not None and max_score is not None:
            if score > max_score:
                supercategory_summary[supercategory] = max_score

    # Calculate the total score for the curation as the sum of the supercategory scores
    total_score = sum(supercategory_summary.values())

    def count_unique_publications(citations):
        unique_publications = set()
        for citation in citations.dropna():
            for publication in str(citation).split(';'):
                normalized_publication = publication.strip().lower()
                if normalized_publication:
                    unique_publications.add(normalized_publication)
        return len(unique_publications)

    # Publication count
    publication_count = count_unique_publications(df['Citation'])
    publication_interval_years = publication_interval(list(set(all_pub_dates)))
    if publication_interval_years is not None:
        publication_interval_years = round(publication_interval_years, 2)
    
    # Determine the classification based on the total score and replication requirement
    classification = None
    # if  publication_interval_years is None, set it to 0 for the purpose of classification, but we will still report it as None in the final summary
    if publication_interval_years is None:
        tmp_publication_interval_years = 0
    else:
        tmp_publication_interval_years = publication_interval_years

    # Round total_score to nearest integer for classification lookup
    total_score_rounded = round(total_score)
    for class_name, class_info in CLASSIFICATION_LOOKUP.items():
        score_range = class_info['score_range']
        if score_range is not None and score_range[0] <= total_score_rounded <= score_range[1]:
            if publication_count >= class_info['min_pubs'] and tmp_publication_interval_years >= class_info['min_years']:
                classification = class_name
                break
    # If there's a value in Manual_evidence_level, use that as the classification instead of the calculated one
    if 'Manual_evidence_level' in locus and locus['Manual_evidence_level']:
        classification = locus['Manual_evidence_level']

    if publication_interval_years == 0:
        publication_interval_years = None

    # Collect the Genetic and Experimental evidence details for the curation
    genetic_evidence_details = df[df['evidence_supercategory'] == 'Genetic Evidence'].to_dict(orient='records')
    experimental_evidence_details = df[df['evidence_supercategory'] == 'Experimental Evidence'].to_dict(orient='records')

    # Return updated dictionary containing the curation summary
    locus["category_summary"] = category_summary
    locus["supercategory_summary"] = supercategory_summary
    locus["total_score"] = total_score
    locus["classification"] = classification
    locus["publication_count"] = publication_count
    locus["publication_interval_years"] = publication_interval_years
    locus["genetic_evidence_details"] = genetic_evidence_details
    locus["experimental_evidence_details"] = experimental_evidence_details
    return locus

def sanitize_for_json(value):
    if isinstance(value, dict):
        return {key: sanitize_for_json(val) for key, val in value.items()}
    if isinstance(value, list):
        return [sanitize_for_json(item) for item in value]
    if isinstance(value, tuple):
        return [sanitize_for_json(item) for item in value]
    if isinstance(value, pd.Timestamp):
        return value.date().isoformat()
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if pd.isna(value):
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            return sanitize_for_json(value.item())
        except (ValueError, TypeError):
            pass
    return value

def main(args):

    if args.out == args.json:
        pause = 5
        sys.stderr.write(f'WARNING: overwriting {args.json} in {pause} seconds\n')
        sys.stderr.write('Press Ctrl+C to cancel\n')
        time.sleep(pause)
    else:
        sys.stderr.write(f'Writing {args.out}\n')

    out_data = []
    with open(args.json, 'r') as in_json:
        data = json.load(in_json)

        for locus in data:
            out_data.append(summarize_curations(locus))

    out_data = sanitize_for_json(out_data)

    # Sort records by gene name then id
    out_data = sorted(out_data, key = lambda x: (x['Gene'], x['Locus_ID']))

    # Make sure json is sorted within each record based on the schema order
    if args.schema:
        with open(args.schema, 'r') as schema_file:
            schema = json.load(schema_file)
            if 'properties' in schema:
                schema_properties = list(schema['properties'].keys())
            elif isinstance(schema.get('items'), dict) and 'properties' in schema['items']:
                schema_properties = list(schema['items']['properties'].keys())
            else:
                raise KeyError("Schema must contain 'properties' at the top level or under 'items'")
            out_data = [
                {key: record.get(key) for key in schema_properties} for record in out_data
            ]

        # Validate the processed data against the schema
        validator = jsonschema.Draft202012Validator(schema)
        any_errors = False
        for error in sorted(validator.iter_errors(out_data), key=str):
            any_errors = True
            sys.stderr.write(f"Schema validation error: {error.message}\n")
            sys.stderr.write(f"  Path: {list(error.absolute_path)}\n")
        if any_errors:
            sys.stderr.write("Schema validation failed.\n")
            sys.exit(1)
        else:
            sys.stderr.write("Schema validation succeeded.\n")

    # Write JSON file
    with open(args.out, 'w') as out_json_file:
        options = jsbeautifier.default_options()
        options.indent_size = 2
        options.brace_style="expand"
        out_json_file.write(jsbeautifier.beautify(json.dumps(out_data, ensure_ascii=False, allow_nan=False, default=str), options))
        out_json_file.write('\n')
        
if __name__ == '__main__':
    doctest.testmod()
    args = parse_args()
    main(args)
