# Extract criTRia curations from a tsv file and convert them to a json format that can be used by the criTRia web applicatio
import sys
import argparse
import doctest
import os
import time
import json
import math
import jsbeautifier
import pandas as pd
import requests

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
    "Case-Control Data": {
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
    "Protein Interaction": {
        "evidence_category": "Function",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Regulatory Impact": {
        "evidence_category": "Function",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Patient Cells": {
        "evidence_category": "Functional Alteration",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Non-patient Cells": {
        "evidence_category": "Functional Alteration",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 1,
        "category_max_score": 2,
    },
    "Non-Human Model Organism": {
        "evidence_category": "Models",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 4,
        "category_max_score": 4,
    },
    "Cell Culture": {
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
    "Rescue in Non-Human Model Organism": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 4,
        "category_max_score": 4,
    },
    "Rescue in Cell Culture": {
        "evidence_category": "Rescue",
        "evidence_supercategory": "Experimental Evidence",
        "evidence_max_score": 2,
        "category_max_score": 2,
    },
    "Rescue in Patient Cells": {
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
        "score_range": (6, 12),
        "replication_requirement": "Some support, no contradicitions",
        "min_pubs": 1,
        "min_years": 0,
    },
    "Limited": {
        "score_range": (0, 6),
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


def lookup_evidence(evidence_type):
    return NORMALIZED_EVIDENCE_LOOKUP.get(_normalize_evidence_type(evidence_type), {
        "evidence_category": None,
        "evidence_supercategory": None,
    })

def parse_args():
    parser = argparse.ArgumentParser(description='Extract criTRia curations from TSV files and convert them to JSON format.')
    parser.add_argument('tsv', nargs='+', help='Input TSV file(s) containing criTRia curations')
    parser.add_argument('--json', default='data/criTRia-curations.json', help='Output JSON file path (default: data/criTRia-curations.json)')
    return parser.parse_args()

def extract_curations_from_tsv(tsv_file):
    """ Expected format of the TSV file:
#Disease_ID: FECD3	Gene: TCF4	Locus_ID: FECD3_TCF4	Inheritance: AD	Curator: Macayla Weiner	Date: 3/2/26
Evidence type	Score	Citation	Values	Evidence detail	Notes
Probands	6	pmid:25168903	probands:46		68 affected individuals + 1 unaffected with expansion but 46 affected with expansion																
    """
    # Extract the metadata from the first line of the TSV file
    sys.stderr.write(f"Extracting curations from {tsv_file}...\n")
    with open(tsv_file, 'r') as f:
        first_line = f.readline().strip('#').strip()
        second_line = f.readline().strip()
    metadata = {}
    for item in first_line.split('\t'):
        key, value = item.split(':', 1)
        metadata[key.strip()] = value.strip()

    # skip cols with missing headers
    header_indices = [i for i, item in enumerate(second_line.split('\t')) if item.strip() != '']

    df = pd.read_csv(tsv_file, sep='\t', comment='#', usecols=header_indices)
    source_columns = df.columns.tolist()
    df = df.replace(r'^\s*$', pd.NA, regex=True)
    df = df.dropna(subset=source_columns, how='all')

    if 'Score' in df.columns:
        invalid_score_mask = df['Score'].notna() & pd.to_numeric(df['Score'], errors='coerce').isna()
        if invalid_score_mask.any():
            invalid_scores = sorted(set(df.loc[invalid_score_mask, 'Score'].astype(str).tolist()))
            sys.stderr.write(
                f"Warning: Non-numeric Score value(s) in {tsv_file}: {', '.join(invalid_scores)}. Treating as missing.\n"
            )
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')

    # use LOOKUP_EVIDENCE to add evidence_category and evidence_supercategory columns to the dataframe
    df['evidence_category'] = df['Evidence type'].apply(lambda x: lookup_evidence(x)['evidence_category'])
    df['evidence_supercategory'] = df['Evidence type'].apply(lambda x: lookup_evidence(x)['evidence_supercategory'])
    df['evidence_max_score'] = df['Evidence type'].apply(lambda x: lookup_evidence(x).get('evidence_max_score'))
    df['category_max_score'] = df['Evidence type'].apply(lambda x: lookup_evidence(x).get('category_max_score'))

    return metadata, df

def fetch_pubmed_publication_date(pmid):
    """Fetch the publication date for a given PubMed ID using the NCBI E-utilities API."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'result' in data and pmid in data['result']:
            pub_date = data['result'][pmid].get('pubdate')
            return pub_date
    return None

def publication_interval(pmids):
    """Calculate the interval in years between the earliest and latest publication dates for a list of PubMed IDs."""
    pmids = [pmid.replace('pmid:', '') for pmid in pmids if isinstance(pmid, str) and pmid.startswith('pmid:')]
    pub_dates = []
    for pmid in pmids:
        pub_date = fetch_pubmed_publication_date(pmid)
        if pub_date:
            pub_dates.append(pub_date)
    if len(pub_dates) < 2:
        return None
    pub_dates.sort()
    earliest_date = pd.to_datetime(pub_dates[0])
    latest_date = pd.to_datetime(pub_dates[-1])
    interval_years = (latest_date - earliest_date).days / 365.25
    return interval_years

def summarize_curations(locus_id, curations):
    """ 
    Expecting curations to be a list of tuples (metadata, df) where metadata is a dictionary containing the metadata for the curation and 
    df is a pandas dataframe containing the evidence for the curation. 
    The function should summarize the evidence for each locus and return a dictionary containing the summary.
    """
    for metadata, df in curations:
        # Report any evidence types that are not recognized in the EVIDENCE_LOOKUP
        for evidence_type in df['Evidence type'].unique():
            if _normalize_evidence_type(evidence_type) not in NORMALIZED_EVIDENCE_LOOKUP:
                sys.stderr.write(f"Warning: Unrecognized evidence type '{evidence_type}' in locus '{locus_id}'\n")


        # Sum the scores for each evidence category and supercategory
        category_summary = df.groupby('evidence_category')['Score'].sum().to_dict()
        # If there are any categories that are missing from the summary, add them with a score of 0
        for category in set(EVIDENCE_LOOKUP[evidence_type]['evidence_category'] for evidence_type in EVIDENCE_LOOKUP):
            if category not in category_summary:
                category_summary[category] = 0

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

        # Publication count
        publication_count = df['Citation'].nunique()
        publication_interval_years = publication_interval(df['Citation'].dropna().unique())
        if publication_interval_years is not None:
            publication_interval_years = round(publication_interval_years, 2)
        
        # Determine the classification based on the total score and replication requirement
        classification = None
        # if  publication_interval_years is None, set it to 0 for the purpose of classification, but we will still report it as None in the final summary
        if publication_interval_years is None:
            tmp_publication_interval_years = 0
        else:
            tmp_publication_interval_years = publication_interval_years

        for class_name, class_info in CLASSIFICATION_LOOKUP.items():
            score_range = class_info['score_range']
            if score_range is not None and score_range[0] < total_score <= score_range[1]:
                if publication_count >= class_info['min_pubs'] and tmp_publication_interval_years >= class_info['min_years']:
                    classification = class_name
                    break
        # If there's a value in Manual_evidence_level, use that as the classification instead of the calculated one
        if 'Manual_evidence_level' in metadata and metadata['Manual_evidence_level']:
            classification = metadata['Manual_evidence_level']

        # Collect the Genetic and Experimental evidence details for the curation
        genetic_evidence_details = df[df['evidence_supercategory'] == 'Genetic Evidence'].to_dict(orient='records')
        experimental_evidence_details = df[df['evidence_supercategory'] == 'Experimental Evidence'].to_dict(orient='records')

        # Return a dictionary containing the summary for the curation
        return {
            **metadata,
            "category_summary": category_summary,
            "supercategory_summary": supercategory_summary,
            "total_score": total_score,
            "classification": classification,
            "publication_count": publication_count,
            "publication_interval_years": publication_interval_years,
            "genetic_evidence_details": genetic_evidence_details,
            "experimental_evidence_details": experimental_evidence_details
        }


def sanitize_for_json(value):
    if isinstance(value, dict):
        return {key: sanitize_for_json(val) for key, val in value.items()}
    if isinstance(value, list):
        return [sanitize_for_json(item) for item in value]
    if isinstance(value, tuple):
        return [sanitize_for_json(item) for item in value]
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
    raw_data = {}
    for tsv_file in args.tsv:
        metadata, df = extract_curations_from_tsv(tsv_file)
        if metadata['Locus_ID'] not in raw_data:
            raw_data[metadata['Locus_ID']] = []
        raw_data[metadata['Locus_ID']].append((metadata, df))

    # Summarize the data for each locus and convert it to the format expected by the criTRia web application
    data = []
    for locus_id, curations in raw_data.items():
        summary = summarize_curations(locus_id, curations)
        data.append(summary)

    data = sanitize_for_json(data)

    # Write JSON file
    with open(args.json, 'w') as out_json_file:
        options = jsbeautifier.default_options()
        options.indent_size = 2
        options.brace_style="expand"
        out_json_file.write(jsbeautifier.beautify(json.dumps(data, ensure_ascii=False, allow_nan=False), options))
        out_json_file.write('\n')
        

if __name__ == '__main__':
    doctest.testmod()
    args = parse_args()
    main(args)
