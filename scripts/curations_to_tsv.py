import json
import csv
import sys
import argparse
from pathlib import Path

def parse_args():
    """Parse command-line arguments."""
    repo_root = Path(__file__).parent.parent
    default_json = repo_root / 'data' / 'criTRia-curations.json'
    default_tsv = repo_root / 'data' / 'criTRia-curations.tsv'

    parser = argparse.ArgumentParser(
        description='Convert curations JSON to TSV format.'
    )
    parser.add_argument(
        '-i', '--input',
        default=default_json,
        type=Path,
        help=f'Path to input curations JSON file (default: {default_json})'
    )
    parser.add_argument(
        '-o', '--output',
        default=default_tsv,
        type=Path,
        help=f'Path to output TSV file (default: {default_tsv})'
    )
    return parser.parse_args()

def load_curations(json_file):
    """Load curations from JSON file."""
    with open(json_file, 'r') as f:
        return json.load(f)

def flatten_curation(curation):
    """Flatten a curation object for TSV output."""
    flat = {}
    
    # Direct fields
    for key in ['Disease_ID', 'Gene', 'Inheritance', 
                # 'Locus_ID', 'Curator', 'Source', 'SOP_version',
                'Date', 
                'Description', 
                'total_score', 'publication_count', 'publication_interval_years', 
                'classification']:
        flat[key] = curation.get(key)
    
    # Category and supercategory summaries
    # category_summary = curation.get('category_summary', {})
    supercategory_summary = curation.get('supercategory_summary', {})
    
    flat['genetic_evidence_score'] = supercategory_summary.get('Genetic Evidence')
    flat['experimental_evidence_score'] = supercategory_summary.get('Experimental Evidence')
    
    # for category in ['Collective Evidence', 'Function', 'Functional Alteration', 
    #                  'Singular Evidence', 'Statistics', 'Rescue', 'Models']:
    #     flat[f'category_{category}'] = category_summary.get(category)
    
    return flat

def curations_to_tsv(json_file, tsv_file):
    """Convert curations JSON to TSV format."""
    curations = load_curations(json_file)
    
    # Filter out incomplete entries (those with only classification)
    #valid_curations = [c for c in curations if 'Locus_ID' in c]
    # Filter to entries where Source is criTRia
    valid_curations = [c for c in curations if c.get('Source') == 'criTRia']
    
    if not valid_curations:
        print("No valid curations found.", file=sys.stderr)
        return
    
    # Flatten all curations
    flattened = [flatten_curation(c) for c in valid_curations]
    
    # Get all fieldnames
    fieldnames = sorted(set().union(*(d.keys() for d in flattened)))
    
    # Put key fields first
    key_fields = ['Gene', 'Disease_ID', 'Inheritance', 
                  'Date', 'classification', 
                  'total_score', 'genetic_evidence_score', 'experimental_evidence_score', 
                  'publication_count', 'publication_interval_years']
    fieldnames = key_fields + [f for f in fieldnames if f not in key_fields]

    # Write TSV
    with open(tsv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(flattened)
    
    print(f"Successfully wrote {len(flattened)} loci to {tsv_file}")

if __name__ == '__main__':
    args = parse_args()
    curations_to_tsv(args.input, args.output)