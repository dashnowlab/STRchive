# Validate the STRchive loci JSON file and overwrite it with changes if needed.

import sys
import argparse
import doctest
import os
import time
import json
from Bio.Seq import Seq
import pandas as pd
import re
#import numpy as np

# in case person runs from data instead of scripts
sys.path.append('../')

# This should be overwritten if schema json file provided
list_fields = [
    "reference_motif_reference_orientation",
    "pathogenic_motif_reference_orientation",
    "pathogenic_motif_gene_orientation",
    "benign_motif_reference_orientation",
    "benign_motif_gene_orientation",
    "unknown_motif_reference_orientation",
    "unknown_motif_gene_orientation",
    "inheritance",
    "source",
    "omim",
    "stripy",
    "gnomad",
    "genereviews",
    "mondo",
    "medgen",
    "orphanet",
    "gard",
    "malacard",
    "webstr_hg38",
    "webstr_hg19",
    "tr_atlas",
    "locus_tags",
    "disease_tags",
]

def parse_args():
    parser = argparse.ArgumentParser(description='Validate the STRchive loci JSON file and overwrite it with changes if needed.')
    parser.add_argument('json', default='data/STRchive-database.json', help='STRchive loci JSON input file path')
    parser.add_argument('--schema', default=None, help='JSON schema file path (optional)')
    parser.add_argument('--out', default=None, help='Defaults to same as input JSON file path (overwriting)')
    parser.add_argument('--pause', default=5, help='Pause time in seconds before overwriting the file. Default: 5')
    return parser.parse_args()

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All circular permutations of x
    >>> circular_permuted('GAG')
    ['GAG', 'AGG', 'GGA']
    >>> circular_permuted('AG')
    ['AG', 'GA']
    >>> circular_permuted('TAGAA')
    ['TAGAA', 'AGAAT', 'GAATA', 'AATAG', 'ATAGA']
    >>> circular_permuted('')
    []
    """
    n = len(x)
    modified_sequences = []
    modified_sequences.extend([x[i:] + x[:i] for i in range(n)])
    return modified_sequences

def normalise_str(in_dna):
    """
    Args:
        in_dna (sequence)
    Returns:
        the normalized output of the string
    Find all possible equivalent STR sequences and return the first alphabetically for each
    >>> normalise_str('ATAG')
    'AGAT'
    >>> normalise_str('NGC')
    'CNG'
    >>> normalise_str('AAG')
    'AAG'
    """
    if in_dna is None or len(in_dna) == 0:
        return ''

    all_possible = []
    for permuted_seq in circular_permuted(in_dna):
        all_possible.append(permuted_seq)

    return min(all_possible)

def get_new_motif(motif, gene_strand):
    """
    Args:
        motif (string)
        gene_strand: either + or -
    Returns:
        the normalized output of the string from ref to gene orientation
    Get the new normalized motif for each row.
    If gene_strand is +, reference orientation = gene orientation
    If gene_strand is -, reverse_complement ref_ori for gene_ori
    >>> get_new_motif('GAG', '+')
    'AGG'
    >>> get_new_motif('GAG', '-')
    'CCT'
    >>> get_new_motif('TCATC', '-')
    'AGATG'
    >>> get_new_motif('TAG', 'plus')
    Traceback (most recent call last):
    ...
    AssertionError: Gene strand plus is not +/-
    """
    if gene_strand == "+":
        normalized_motif = normalise_str(motif)
    elif gene_strand == "-":
        seq = Seq(motif)
        reverse_comp = str(seq.reverse_complement())
        normalized_motif = normalise_str(reverse_comp)
    else:
        raise AssertionError(f'Gene strand {gene_strand} is not +/-')
    return normalized_motif
    
def check_motif_orientation(record):
    """
    Args:
        record (dict): a dictionary containing a single locus from the STRchive json
    Returns:
        record (dict): the record with any motif fields with incorrect orientation updated
    """
    field_pairs = [
        ('pathogenic_motif_reference_orientation', 'pathogenic_motif_gene_orientation'),
        ('benign_motif_reference_orientation', 'benign_motif_gene_orientation'),
        ('unknown_motif_reference_orientation', 'unknown_motif_gene_orientation')
    ]
    for ref_field, gene_field in field_pairs:
        if record[ref_field] is None:
            continue
        old = record[gene_field]
        new = [get_new_motif(x, record['gene_strand']) for x in record[ref_field]]
        if old != new:
            for old_motif, new_motif in zip(old, new):
                if old_motif != new_motif:
                    sys.stderr.write(f'Updating {record['id']} {gene_field} from {old_motif} to {new_motif}\n')
        record[gene_field] = new
    return record

def check_list_fields(record, list_fields = list_fields):
    """
    Args:
        record (dict): a dictionary containing a single locus from the STRchive json
        list_fields (list): a list of names of fields that should be lists
    Returns:
        record (dict): the record with any fields that should be lists converted to lists
    """
    for field in list_fields:
        old = record[field]
        if record[field] == "" or record[field] is None or record[field] == []:
            record[field] = [] #XXX need to decide if Null or empty list is preferred
        elif isinstance(record[field], str):
            if ";" in record[field]:
                record[field] = [x.strip() for x in record[field].split(';')]
            elif "," in record[field]:
                record[field] = [x.strip() for x in record[field].split(',')]
            else:
                record[field] = [record[field].strip()]
        if old != record[field]:
            sys.stderr.write(f'Updating {record['id']} {field} from {old} to {record[field]}\n')
    return record

def main(json_fname, json_schema = None, out_json = None, pause = 5):
    sys.stderr.write('Reading JSON and overwriting it with fixed version\n')
    if out_json == json_fname:
        sys.stderr.write(f'WARNING: overwriting {json_fname} in {pause} seconds\n')
        sys.stderr.write('Press Ctrl+C to cancel\n')
        time.sleep(pause)
    else:
        sys.stderr.write(f'Writing {out_json}\n')

    # Check if file exists
    if not os.path.exists(json_fname):
        sys.stderr.write(f'Error: {json_fname} does not exist\n')
        sys.exit(1)
    # Read JSON file
    with open(json_fname, 'r') as json_file:
        data = json.load(json_file)

        # Check if the field contains a string that should be a list
        for record in data:
            record = check_list_fields(record)
            record = check_motif_orientation(record)
        
        # Write JSON file
        with open(out_json, 'w') as out_json_file:
            json.dump(data, out_json_file, indent=4, separators=(',', ':'))

if __name__ == '__main__':
    doctest.testmod()
    args = parse_args()
    if args.out is None:
        args.out = args.json
    main(args.json, args.schema, args.out, args.pause)
