from Bio.Seq import Seq
import pandas as pd
from itertools import product
import re
import numpy as np
import math
import doctest
import sys

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

def compare_and_print(old_results, new_results, name):
    """
    Args:
        old_results: any values in the gene_orientation columns from input
        new_results: the gene_ori results based on ref_ori and gene_strand
        name: the specific result type (pathogenic, benign, unknown)
    Returns:
        print statements of differences for each result type, if any
    >>> old_results = ['A', 'C', 'G']
    >>> new_results = ['A', 'T', 'G']
    >>> compare_and_print(old_results, new_results, 'pathogenic_results')
    There are differences in pathogenic_results: [('C', 'T')]

    >>> old_results = [None, 'C', 'G']
    >>> new_results = ['A', 'C', None]
    >>> compare_and_print(old_results, new_results, 'benign_results')

    >>> old_results = ['A', 'C', 'G']
    >>> new_results = ['A', 'C', 'G']
    >>> compare_and_print(old_results, new_results, 'unknown_results')

    """
    differences = []
    for old_item, new_item in zip(old_results, new_results):
        if not pd.isna(old_item) and not pd.isna(new_item) and old_item != new_item:
            differences.append((old_item, new_item))

    if differences:
        sys.stderr.write(f"The following differences in {name} were replaced: {differences}\n")

def process_csv(in_csv, out_csv):
    """
    Args:
        in_csv: the input csv
    Returns:
        out_csv: the output csv, processed with the new columns
   Add the gene orientation columns with normalized motifs
    """
    # Read the CSV file into a DataFrame
    df = pd.read_csv(in_csv, dtype=str)

    # we need these empty bits to fill with the new column results!
    pathogenic_results = []
    benign_results = []
    unknown_results = []

    for index, row in df.iterrows():
        # input the rows we need: the gene strand and reference orientation
        gene_strand = row['gene_strand']
        pathogenic_reference_orientation = row['pathogenic_motif_reference_orientation']
        benign_reference_orientation = row['benign_motif_reference_orientation']
        unknown_reference_orientation = row['unknown_motif_reference_orientation']

        # because there can be multiple motifs, we need to check for commas/lists
        # and then we can get our new motifs
        pathogenic_motifs = [motif.strip() for motif in re.split(r',',
                                                    pathogenic_reference_orientation) if motif.strip()]
        normalized_pathogenic_motifs = [get_new_motif(motif, gene_strand) for motif in pathogenic_motifs]

        try:
            benign_motifs = [motif.strip() for motif in re.split(r',',
                                                    benign_reference_orientation) if motif.strip()]
        except TypeError as e:
            sys.stderr.write(f"Unexpected value in benign_reference_orientation: {benign_reference_orientation}\n")
            benign_motifs = []

        normalized_benign_motifs = [get_new_motif(motif, gene_strand) for motif in benign_motifs]

        # we don't always have unknown motifs, so we need to check the values before normalizing
        if isinstance(unknown_reference_orientation, str) and unknown_reference_orientation.strip():
            unknown_motifs = [motif.strip() for motif in re.split(r',', unknown_reference_orientation)]
            normalized_unknown_motifs = [get_new_motif(motif, gene_strand) for motif in unknown_motifs]
        else:
            normalized_unknown_motifs = []

        # Append the list of normalized motifs to results
        pathogenic_results.append(normalized_pathogenic_motifs)
        benign_results.append(normalized_benign_motifs)
        unknown_results.append(normalized_unknown_motifs)


    # Initialize old results lists, for comparison of input and output data
    old_pathogenic_results = []
    old_benign_results = []
    old_unknown_results = []

    # for each new column, let's compare if there are any prev values
    if 'pathogenic_motif_gene_orientation' in df.columns:
        old_pathogenic_results = df['pathogenic_motif_gene_orientation'].apply(lambda x: x.split(',') if isinstance(x, str) else None).tolist()
        compare_and_print(old_pathogenic_results, pathogenic_results, 'pathogenic_motif_gene_orientation')
        df['pathogenic_motif_gene_orientation'] = [','.join(x) for x in pathogenic_results]

    if 'benign_motif_gene_orientation' in df.columns:
        old_benign_results = df['benign_motif_gene_orientation'].apply(lambda x: x.split(',') if isinstance(x, str) else None).tolist()
        compare_and_print(old_benign_results, benign_results, 'benign_motif_gene_orientation')
        df['benign_motif_gene_orientation'] = [','.join(x) for x in benign_results]

    if 'unknown_motif_gene_orientation' in df.columns:
        old_unknown_results = df['unknown_motif_gene_orientation'].apply(lambda x: x.split(',') if isinstance(x, str) else None).tolist()
        compare_and_print(old_unknown_results, unknown_results, 'unknown_motif_gene_orientation')
        df['unknown_motif_gene_orientation'] = [','.join(x) for x in unknown_results]

    # Save the updated
    df.to_csv(out_csv, index=False)

if __name__ == "__main__":
    doctest.testmod()