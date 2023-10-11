from Bio.Seq import Seq
import pandas as pd
from itertools import product
import re
import numpy as np
import math

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All circular permutations of x, with all 'Ns' replaced by 'A', 'T', 'G', and 'C'
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
    Find all possible equivalent STR sequences and return the first alphabetically for each, with
    multiple possible outputs for each replacement of 'N'.
    """
    if in_dna is None or len(in_dna) == 0:
        return ''

    all_possible = []

    # # Find all positions of 'N' in the input sequence
    # n_positions = [i for i, nucleotide in enumerate(in_dna) if nucleotide == 'N']

    # # If 'N' is present, replace each 'N' with 'A', 'T', 'G', and 'C' separately
    # if n_positions:
    #     results = []
    #     for replacement_combination in product('ATGC', repeat=len(n_positions)):
    #         modified_dna = list(in_dna)
    #         for i, replacement in zip(n_positions, replacement_combination):
    #             modified_dna[i] = replacement
    #         modified_results = circular_permuted("".join(modified_dna))
    #         if modified_results:
    #             results.append(min(modified_results))
    #     return results
    # else:
        # Circularly permute the original sequence and reverse complement
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
        benign_motifs = [motif.strip() for motif in re.split(r',',
                                                    benign_reference_orientation) if motif.strip()]
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

    #update the dataframe
    df['pathogenic_motif_gene_orientation'] = [', '.join(x) for x in pathogenic_results]
    df['benign_motif_gene_orientation'] = [', '.join(x) for x in benign_results]
    df['unknown_motif_gene_orientation'] = [', '.join(x) for x in unknown_results]

    # Save the updated
    df.to_csv(out_csv, index=False)
