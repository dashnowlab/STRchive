from Bio.Seq import Seq
import pandas as pd
from itertools import product
import re
import numpy as np

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All circular permutations of x, with all 'Ns' replaced by 'A', 'T', 'G', and 'C'
    """
    n = len(x)
    modified_sequences = []

    # Replace 'N' with 'A', 'T', 'G', and 'C' separately and generate circular permutations for each
    for replacement in ['A', 'T', 'G', 'C']:
        modified_x = x.replace('N', replacement)
        modified_sequences.extend([modified_x[i:] + modified_x[:i] for i in range(n)])

    return modified_sequences

def normalise_str(in_dna):
    """
    Args:
        in_dna (sequence)
    Returns:
        the normalized output of the string
    Find all possible equivalent STR sequences and return the first alphabetically for each, with multiple possible
    outputs for each replacement of 'N'.
    """
    if in_dna is None or len(in_dna) == 0:
        return ''

    all_possible = []

    # Find all positions of 'N' in the input sequence
    n_positions = [i for i, nucleotide in enumerate(in_dna) if nucleotide == 'N']

    # If 'N' is present, replace each 'N' with 'A', 'T', 'G', and 'C' separately
    if n_positions:
        results = []
        for replacement_combination in product('ATGC', repeat=len(n_positions)):
            modified_dna = list(in_dna)
            for i, replacement in zip(n_positions, replacement_combination):
                modified_dna[i] = replacement
            modified_results = circular_permuted("".join(modified_dna))
            if modified_results:
                results.append(min(modified_results))
        return results
    else:
        # Circularly permute the original sequence and reverse complement
        for permuted_seq in circular_permuted(in_dna):
            all_possible.append(permuted_seq)

    return [min(all_possible)]

def process_csv(in_csv, out_csv):

    # Read the CSV file into a DataFrame
    df = pd.read_csv(in_csv, dtype=str)

    # Create an empty list to store the results
    results = []

    # Iterate through the rows of the DataFrame
    for index, row in df.iterrows():
        gene_strand = row['gene_strand']
        reference_orientation = row['pathogenic_motif_reference_orientation']

        # Check if gene_strand is "+"
        if gene_strand == "+":
            pathogenic_motif_gene_orientation = reference_orientation
        # Deal with the commas
        elif ',' in reference_orientation:
            normalized_reference_rc = []
            reference_motifs = [motif.strip() for motif in re.split(r',', reference_orientation)]
            # print("reference motifs", reference_motifs)
            # get the reverse complement of the reference orientation pathogenic motifs, and normalize
            for motif in reference_motifs:
                seq = Seq(motif)
                reverse_comp = str(seq.reverse_complement())
                normalized_reference_motifs_rc = normalise_str(reverse_comp)
                normalized_reference_rc.append(str(normalized_reference_motifs_rc))  # Append reverse complement as a string to the list
            pathogenic_motif_gene_orientation = normalized_reference_rc
        else:
            # Reverse complement for gene_strand == "-"
            reference_orientation_rc = str(Seq(reference_orientation).reverse_complement())

            # Run normalise_str for reference_orientation and its reverse complement
            normalized_reference_rc = normalise_str(reference_orientation_rc)

            pathogenic_motif_gene_orientation = normalized_reference_rc

        # Append the result to the results list
        results.append(pathogenic_motif_gene_orientation)

    df['pathogenic_motif_gene_orientation'] = results

    # Apply the normalization function to the filtered rows
    # print(df['pathogenic_motif_gene_orientation'])

    # df['pathogenic_motif_gene_orientation'] = df['pathogenic_motif_gene_orientation'].apply(
    #     lambda x: ','.join(map(str, [normalise_str(m) for m in x.split(',')])) if ',' in x else normalise_str(x)
    # )

    # ### Pathogenic
    #     df['pathogenic_motif_gene_orientation'] = df['pathogenic_motif_reference_orientation'].apply(lambda x: str(Seq(x).reverse_complement()))
    # # Standardize the new column, handling multiple values separated by commas
    #     df['pathogenic_motif_gene_orientation'] = df['pathogenic_motif_reference_orientation'].apply(lambda x: ','.join(map(str, [normalise_str(m) for m in x.split(',')])) if ',' in x else normalise_str(x))

    ### Benign
    df['benign_motif_gene_orientation'] = df['benign_motif_reference_orientation'].apply(lambda x: str(Seq(x).reverse_complement()))
    # Standardize the new column, handling multiple values separated by commas
    df['benign_motif_gene_orientation'] = df['benign_motif_gene_orientation'].apply(lambda x: ','.join(map(str, [normalise_str(m) for m in x.split(',')])) if ',' in x else normalise_str(x))

    ### Unknown
    df['unknown_motif_reference_orientation'] = df['unknown_motif_reference_orientation'].apply(lambda x: str(x) if isinstance(x, str) else (str(x) if x is not None else None))
    # Apply transformations
    df['unknown_motif_gene_orientation'] = df['unknown_motif_reference_orientation'].apply(lambda x: str(Seq(x).reverse_complement()) if isinstance(x, str) else None)
    # Standardize the new column, handling multiple values separated by commas
    df['unknown_motif_gene_orientation'] = df['unknown_motif_gene_orientation'].apply(lambda x: ','.join(map(str, [normalise_str(m) for m in x.split(',')])) if isinstance(x, str) and ',' in x else (normalise_str(x) if isinstance(x, str) else None))

    # df = reverse_complement_and_standardize(df, df['benign_motif_reference_orientation'], df['benign_motif_gene_orientation'])
    df = df.applymap(lambda x: str(x).replace('[', '').replace(']', '').replace("'", '').replace('"',''))

    # Replace 'NNT' values with blank (empty) values in the entire DataFrame
    df = df.applymap(lambda x: '' if x == 'nnt' or x == 'nan' else x)
    # Save the updated DataFrame to a CSV file
    df.to_csv(out_csv, index=False)
