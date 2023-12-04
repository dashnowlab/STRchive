import csv
import json
import sys
import re

def chrom_to_int(chrom):
    """
    Convert to int, such that order should be chr 1-22, X, Y, M
    :param chrom: chromosome name
    :return: chromosome number
    """
    chrom = chrom.strip('chr')
    if chrom == 'X':
        return 100
    elif chrom == 'Y':
        return 101
    elif chrom == 'M':
        return 102
    else:
        return int(chrom)
    
def clean_loci(data, genome):
    """Exclude incomplete entries. Warn about multiple nearby entries in same gene.
    :param data: list of dictionaries with STR data
    :param genome: genome build (hg19, hg38 or t2t)
    """
    keep_rows = []
    for row in data:
        complete = True
        for field in ['chrom', 'start_' + genome, 'stop_' + genome, 'pathogenic_motif_reference_orientation', 'gene', 'id']:
            if field not in row:
                raise ValueError(f'Missing field {field} in input file.')
            if row[field].strip() == '':
                sys.stderr.write(f'Missing value for field {field} in input file. Skipping {row["id"]}\n')
                complete = False
        if complete:
            keep_rows.append(row)

    # report multiple entries in same gene
    genes = [row['gene'] for row in keep_rows]
    duplicates = ', '.join(set([x for x in genes if genes.count(x) > 1]))
    sys.stderr.write(f'Warning: multiple loci found in the same gene, keeping all: {duplicates}\n')

    return keep_rows

def trgt_catalog(row, genome = 'hg38'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :return: TRGT format catalog string

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': 'CAG', 'gene': 'mygene', 'id': 'myid', 'locus_structure': '', 'flank_motif': ''})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG;STRUC=(CAG)n'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': 'CAG,CCG', 'gene': 'mygene', 'id': 'myid', 'locus_structure': '', 'flank_motif': ''})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG,CCG;STRUC=(CAG)n(CCG)n'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': 'CAGG', 'gene': 'CNBP', 'id': 'DM2_CNBP', 'locus_structure': '(CAGG)*(CAGA)*(CA)*', 'flank_motif': '(CAGG)n(CAGA)10(CA)19'})
    'chr1\t100\t278\tID=DM2_CNBP;MOTIFS=CAGG,CAGA,CA;STRUC=(CAGG)n(CAGA)n(CA)n'
    """
    start = int(row['start_' + genome])
    stop = int(row['stop_' + genome])
    struc = ''

    if row['flank_motif'] != '':
        # get motifs in parentheses using regex
        flank_motif = row['flank_motif']
        motifs = re.findall('\((.*?)\)', flank_motif)
        counts = re.findall('\)(.*?)[\(|$]', flank_motif.replace('n', 'n(') + '$')
        n_found = False
        for motif, count in zip(motifs, counts):
            struc += f'({motif})n'
            if count == 'n':
                n_found = True
                continue
            else:
                if n_found:
                    stop += int(count) * len(motif)
                else:
                    start -= int(count) * len(motif)
    elif row['locus_structure'] != '':
        locus_structure = row['locus_structure'].strip()
        motifs = re.findall('\((.*?)\)', locus_structure)
        # Substitute * and + with n
        struc = locus_structure.replace('*', 'n').replace('+', 'n')
    else:
        motifs = []
        for motif in row['pathogenic_motif_reference_orientation'].split(','):
            motif = motif.strip() # remove leading and trailing whitespace
            struc += f'({motif})n'
            motifs.append(motif)
    if row['gene'] == 'RFC1':
        struc = '<RFC1>' # special case for RFC1 coded as HMM by TRGT. May be removed in future versions of TRGT.
    # unique motifs mainitaining order
    motifs = list(dict.fromkeys(motifs))
    definition = f"{row['chrom']}\t{start}\t{stop}\tID={row['id']};MOTIFS={','.join(motifs)};STRUC={struc}"

    return definition

def main(input: str, output: str, *, format: str = 'TRGT', genome: str = 'hg38'):
    """
    :param input: STRchive database file name (CSV or JSON)
    :param output: Output file name
    :param genome: Genome build (hg19, hg38 or T2T)
    :param format: Variant caller catalog file format (TRGT)
    """

    genome = genome.lower()

    if input.lower().endswith('.csv'):
        with open(input, 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            # convert to dictionary
            data = []
            for row in csv_reader:
                row_dict = {}
                for k,v in row.items():
                    row_dict[k] = v
                data.append(row_dict)
    elif input.lower().endswith('.json'):
        with open(input, 'r') as json_file:
            data = json.load(json_file)
            
    else:
        raise ValueError('Unknown input file format. Expected .csv or .json.')

    data = clean_loci(data, genome)

    # sort by chromosome and start position
    data = sorted(data, key = lambda x: (chrom_to_int(x['chrom']), int(x['start_' + genome])))

    if format.upper() == 'TRGT':
        with open(output, 'w') as out_file:
            for row in data:
                out_file.write(trgt_catalog(row, genome) + '\n')
    else:
        raise ValueError('Unknown output file format. Expected TRGT.')
    


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)