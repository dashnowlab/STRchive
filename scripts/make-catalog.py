import json
import sys
import re
import decimal

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

def remove_exponent(d):
    """
    >>> remove_exponent(5.00)
    Decimal('5')

    >>> remove_exponent(5.500)
    Decimal('5.5')

    >>> remove_exponent(10000.0)
    Decimal('10000')
    """
    d = decimal.Decimal(d)
    return d.quantize(decimal.Decimal(1)) if d == d.to_integral() else d.normalize()

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
            if row[field] is None or row[field] == '':
                sys.stderr.write(f'Missing value for field {field} in input file. Skipping {row["id"]}\n')
                complete = False
        if complete:
            keep_rows.append(row)

    # report multiple entries in same gene
    genes = [row['gene'] for row in keep_rows]
    duplicates = ', '.join(set([x for x in genes if genes.count(x) > 1]))
    sys.stderr.write(f'Warning: multiple loci found in the same gene, keeping all: {duplicates}\n')

    return keep_rows

def trgt_catalog(row, genome = 'hg38', struc_type = 'default'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :param struc_type: options: 'motif', 'default' or 'none'. If 'motif', use pathogenic_motif_reference_orientation as locus structure. If 'default', use <TR>. If 'none', do not include locus structure.
    :return: TRGT format catalog string

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': ['CAG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': '', 'flank_motif': '', 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG;STRUC=<TR>'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': ['CAG', 'CCG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': '', 'flank_motif': '', 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG,CCG;STRUC=<TR>'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'period': '3', 'pathogenic_motif_reference_orientation': ['CAGG'], 'gene': 'CNBP', 'id': 'DM2_CNBP', 'locus_structure': '(CAGG)*(CAGA)*(CA)*', 'flank_motif': '(CAGG)n(CAGA)10(CA)19', 'reference_motif_reference_orientation': [], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t278\tID=DM2_CNBP;MOTIFS=CAGG,CAGA,CA;STRUC=<TR>'
    """
    start = int(row['start_' + genome])
    stop = int(row['stop_' + genome])
    struc = ''

    if row['flank_motif'] != '' and row['flank_motif'] is not None:
        # get motifs in parentheses using regex
        flank_motif = row['flank_motif']
        motifs = re.findall(r'\((.*?)\)', flank_motif)
        counts = re.findall(r'\)(.*?)[\(|$]', flank_motif.replace('n', 'n(') + '$')
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
    elif row['locus_structure'] != '' and row['locus_structure'] != None:
        locus_structure = row['locus_structure'].strip()
        motifs = re.findall(r'\((.*?)\)', locus_structure)
        # Substitute * and + with n
        struc = locus_structure.replace('*', 'n').replace('+', 'n')
    else:
        motifs = []
        for motif in row['pathogenic_motif_reference_orientation']:
            motif = motif.strip() # remove leading and trailing whitespace
            struc += f'({motif})n'
            motifs.append(motif)
    if row['gene'] == 'RFC1':
        struc = '<RFC1>' # special case for RFC1 coded as HMM by TRGT. May be removed in future versions of TRGT.
    # unique motifs mainitaining order

    if struc_type == 'motif':
        full_struc = f";STRUC{struc}"
    elif struc_type == 'default':
        # Add all motifs with known function
        for motif_field in ['pathogenic_motif_reference_orientation', 'reference_motif_reference_orientation', 'benign_motif_reference_orientation']:
            motifs.extend(row[motif_field])
        motifs = list(dict.fromkeys(motifs)) # remove duplicates
        full_struc = ';STRUC=<TR>'
    elif struc_type == 'none':
        full_struc = ''
    else:
        raise ValueError(f'Unknown structure type: {struc_type}. Expected motif, default or none.')

    motifs = list(dict.fromkeys(motifs))
    definition = f"{row['chrom']}\t{start}\t{stop}\tID={row['id']};MOTIFS={','.join(motifs)}{full_struc}"

    return definition

def atarva_catalog(row, genome = 'hg38'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :return: atarva format catalog string which is a modified BED format with fields: chrom start stop motif motif_len [id]

    Note, compound loci will be split into multiple entries, one for each motif. Overlapping loci are okay.
    For loci with multiple pathogenic motifs, only the first motif will be used. Atarva does motif decomposition, so alternate motifs should be detected by the caller.

    >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'flank_motif': '', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    'chr1\t100\t200\tCAG\t3\tmyid'

    >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['AAGGG', 'ACAGG'], 'flank_motif': '', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    'chr1\t100\t200\tAAGGG\t5\tmyid'

    >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'flank_motif': '(CAG)nCAACAG(CCG)12', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    'chr1\t100\t200\tCAG\t3\tmyid\nchr1\t200\t206\tCAACAG\t6\tmyid_flank\nchr1\t206\t242\tCCG\t3\tmyid_flank'

    >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'flank_motif': '(CAG)n(CCG)10(CAA)10', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    'chr1\t100\t200\tCAG\t3\tmyid\nchr1\t200\t230\tCCG\t3\tmyid_flank\nchr1\t230\t260\tCAA\t3\tmyid_flank'
    """
    bed_string = ''

    motif_field = 'pathogenic_motif_reference_orientation'
    id_field = 'id'
    start = int(row['start_' + genome])
    stop = int(row['stop_' + genome])

    motifs = row[motif_field]
    this_id = row[id_field]

    motif = motifs[0] # use first motif only
    motif_len = len(motif)
    bed_string += f"{row['chrom']}\t{start}\t{stop}\t{motif}\t{motif_len}\t{this_id}\n"

    # check for flanking motif(s)
    if row['flank_motif'] != '' and row['flank_motif'] is not None:
    # get motifs in parentheses using regex
        flank_motif = row['flank_motif']
    
        split_flank_counts = re.split(r"[ATCGN]+", flank_motif, )
        all_flank_motifs_counts = []
        for i in range(1, len(split_flank_counts)):
            all_flank_motifs_counts.append(split_flank_counts[i].replace('(', '').replace(')', ''))

        all_flank_motifs = ''
        for char in flank_motif:
            if char in ['(', ')', 'n', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                all_flank_motifs += ' '
            else:
                all_flank_motifs += char

        all_flank_motifs = all_flank_motifs.split()

        flank_start = stop
        flank_stop = stop
        for motif, count in zip(all_flank_motifs, all_flank_motifs_counts):
            if count == '':
                count = 1
            if count == 'n':
                continue
            else:
                flank_stop += int(count) * len(motif)
            bed_string += f"{row['chrom']}\t{flank_start}\t{flank_stop}\t{motif}\t{len(motif)}\t{this_id}_flank\n"
            flank_start = flank_stop

    return bed_string.rstrip('\n')

def longtr_catalog(row, genome = 'hg38'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :return: LongTR format catalog string

    Note, LongTR uses 1-based coordinates (i.e. is non-standard BED format)
    """
    start = int(row['start_' + genome]) + 1
    stop = int(row['stop_' + genome])
    motifs = row['pathogenic_motif_reference_orientation'] + row['benign_motif_reference_orientation'] + row['reference_motif_reference_orientation']
    # remove duplicates
    motifs = list(dict.fromkeys([x.strip() for x in motifs]))
    motifs_string = ','.join(motifs)

    definition = f"{row['chrom']}\t{start}\t{stop}\t{motifs_string}\t{row['id']}"

    return definition


def extended_bed(row, fields = [], genome = 'hg38'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param fields: list of fields to include in the extended BED format beyond chrom, start, stop
    :param genome: genome build (hg19, hg38 or T2T)
    :return: BED format catalog string

    >>> extended_bed({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'pathogenic_motif_reference_orientation': 'CAG', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': '10', 'inheritance': 'AD', 'disease': 'Disease Name'}, ['id', 'gene', 'pathogenic_motif_reference_orientation', 'pathogenic_min', 'inheritance', 'disease'])
    'chr1\t100\t200\tmyid\tmygene\tCAG\t10\tAD\tDisease Name'

    >>> extended_bed({'chrom': 'chr1', 'start_hg38': '100', 'stop_hg38': '200', 'pathogenic_motif_reference_orientation': 'CAG', 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': '10', 'inheritance': 'AD', 'disease': 'Disease Name'})
    'chr1\t100\t200'
    """
    start = int(row['start_' + genome])
    stop = int(row['stop_' + genome])

    bed_string = f"{row['chrom']}\t{start}\t{stop}"
    if len(fields) > 0:
        for field in fields:
            if isinstance(row[field], list):
                row[field] = ','.join(row[field])
            if isinstance(row[field], float):
                row[field] = remove_exponent(row[field])
            # ensure field does not contain tabs
            if isinstance(row[field], str) and '\t' in row[field]:
                raise ValueError(f'Tab character found in field {field} value: {row[field]}')
            bed_string += f"\t{row[field]}" 
    return bed_string

default_fields = ','.join(['id', 'gene', 'reference_motif_reference_orientation', 'pathogenic_motif_reference_orientation', 'pathogenic_min', 'inheritance', 'disease'])

def main(input: str, output: str, *, format: str = 'TRGT', genome: str = 'hg38', cols: str = default_fields):
    """
    :param input: STRchive database file name in JSON format
    :param output: Output file name in bed format
    :param genome: Genome build: hg19, hg38, T2T (also accepted: chm13, chm13-T2T, T2T-CHM13)
    :param format: Variant caller catalog file format BED format (TRGT, LongTR or BED)
    :param cols: Comma separated list of columns to include in the extended BED format beyond chrom,start,stop (no spaces in list). Can be any valid STRchive json field.
    """

    genome = genome.lower()
    if genome not in ['hg19', 'hg38', 't2t', 'chm13', 'chm13-t2t', 't2t-chm13']:
        raise ValueError(f'Unknown genome build: {genome}\nExpected hg19, hg38, T2T, chm13, chm13-T2T or T2T-chm13.')
    if genome in ['chm13', 'chm13-t2t', 't2t-chm13']:
        genome = 't2t'

    if input.lower().endswith('.json'):
        with open(input, 'r') as json_file:
            data = json.load(json_file)    
    else:
        raise ValueError(f'Unknown input file extension: {input} \nExpected .json')

    fields = cols
    if fields != default_fields and format.lower() != 'bed':
        raise ValueError('Fields option is only available for BED format output.')

    data = clean_loci(data, genome)

    # sort by chromosome and start position
    data = sorted(data, key = lambda x: (chrom_to_int(x['chrom']), int(x['start_' + genome])))

    if format.lower() == 'trgt':
        with open(output, 'w') as out_file:
            for row in data:
                out_file.write(trgt_catalog(row, genome) + '\n')
    elif format.lower() == 'atarva':
        with open(output, 'w') as out_file:
            header = '#' + '\t'.join(['chrom', 'start', 'stop', 'motif', 'motif_len', 'id']) + '\n'
            out_file.write(header)
            for row in data:
                out_file.write(atarva_catalog(row, genome) + '\n')
    elif format.lower() == 'longtr':
        with open(output, 'w') as out_file:
            for row in data:
                out_file.write(longtr_catalog(row, genome) + '\n')
    elif format.lower() == 'bed':
        fields_list = fields.split(',')
        header = '#' + '\t'.join(['chrom', 'start', 'stop'] + fields_list) + '\n'
        with open(output, 'w') as out_file:
            out_file.write(header)
            for row in data:
                out_file.write(extended_bed(row, fields_list, genome) + '\n')
    else:
        raise ValueError('Unknown output file format. Expected TRGT or BED.')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)
