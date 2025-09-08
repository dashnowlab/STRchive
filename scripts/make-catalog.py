import json
import sys
import re
import decimal
import jsbeautifier
from copy import deepcopy

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

def add_flank_coordinates(row, genome = 'hg38'):
    """
    Get the start and stop coordinates of the flanking motif(s) in the locus structure. Only one of the motifs can be missing a count, the rest must have a count.

    example input

    "locus_structure": [
        {
        "motif": "CAG",
        "count": null,
        "type": "pathogenic_repeat"
        },
        {
        "motif": "CAACAG",
        "count": 1,
        "type": "interruption"
        },
        {
        "motif": "CCG",
        "count": 12,
        "type": "flank_repeat"
        }
    ]

    example output
        "locus_structure": [
        {
        "motif": "CAG",
        "count": null,
        "type": "pathogenic_repeat",
        "position": "chr1:100-200",
        "location": "internal",
        "length": 100
        },
        {
        "motif": "CAACAG",
        "count": 1,
        "type": "interruption",
        "position": "chr1:200-206",
        "location": "internal",
        "length": 6
        },
        {
        "motif": "CCG",
        "count": 12,
        "type": "flank_repeat",
        "position": "chr1:206-242",
        "location": "right",
        "length": 36
        }
    ]

    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :return: updated locus_structure with positions for each motif


    >>> add_flank_coordinates({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'locus_structure': []}, genome='hg38')
    []

    >>> add_flank_coordinates({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAACAG', 'count': 1, 'type': 'interruption'}, {'motif': 'CCG', 'count': 12, 'type': 'flank_repeat'}]}, genome='hg38')
    [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat', 'length': 94, 'location': 'internal', 'start_hg38': 100, 'stop_hg38': 194}, {'motif': 'CAACAG', 'count': 1, 'type': 'interruption', 'length': 6, 'location': 'internal', 'start_hg38': 194, 'stop_hg38': 200}, {'motif': 'CCG', 'count': 12, 'type': 'flank_repeat', 'length': 36, 'location': 'right', 'start_hg38': 200, 'stop_hg38': 236}]
    
    """

    if len(row['locus_structure']) == 0:
        # No locus structure, return empty list
        return []

    new_locus_structure = deepcopy(row['locus_structure'])

    # Locus structure includes:
    # - pathogenic_repeat: unknown size, inside the STRchive coordinates - the ref size can be inferred from the coordinates of everything else
    # - internal_repeat: defined count, inside the STRchive coordinates
    # - interruption: defined count, inside the STRchive coordinates
    # - flank_repeat: defined count, OUTSIDE the STRchive coordinates - can be left or right flanks (determined from the order)
    # The logic:
    # - Calculate the total length of the STR region by adding the left and right flank lengths if applicable
    # - For each motif in the locus structure:
    #   - If it's a pathogenic_repeat, we can't determine its length yet
    #   - If it's an internal_repeat or interruption, we can use its count to determine its length
    #   - If it's a flank_repeat, we need to determine if it's on the left or right flank

    internal_length = row['stop_' + genome] - row['start_' + genome]
    pathogenic_repeat_found = False
    for struct_dict in new_locus_structure:
        if struct_dict['type'] == 'pathogenic_repeat' and struct_dict['count'] is None: # Do I need to assume it is none?
            pathogenic_repeat_found = True
            struct_dict['length'] = struct_dict['count']
            struct_dict['location'] = 'internal'
        elif struct_dict['type'] == 'internal_repeat' or struct_dict['type'] == 'interruption':
            struct_dict['length'] = len(struct_dict['motif']) * struct_dict['count']
            struct_dict['location'] = 'internal'
        elif struct_dict['type'] == 'flank_repeat':
            struct_dict['length'] = len(struct_dict['motif']) * struct_dict['count']
            if pathogenic_repeat_found:
                # Add to right flank
                struct_dict['location'] = 'right'
            else:
                # Add to left flank
                struct_dict['location'] = 'left'
        else:
            # Error case: if a motif is not a main repeat and has no count, it should not be in the locus structure
            raise ValueError(f"Motif {struct_dict['motif']} type is not defined. Please check the input data.")

    # Check if more than one unknown length is present
    known_lengths = [struct_dict['length'] for struct_dict in new_locus_structure if struct_dict['location'] == 'internal']
    left_flank = sum([struct_dict['length'] for struct_dict in new_locus_structure if struct_dict['location'] == 'left'])
    right_flank = sum([struct_dict['length'] for struct_dict in new_locus_structure if struct_dict['location'] == 'right'])
    if known_lengths.count(None) > 1:
        raise ValueError(f"Multiple unknown lengths found in locus structure for {row['id']}. Please check the input data.")
    unknown_length = internal_length - sum([x for x in known_lengths if x is not None])
    #print(f"Total length of locus structure for {row['id']} is {total_length}, known lengths are {known_lengths}, flank lengths are {flank_length}, unknown length is {unknown_length}.")
    if unknown_length < 0:
        # sys.stdout.write(f'{row}\n\n')
        # sys.stdout.write(f'{new_locus_structure}\n\n')
        raise ValueError(f"Total length of locus structure for {row['id']} is less than the sum of known lengths. Please check the input data.")

    # Iterate over the new locus structure again
    # Calculate start and stop coordinates for each motif in the locus structure
    this_start = row['start_' + genome] - left_flank
    for struct_dict in new_locus_structure:
        if struct_dict['length'] is None:
            struct_dict['length'] = unknown_length
        this_stop = this_start + struct_dict['length']
        struct_dict['start_' + genome] = this_start
        struct_dict['stop_' + genome] = this_stop
        this_start = this_stop
    
    # Check that the end coordinate of the last item matches the right flank
    if this_start != row['stop_' + genome] + right_flank:
        raise ValueError(f"End coordinate of last item in locus structure for {row['id']} does not match right flank. Probably a logic error.")

    return new_locus_structure

def trgt_catalog(row, genome = 'hg38', struc_type = 'default'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :param struc_type: options: 'motif', 'default' or 'none'. If 'motif', use pathogenic_motif_reference_orientation as locus structure. If 'default', use <TR>. If 'none', do not include locus structure.
    :return: TRGT format catalog string

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'period': 3, 'pathogenic_motif_reference_orientation': ['CAG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [], 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG;STRUC=<TR>'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'period': 3, 'pathogenic_motif_reference_orientation': ['CAG', 'CCG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [], 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG,CCG;STRUC=<TR>'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'period': 3, 'pathogenic_motif_reference_orientation': ['CAG', 'CCG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAACAG', 'count': 1, 'type': 'interruption'}, {'motif': 'CCG', 'count': 12, 'type': 'flank_repeat'}], 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG,CCG;STRUC=<TR>'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'period': 3, 'pathogenic_motif_reference_orientation': ['CAGG'], 'gene': 'CNBP', 'id': 'DM2_CNBP', 'locus_structure': [{'motif': 'CAGG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAGA', 'count': 10, 'type': 'flank_repeat'}, {'motif': 'CA', 'count': 19, 'type': 'flank_repeat'}], 'reference_motif_reference_orientation': [], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []}, struc_type='motif')
    'chr1\t100\t200\tID=DM2_CNBP;MOTIFS=CAGG,CAGA,CA;STRUC=(CAGG)n(CAGA)10(CA)19'

    >>> trgt_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'period': 3, 'pathogenic_motif_reference_orientation': ['CAG', 'CCG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [], 'reference_motif_reference_orientation': ['CAG'], 'benign_motif_reference_orientation': [], 'unknown_motif_reference_orientation': []})
    'chr1\t100\t200\tID=myid;MOTIFS=CAG,CCG;STRUC=<TR>'
    """

    # Add flank coordinates to locus structure
    row['locus_structure'] = add_flank_coordinates(row, genome)

    start = row['start_' + genome]
    stop = row['stop_' + genome]
    struc = ''
    motifs = []

    if len(row['locus_structure']) > 0:
        for struct_dict in row['locus_structure']:
            if struct_dict['type'] == 'interruption':
                # interruptions are not included in the structure
                continue
            motifs.append(struct_dict['motif'])
            if struct_dict['count'] is None:
                struc += f"({struct_dict['motif']})n"
            else:
                struc += f"({struct_dict['motif']}){struct_dict['count']}"

    else: # should just do this always?
        for motif in row['pathogenic_motif_reference_orientation']:
            motif = motif.strip() # remove leading and trailing whitespace
            struc += f'({motif})n'
            motifs.append(motif)
    if row['gene'] == 'RFC1':
        struc = '<RFC1>' # special case for RFC1 coded as HMM by TRGT. May be removed in future versions of TRGT.
    # unique motifs mainitaining order

    if struc_type == 'motif':
        full_struc = f";STRUC={struc}"
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

    # >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tCAG\t3\tmyid'

    # >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['AAGGG', 'ACAGG'], 'locus_structure': [], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tAAGGG\t5\tmyid'

    # >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAACAG', 'count': 2, 'type': 'flank_repeat'}, {'motif': 'CCG', 'count': 3, 'type': 'flank_repeat'}], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tCAG\t3\tmyid\nchr1\t200\t212\tCAACAG\t6\tmyid_flank\nchr1\t212\t221\tCCG\t3\tmyid_flank'

    # >>> atarva_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAACAG', 'count': 6, 'type': 'flank_repeat'}, {'motif': 'CCG', 'count': 3, 'type': 'flank_repeat'}], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tCAG\t3\tmyid\nchr1\t200\t236\tCAACAG\t6\tmyid_flank\nchr1\t236\t245\tCCG\t3\tmyid_flank'
    """

    # Add flank coordinates to locus structure
    row['locus_structure'] = add_flank_coordinates(row, genome)

    bed_string = ''

    motif_field = 'pathogenic_motif_reference_orientation'
    id_field = 'id'
    start = int(row['start_' + genome])
    stop = int(row['stop_' + genome])

    motifs = row[motif_field]
    this_id = row[id_field]

    # check for flanking motif(s)
    if len(row['locus_structure']) > 0:
        for struct_dict in row['locus_structure']:
            motif = struct_dict['motif']
            motif_len = len(motif)
            start = struct_dict['start_' + genome]
            stop = struct_dict['stop_' + genome]
            if struct_dict['type'] == 'pathogenic_repeat':
                # this is the main repeat
                bed_string += f"{row['chrom']}\t{start}\t{stop}\t{motif}\t{motif_len}\t{this_id}\n"
            elif struct_dict['type'] == 'interruption':
                # interruptions are not included in the structure
                continue
            else:
                # this is a flank repeat
                bed_string += f"{row['chrom']}\t{start}\t{stop}\t{motif}\t{motif_len}\t{this_id}_flank\n"

    else:
        motif = motifs[0] # use first motif only
        motif_len = len(motif)
        bed_string += f"{row['chrom']}\t{start}\t{stop}\t{motif}\t{motif_len}\t{this_id}\n"

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

def expansionhunter_catalog(row, genome):
    """
    See format description at https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md

    Example from gnomAD:
    {
        "LocusId": "ABCD3",
        "ReferenceRegion": "chr1:94418421-94418442",
        "LocusStructure": "(GCC)*",
        "VariantType": "Repeat",
        "RepeatUnit": "GCC",
        "Gene": "ABCD3",
        "GeneRegion": "5'-UTR",
        "GeneId": "ENSG00000117528",
        "DiscoveryMethod": "WGS",
        "DiscoveryYear": 2023,
        "Diseases": [
            {
                "Symbol": "OPDM",
                "Name": "Oculopharyngodistal myopathy",
                "Inheritance": "AD",
                "NormalMax": 44,
                "PathogenicMin": 118
            }
        ],
        "MainReferenceRegion": "chr1:94418421-94418442",
        "Inheritance": "AD"
    }

    Simple example from ExpansionHunter:
    {
    "LocusId": "DMPK",
    "LocusStructure": "(CAG)*",
    "ReferenceRegion": "19:46273462-46273522", # 0-based coordinates
    "VariantType": "Repeat"
    },
    {
    "LocusId": "HTT",
    "LocusStructure": "(CAG)*CAACAG(CCG)*",
    "ReferenceRegion": ["4:3076604-3076660", "4:3076666-3076693"],
    "VariantType": ["Repeat", "Repeat"]
    }
    """
    raise NotImplementedError

    # Optional fields used by gnomAD:
    # locus_dict['MainReferenceRegion'] = f"{row['chrom']}:{row['start_' + genome]}-{row['stop_' + genome]}"
    # locus_dict['Inheritance'] = row['inheritance']
    # locus_dict['Gene'] = row['gene']
    # locus_dict['GeneRegion'] = row['type']
    # locus_dict['DiscoveryYear'] = row['year']
    # locus_dict['Diseases'] = [{
    #     'Symbol': row['disease_id'],
    #     'Name': row['disease'],
    #     'Inheritance': row['inheritance'],
    #     'NormalMax': row['benign_max'],
    #     'PathogenicMin': row['pathogenic_min']
    # }]

def stranger_catalog(row, genome = 'hg38'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :return: STRanger format catalog dictionary for a single locus

    Note, the stranger catalog is similar to the ExpansionHunter catalog and in some cases they are used for both purposes.
    However, the STRanger catalog is for annotation, and in this case is designed to be used with long-read genotype data.
    It is critical that the start coordinates of the stranger catalog match the start coordinates of the catalog used for genotyping.
    For this reason, no flanking coordinates are added to the locus structure.

    Example from STRanger:
    {
        "LocusId": "ABCD3",
        "HGNCId": 67,                     # Currently not implemented
        "InheritanceMode": "AD",
        "DisplayRU": "CCG",
        "LocusStructure": "(CCG)*",
        "ReferenceRegion": "1:94418422-94418444",
        "VariantType": "Repeat",
        "Disease": "OPDM",
        "NormalMax": 50,
        "PathologicMin": 118
    },
    {
        "LocusId": "FXN",
        "HGNCId": 3951,
        "InheritanceMode": "AR",
        "DisplayRU": "GAA",
        "SourceDisplay": "GeneReviews Internet 2019-11-07",
        "Source": "GeneReviews",
        "SourceId": "NBK535148",
        "LocusStructure": "(A)*(GAA)*",
        "ReferenceRegion": "chr9:69037286-69037304",
        "VariantId": "FXN",
        "VariantType": "Repeat",
        "Disease": "FRDA",
        "NormalMax": 35,
        "PathologicMin": 51,
        "PathologicRegion": "chr9:69037286-69037304"
    },

    # >>> stranger_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [], 'benign_max': 5, 'pathogenic_min': 10, 'inheritance': 'AD', 'disease_id': 'disease_id', 'disease': 'Disease Name', 'year': 2023}, genome='hg38')
    # {'LocusId': 'myid', 'ReferenceRegion': 'chr1:100-200', 'LocusStructure': '(CAG)*', 'VariantType': 'Repeat', 'HGNCId': None, 'InheritanceMode': 'AD', 'DisplayRU': 'CAG', 'Disease': 'disease_id', 'NormalMax': 5, 'PathologicMin': 10, 'Gene': 'mygene'}

    # >>> stranger_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'gene': 'mygene', 'id': 'myid', 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CAACAG', 'count': 1, 'type': 'interruption'}, {'motif': 'CCG', 'count': 3, 'type': 'flank_repeat'}], 'benign_max': 5, 'pathogenic_min': 10, 'inheritance': 'AD', 'disease_id': 'disease_id', 'disease': 'Disease Name', 'year': 2023}, genome='hg38')
    # {'LocusId': 'myid', 'ReferenceRegion': ['chr1:100-200', 'chr1:206-215'], 'LocusStructure': '(CAG)*CAACAG(CCG)*', 'VariantType': ['Repeat', 'Repeat'], 'VariantId': ['myid', 'myid_CCG'], 'PathologicRegion': 'chr1:100-200', 'HGNCId': None, 'InheritanceMode': 'AD', 'DisplayRU': 'CAG', 'Disease': 'disease_id', 'NormalMax': 5, 'PathologicMin': 10, 'Gene': 'mygene'}
    """

    row['locus_structure'] = add_flank_coordinates(row, genome)

    locus_dict = {}

    # Required/standard fields from ExpansionHunter:
    locus_dict['LocusId'] = row['id']
    locus_dict['ReferenceRegion'] = []
    if len(row['locus_structure']) > 0:
        locus_dict['LocusStructure'] = ''
        locus_dict['VariantType'] = []
        locus_dict['VariantId'] = [] # used to store the ID of the variant, e.g. myid_CCG for the CCG motif in the locus structure
        for struct_dict in row['locus_structure']:
            
            if struct_dict['type'] == 'interruption':
                locus_dict['LocusStructure'] += f"{struct_dict['motif']*struct_dict['count']}" # interruptions are included in the structure but not in the variant list
            else:
                locus_dict['ReferenceRegion'].append(f"{row['chrom']}:{struct_dict['start_' + genome]}-{struct_dict['stop_' + genome]}") # 0-based coordinates
                locus_dict['LocusStructure'] += f"({struct_dict['motif']})*"
                locus_dict['VariantType'].append('Repeat')
                if struct_dict['type'] == 'pathogenic_repeat':
                    locus_dict['VariantId'].append(row['id'])
                    locus_dict['PathologicRegion'] = f"{row['chrom']}:{struct_dict['start_' + genome]}-{struct_dict['stop_' + genome]}"
                else:
                    locus_dict['VariantId'].append(f"{row['id']}_{struct_dict['motif']}")
    else:
        locus_dict['ReferenceRegion'] = f"{row['chrom']}:{row['start_' + genome]}-{row['stop_' + genome]}" # 0-based coordinates
        locus_dict['LocusStructure'] = f"({row['pathogenic_motif_reference_orientation'][0]})*" # Just use the first pathogenic motif for the structure
        locus_dict['VariantType'] = 'Repeat'

    # Fields used by STRanger:
    locus_dict['HGNCId'] = None # Currently not implemented, would require mapping of gene names to HGNC IDs
    locus_dict['InheritanceMode'] = row['inheritance']
    locus_dict['DisplayRU'] = row['pathogenic_motif_reference_orientation'][0] # use first pathogenic motif for display
    locus_dict['Disease'] = row['disease_id'] # disease ID
    

    # Do some special handling for missing values and unusual cases:
    # If both pathogenic min and benign max are missing, skip the locus
    if row['pathogenic_min'] is None and row['benign_max'] is None:
        return None

    # If normal max is missing, set it to pathogenic min - 1
    if row['benign_max'] is None:
        locus_dict['NormalMax'] = row['pathogenic_min'] - 1
    else:
        locus_dict['NormalMax'] = row['benign_max']

    # If pathogenic min is <= benign max or missing, set pathologic min to benign max + 1
    if row['pathogenic_min'] is None or row['pathogenic_min'] <= locus_dict['NormalMax']:
        locus_dict['PathologicMin'] = locus_dict['NormalMax'] + 1
    else:
        locus_dict['PathologicMin'] = row['pathogenic_min']

    # Optional extra fields:
    locus_dict['Gene'] = row['gene']

    return locus_dict

def straglr_catalog(row, genome = 'hg38', format = 'default'):
    r"""
    :param row: dictionary with STR data for a single locus
    :param genome: genome build (hg19, hg38 or T2T)
    :param format: options: 'default' or 'wf-human-variation'. If 'default', use chrom, start, stop, motif. If 'wf-human-variation', use chrom, start, stop, motif, gene, id.
    :return: straglr format catalog string which is a modified BED format with fields: chrom start stop motif [gene id]

    Example format:
    https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed

    chr1	149390802	149390841	GGC	NOTCH2NLC	NOTCH2NLC
    chr2	190880872	190880920	GCA	GLS	GLS
    chr3	63912684	63912714	GCA	ATXN7	ATXN7

    # >>> straglr_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CCG', 'count': 10, 'type': 'flank_repeat'}, {'motif': 'CAA', 'count': 10, 'type': 'flank_repeat'}], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tCAG\nchr1\t200\t230\tCCG\nchr1\t230\t260\tCAA'

    # >>> straglr_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CCG', 'count': 10, 'type': 'flank_repeat'}, {'motif': 'CAA', 'count': 10, 'type': 'flank_repeat'}], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38')
    # 'chr1\t100\t200\tCAG\nchr1\t200\t230\tCCG\nchr1\t230\t260\tCAA'

    # >>> straglr_catalog({'chrom': 'chr1', 'start_hg38': 100, 'stop_hg38': 200, 'pathogenic_motif_reference_orientation': ['CAG'], 'locus_structure': [{'motif': 'CAG', 'count': None, 'type': 'pathogenic_repeat'}, {'motif': 'CCG', 'count': 10, 'type': 'flank_repeat'}, {'motif': 'CAA', 'count': 10, 'type': 'flank_repeat'}], 'gene': 'mygene', 'id': 'myid', 'pathogenic_min': 10, 'inheritance': 'AD', 'disease': 'Disease Name'}, 'hg38', format='wf-human-variation')
    # 'chr1\t100\t200\tCAG\tmyid\tmyid\nchr1\t200\t230\tCCG\tmyid\tmyid_CCG\nchr1\t230\t260\tCAA\tmyid\tmyid_CAA'
    """

    # Do some special handling for missing values and unusual cases:
    # If both pathogenic min and benign max are missing, skip the locus
    if row['pathogenic_min'] is None and row['benign_max'] is None:
        return None
    
    bed_list = []
    # Use the same approach as atarva_catalog, but only return the first 4 columns (chrom, start, stop)
    atarva_list = [x.split('\t') for x in atarva_catalog(row, genome).split('\n')]  
    if format == 'default':
        for bed_row in atarva_list:
            bed_list.append('\t'.join(bed_row[0:4]))
    elif format == 'wf-human-variation':
        for bed_row in atarva_list:
            motif = bed_row[3]
            # replace the word "flank" with the motif name in the id
            bed_row[5] = bed_row[5].replace('_flank', f'_{motif}')
            bed_list.append('\t'.join(bed_row[0:4] + [row['id'], bed_row[5]]))
    else:
        raise ValueError(f'Unknown format: {format}. Expected default or wf-human-variation.')

    bed_string = '\n'.join(bed_list)

    return bed_string

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
    :param format: Variant caller catalog file format or BED format (TRGT, atarva, LongTR, straglr, stranger, ExpansionHunter or BED)
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
    elif format.lower() == 'expansionhunter':
        eh_loci = []
        for row in data:
            locus = expansionhunter_catalog(row, genome)
            eh_loci.append(locus)
        # Write the catalog as a JSON array
        output = output if output.endswith('.json') else output + '.json'
        with open(output, 'w') as out_json_file:
            options = jsbeautifier.default_options()
            options.indent_size = 2
            options.brace_style="expand"
            out_json_file.write(jsbeautifier.beautify(json.dumps(eh_loci, ensure_ascii=False), options))
            out_json_file.write('\n')
    elif format.lower() == 'stranger':
        stranger_loci = []
        for row in data:
            locus = stranger_catalog(row, genome)
            if locus is not None:
                stranger_loci.append(locus)
        # Write the catalog as a JSON array
        output = output if output.endswith('.json') else output + '.json'
        with open(output, 'w') as out_json_file:
            options = jsbeautifier.default_options()
            options.indent_size = 2
            options.brace_style="expand"
            out_json_file.write(jsbeautifier.beautify(json.dumps(stranger_loci, ensure_ascii=False), options))
            out_json_file.write('\n')
    elif format.lower() == 'straglr':
        with open(output, 'w') as out_file:
            # No header for straglr format
            for row in data:
                straglr_string = straglr_catalog(row, genome, format = 'wf-human-variation')
                if straglr_string is not None:  # Check if the string is not None
                    out_file.write(straglr_string + '\n')
    elif format.lower() == 'bed':
        fields_list = fields.split(',')
        header = '#' + '\t'.join(['chrom', 'start', 'stop'] + fields_list) + '\n'
        with open(output, 'w') as out_file:
            out_file.write(header)
            for row in data:
                out_file.write(extended_bed(row, fields_list, genome) + '\n')
    else:
        raise ValueError('Unknown output file format. Expected TRGT, atarva, straglr or BED.')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)
