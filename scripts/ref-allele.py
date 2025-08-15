import pysam
import re
import sys
import os
import urllib.request

def split_repeat_sequence(motifs, sequence):
    """
    Splits repeat sequence based on the defined motifs

    :param motifs: the list of motifs defined for the locus
    :param sequence: sequence of the locus
    :return: the formatted sequence string split based on the motifs
    """
    
    # initialisers
    position = 0    # search checkpoint index within the sequence
    split_seq = ""

    # regular expressions to search for the motifs
    pattern_str = f'({"|".join([x.replace("N", "[ACGT]") for x in motifs])})'
    pattern = re.compile(pattern_str)
    seq_len = len(sequence)

    while position < seq_len:   # search till the end of the sequence

        mpos = pattern.search(sequence[position:])
        if mpos is None: # if no more motifs are found break
            break
        else:
            
            if sequence[position : mpos.start() + position] != "":
                 # sequence before the identified motif is non-empty add it as a split
                split_seq += f"{sequence[position:mpos.start()+position]} "
            
            # adding the identified motif pattern to the formatted string
            split_seq += f"{sequence[position + mpos.start(): position + mpos.start() + (mpos.end() - mpos.start())]} "
            
            # update the search check point
            position += mpos.start() + mpos.end() - mpos.start()

    # add the remaining stretch of the sequence
    if position < seq_len:
        split_seq += f"{sequence[position:]}"
    
    return split_seq

def read_bed(bed, format):
    """
    Reads the bed file and returns the locus information as a dictionary

    :param bed: the bed file
    :param format: the format of the bed file (strchive or pacbio)
    :return: the list of loci information
    """
    
    with open(bed) as fh:
        for line in fh:
            if line.startswith('#'):
                header = line.lstrip('#').strip().split('\t')
                ref_motif_index = header.index('reference_motif_reference_orientation')
                path_motif_index = header.index('pathogenic_motif_reference_orientation')
                locusid_index = header.index('id')
                continue
            line = line.strip().split('\t')

            if format.lower() == 'strchive':
                motifs = list(dict.fromkeys([x.strip().upper() for x in line[ref_motif_index].split(',')] + [x.strip().upper() for x in line[path_motif_index].split(',')]))
                locusid = line[locusid_index]
            elif format.lower() == 'pacbio' or format.lower() == 'trgt':
                annotations = line[3].split(';')
                motifs = [x.strip().upper() for x in annotations[1].split('=')[1].split(',')]
                locusid = annotations[0].split('=')[1]
            elif format.lower() == 'atarva':
                # motifs are in the 4th column, locusid in the 6th column
                motifs = [x.strip().upper() for x in line[3].split(',')]
                locusid = line[5]

            else:
                raise ValueError(f"Unknown format: {format}. Supported formats are 'strchive', 'pacbio', 'trgt', and 'atarva'.")

            locus = {
                'chrom': line[0],
                'start': int(line[1]),
                'end': int(line[2]),
                'motifs': motifs,
                'id': locusid
            }
            yield locus

def get_ref(fasta, ref_directory='.'):
    """
    Check if the reference genome can be used remotely or if file exists locally, if not download it
    Returns the reference genome object
    """
    # If directory doesn't exist, create it
    if not os.path.isdir(ref_directory):
        os.makedirs(ref_directory)
    # Try to use the ref genome remotely
    try:
        ref = pysam.Fastafile(fasta)
        return ref
    except:
        sys.stderr.write(f"Warning: Couldn't use open the fasta file remotely: {fasta}. Checking for a local copy\n")
        # Download file if it doesn't exist
    
    # Check if the reference genome file exists locally
    ref_file = os.path.basename(fasta)
    ref_path = ref_directory + ref_file
    if os.path.isfile(ref_path) or os.path.isfile(ref_path.replace('.gz', '')):
        try:
            ref = pysam.Fastafile(ref_path)
            return ref
        except:
            try:
                ref = pysam.Fastafile(ref_path.replace('.gz', ''))
                return ref
            except:
                pass
    # Download the reference genome file
    sys.stderr.write(f"Warning: couldn't find a local copy of reference genome file: {ref_path}. Downloading.\n")
    urllib.request.urlretrieve(fasta, ref_path)
    try:
        ref = pysam.Fastafile(ref_path)
        return ref
    except:
        if ref_path.endswith('.gz'):
            # try to unzip the file
            os.system(f"gzip -d {ref_path}")
            ref = pysam.Fastafile(ref_path.replace('.gz', ''))
            return ref
    # If the file is not found, raise an error
    raise FileNotFoundError(f"Reference genome file not found and couldn't be downloaded: {fasta}")

def compare_beds(bed1, bed2, name1, name2, ref, flank=10):
    for bed1_info, bed2_info in zip(read_bed(bed1, name1), read_bed(bed2, name2)):

        if bed1_info['id'] != bed2_info['id']:
            raise ValueError(f"IDs do not match: {bed1_info['id']} != {bed2_info['id']}")

        # choose the start of the left flank based on the smallest start coordinate
        lflank_start = bed1_info['start'] - flank if bed1_info['start'] < bed2_info['start'] else bed2_info['start'] - flank
        
        # choose the end of the right flank based on the largest end coordinate
        rflank_end = bed1_info['end'] + flank if bed1_info['end'] > bed2_info['end'] else bed2_info['end'] + flank
        
        # fetching the sequences
        # NOTE: the flank sequences are extracted from the start and end of bed1_info locus
        lflank            = ref.fetch(bed1_info['chrom'], lflank_start, bed1_info['start']).upper()
        rflank            = ref.fetch(bed1_info['chrom'], bed1_info['end'], rflank_end).upper()
        bed1_info_seq = ref.fetch(bed1_info['chrom'], bed1_info['start'], bed1_info['end']).upper()
        bed2_info_seq   = ""

        # flank sequence initialisers
        bed1_info_lflank = ""; bed1_info_rflank = ""
        bed2_info_lflank = ""; bed2_info_rflank = ""
            
        if bed1_info['start'] > bed2_info['start']:
            # if bed1_info start is greater than bed2_info start
            diff = bed1_info['start'] - bed2_info['start']
            # pull difference bases into the repeat for bed2_info and add gaps in the flank
            bed2_info_seq    = lflank[-diff:] + bed1_info_seq
            bed2_info_lflank = lflank[:-diff] + " "*diff
            # add difference gaps in the repeat for bed2_info and unchanged flank
            bed1_info_seq  = " "*diff + bed1_info_seq
            bed1_info_lflank = lflank
            
        elif bed1_info['start'] < bed2_info['start']:
            # if bed2_info start is greater than bed1_info start
            diff = bed2_info['start'] - bed1_info['start']
            # add difference gaps to bed2_info repeat sequence and different bases to flank
            bed2_info_seq = " "*diff + bed1_info_seq[diff:]
            bed2_info_lflank = lflank + bed1_info_seq[:diff]
            # add difference gaps to bed1_info flank 
            bed1_info_lflank = lflank[:-diff] + " "*diff
        else:
            # start coordinates are same
            bed2_info_lflank   = lflank
            bed1_info_lflank = lflank
            bed2_info_seq = bed1_info_seq
            
        if bed1_info['end'] < bed2_info['end']:
            # if bed1_info end is smaller than bed2_info end
            diff = bed2_info['end'] - bed1_info['end']
            # add the difference based from flank to bed2_info repeat and gaps in flank
            bed2_info_seq = bed2_info_seq + rflank[:diff]
            bed2_info_rflank = " "*diff + rflank[diff:]
            # add the difference gaps bed1_info repeat and unchanged flank
            bed1_info_seq += " "*diff
            bed1_info_rflank = rflank
        
        elif bed2_info['end'] < bed1_info['end']:
            # if bed2_info end is smalled than bed1_info end
            diff = bed1_info['end'] - bed2_info['end']
            # remove difference bases from bed2_info repeat and add gaps and add difference bases to the flank
            bed2_info_seq = bed2_info_seq[:-diff] + " "*diff
            bed2_info_rflank = bed2_info_seq[-diff:] + rflank
            # unchanged bed1_info repeat sequence and add difference gaps to bed1_info flank
            bed1_info_rflank = " "*diff + rflank
        else:
            # if the end coordinates are the same
            bed1_info_rflank = rflank
            bed2_info_rflank = rflank

        bed1_info_seq = split_repeat_sequence(bed1_info['motifs'], bed1_info_seq)
        bed2_info_seq = split_repeat_sequence(bed2_info['motifs'], bed2_info_seq)

        yield f"{bed1_info['id']}\n"
        yield f"{bed1_info['chrom']}\t{bed1_info['start']}\t{bed1_info['end']}\t{','.join(bed1_info['motifs'])}\t{name1}\n"
        yield f"{bed2_info['chrom']}\t{bed2_info['start']}\t{bed2_info['end']}\t{','.join(bed2_info['motifs'])}\t{name2}\n"
        yield f"{bed1_info_lflank}\t{bed1_info_seq}\t{bed1_info_rflank}\n"
        yield f"{bed2_info_lflank}\t{bed2_info_seq}\t{bed2_info_rflank}\n"
        yield "\n"

def main(beds: list[str], names: list[str], fasta: str, out: str, storage: str = '.', flank: int = 10):
    """
    Extracts the flanking sequences and the repeat sequence for tandem repeat loci from the reference genome

    # Accept one or more bedfiles and a corresponding list of names.
    :param list[str] beds: List of input bed files containing the tandem repeat loci
    :param list[str] names: List of names corresponding to each bed file
    :param fasta: Reference genome fasta file (required)
    :param output: Output file name for the extracted sequences
    :param storage: Reference genome storage directory (default: current directory)
    :param flank: Flanking sequence length to extract (default: 10)
    """
    ref = get_ref(fasta, storage)
    sys.stderr.write(f"Successfully accessed reference genome: {os.path.basename(fasta)}\n")

    # Assuming the bed files are sorted and contain the same loci
    with open(out, 'w') as outfh:
        
        for line in compare_beds(beds[0], beds[1], names[0], names[1], ref, flank):
            outfh.write(line)


    ref.close()
     
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main, cli_options='all')