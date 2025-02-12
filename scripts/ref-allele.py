import pysam
import re
import sys
import os
import urllib.request

def split_repeat_sequence(motifs, sequence):
    """
    Splits repeat sequence based on the defined motifs

    Args:
        motifs  : the list of motifs defined for the locus
        sequence: sequence of the locus
    
    Returns:
        the formatted sequence string split based on the motifs
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

    Args:
        bed: the bed file
        format: the format of the bed file (strchive or pacbio)
    
    Returns:
        the list of loci information
    """
    with open(bed) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')

            if format == 'strchive':
                motifs = [x.strip().upper() for x in line[5].split(',')]
                locusid = line[3]
            elif format == 'pacbio':
                annotations = line[3].split(';')
                motifs = [x.strip().upper() for x in annotations[1].split('=')[1].split(',')]
                locusid = annotations[0].split('=')[1]

            locus = {
                'chrom': line[0],
                'start': int(line[1]),
                'end': int(line[2]),
                'motifs': motifs,
                'id': locusid
            }
            yield locus

def get_ref(fasta):
    """
    Check if the reference genome can be used remotely or if file exists locally, if not download it
    Returns the reference genome object
    """
    # Try to use the ref genome remotely
    try:
        ref = pysam.Fastafile(fasta)
        return ref
    except:
        sys.stderr.write(f"Warning: Couldn't use open the fasta file remotely: {fasta}. Checking for a local copy.\n")
        # Download file if it doesn't exist
    
    # Check if the reference genome file exists locally
    ref_file = os.path.basename(fasta)
    if os.path.isfile(ref_file) or os.path.isfile(ref_file.replace('.gz', '')):
        try:
            ref = pysam.Fastafile(ref_file)
            return ref
        except:
            try:
                ref = pysam.Fastafile(ref_file.replace('.gz', ''))
                return ref
            except:
                pass
    # Download the reference genome file
    sys.stderr.write(f"Warning: couldn't find a local copy of reference genome file: {fasta}. Downloading.\n")
    urllib.request.urlretrieve(fasta, ref_file)
    try:
        ref = pysam.Fastafile(ref_file)
        return ref
    except:
        if ref_file.endswith('.gz'):
            # try to unzip the file
            os.system(f"gzip -d {ref_file}")
            ref_file = ref_file.replace('.gz', '')
            ref = pysam.Fastafile(ref_file)
            return ref
    # If the file is not found, raise an error
    raise FileNotFoundError(f"Reference genome file not found and couldn't be downloaded: {fasta}")

def main(bed1: str, bed2: str, fasta: str, out: str, ref: str = None, flank: int = 10):
    """
    Extracts the flanking sequences and the repeat sequence for tandem repeat loci from the reference genome

    :param bed1: Input bed file containing the tandem repeat loci from STRchive
    :param bed2: Input bed file containing the tandem repeat loci from PacBio
    :param fasta: Reference genome fasta file
    :param output: Output file name for the extracted sequences
    :param ref: Reference genome name (hg19, hg38, T2T-chm13)
    :param flank: Flanking sequence length to extract (default: 10)
    """
    ref = get_ref(fasta)

    # Assuming the bed files are sorted and contain the same loci
    with open(out, 'w') as outfh:
        for strchive, pacbio in zip(read_bed(bed1, 'strchive'), read_bed(bed2, 'pacbio')):
            if strchive['id'] != pacbio['id']:
                raise ValueError(f"IDs do not match: {strchive['id']} != {pacbio['id']}")

            # choose the start of the left flank based on the smallest start coordinate
            lflank_start = strchive['start'] - flank if strchive['start'] < pacbio['start'] else pacbio['start'] - flank
            
            # choose the end of the right flank based on the largest end coordinate
            rflank_end = strchive['end'] + flank if strchive['end'] > pacbio['end'] else pacbio['end'] + flank
            
            # fetching the sequences
            # NOTE: the flank sequences are extracted from the start and end of strchive locus
            lflank            = ref.fetch(strchive['chrom'], lflank_start, strchive['start']).upper()
            rflank            = ref.fetch(strchive['chrom'], strchive['end'], rflank_end).upper()
            strchive_seq = ref.fetch(strchive['chrom'], strchive['start'], strchive['end']).upper()
            pacbio_seq   = ""

            # flank sequence initialisers
            strchive_lflank = ""; strchive_rflank = ""
            pacbio_lflank = ""; pacbio_rflank = ""
                
            if strchive['start'] > pacbio['start']:
                # if strchive start is greater than pacbio start
                diff = strchive['start'] - pacbio['start']
                # pull difference bases into the repeat for pacbio and add gaps in the flank
                pacbio_seq    = lflank[-diff:] + strchive_seq
                pacbio_lflank = lflank[:-diff] + " "*diff
                # add difference gaps in the repeat for pacbio and unchanged flank
                strchive_seq  = " "*diff + strchive_seq
                strchive_lflank = lflank
                
            elif strchive['start'] < pacbio['start']:
                # if pacbio start is greater than strchive start
                diff = pacbio['start'] - strchive['start']
                # add difference gaps to pacbio repeat sequence and different bases to flank
                pacbio_seq = " "*diff + strchive_seq[diff:]
                pacbio_lflank = lflank + strchive_seq[:diff]
                # add difference gaps to strchive flank 
                strchive_lflank = lflank[:-diff] + " "*diff
            else:
                # start coordinates are same
                pacbio_lflank   = lflank
                strchive_lflank = lflank
                pacbio_seq = strchive_seq
                
            if strchive['end'] < pacbio['end']:
                # if strchive end is smaller than pacbio end
                diff = pacbio['end'] - strchive['end']
                # add the difference based from flank to pacbio repeat and gaps in flank
                pacbio_seq = pacbio_seq + rflank[:diff]
                pacbio_rflank = " "*diff + rflank[diff:]
                # add the difference gaps strchive repeat and unchanged flank
                strchive_seq += " "*diff
                strchive_rflank = rflank
            
            elif pacbio['end'] < strchive['end']:
                # if pacbio end is smalled than strchive end
                diff = strchive['end'] - pacbio['end']
                # remove difference bases from pacbio repeat and add gaps and add difference bases to the flank
                pacbio_seq = pacbio_seq[:-diff] + " "*diff
                pacbio_rflank = pacbio_seq[-diff:] + rflank
                # unchanged strchive repeat sequence and add difference gaps to strchive flank
                strchive_rflank = " "*diff + rflank
            else:
                # if the end coordinates are the same
                strchive_rflank = rflank
                pacbio_rflank = rflank


            strchive_seq = split_repeat_sequence(strchive['motifs'], strchive_seq)
            pacbio_seq = split_repeat_sequence(pacbio['motifs'], pacbio_seq)
            
            outfh.write(f"{strchive['id']}\n")
            outfh.write(f"{strchive['chrom']}\t{strchive['start']}\t{strchive['end']}\t{','.join(strchive['motifs'])}\tSTRchive\n")
            outfh.write(f"{pacbio['chrom']}\t{pacbio['start']}\t{pacbio['end']}\t{','.join(pacbio['motifs'])}\tTRGT\n")
            outfh.write(f"{strchive_lflank}\t{strchive_seq}\t{strchive_rflank}\n")
            outfh.write(f"{pacbio_lflank}\t{pacbio_seq}\t{pacbio_rflank}\n")
            outfh.write("\n")

    ref.close()

# To do
# Check the coordintes in the STRchive json (or bed?) for each reference genome
# Also check the TRGT bed file for the coordinates


# main('/Users/dashnowh/Downloads/STRchive-disease-loci.v2.2.1.hg38.bed',
#      '/Users/dashnowh/Downloads/STRchive-disease-loci.v2.2.1.hg38.TRGT.bed',
#      'https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa',
#      'ref-alleles.hg38.txt')

# main('/Users/dashnowh/Documents/git/STRchive/data/STRchive-disease-loci.hg19.bed',
#      '/Users/dashnowh/Documents/git/STRchive/data/STRchive-disease-loci.hg19.TRGT.bed',
#      'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
#      'ref-alleles.hg19.txt')

# main('/Users/dashnowh/Documents/git/STRchive/data/STRchive-disease-loci.T2T-chm13.bed',
#      '/Users/dashnowh/Documents/git/STRchive/data/STRchive-disease-loci.T2T-chm13.TRGT.bed',
#      'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz',
#      'ref-alleles.T2T-chm13.txt')

# Next: run this from snakemake instead of here. 
# Maybe host these genomes somewhere so I don't have to download them?
# https://zenodo.org/records/14853928/files/chm13v2.0_maskedY_rCRS.fa.gz
     
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)