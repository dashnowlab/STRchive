import pysam
import re
import sys

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
                motifs = [x.strip() for x in line[5].split(',')]
                locusid = line[3]
            elif format == 'pacbio':
                annotations = line[3].split(';')
                motifs = [x.strip() for x in annotations[1].split('=')[1].split(',')]
                locusid = annotations[0].split('=')[1]

            locus = {
                'chrom': line[0],
                'start': int(line[1]),
                'end': int(line[2]),
                'motifs': motifs,
                'id': locusid
            }
            yield locus

def main(bed1: str, bed2: str, fasta: str, out: str):
    """
    Extracts the flanking sequences and the repeat sequence for tandem repeat loci from the reference genome

    :param bed1: Input bed file containing the tandem repeat loci from STRchive
    :param bed2: Input bed file containing the tandem repeat loci from PacBio
    :param fasta: Reference genome fasta file
    :param output: Output file name for the extracted sequences
    """
    fasta = pysam.Fastafile(fasta)

    flank = 10

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
            lflank            = fasta.fetch(strchive['chrom'], lflank_start, strchive['start'])
            rflank            = fasta.fetch(strchive['chrom'], strchive['end'], rflank_end)
            strchive_seq = fasta.fetch(strchive['chrom'], strchive['start'], strchive['end'])
            pacbio_seq   = ""

            # if there is coordinate/motif change between strchive and pacbio
            if strchive['chrom'] != pacbio['chrom'] or strchive['start'] != pacbio['start'] or strchive['end'] != pacbio['end']:
                coord_change = True
            else:
                coord_change = False
            # if strchive['motifs'] != pacbio['motifs']:
            #     motif_change = True
            # else:
            #     motif_change = False

            # flank sequence initialisers
            strchive_lflank = ""; strchive_rflank = ""
            pacbio_lflank = ""; pacbio_rflank = ""
            
            if coord_change:
                
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
                outfh.write(f"{strchive['chrom']}\t{strchive['start']}\t{strchive['end']}\t{','.join(strchive['motifs'])}\tSTRCHIVE\n")
                outfh.write(f"{pacbio['chrom']}\t{pacbio['start']}\t{pacbio['end']}\t{','.join(pacbio['motifs'])}\tPACBIO\n")
                outfh.write(f"{strchive_lflank}\t{strchive_seq}\t{strchive_rflank}\n")
                outfh.write(f"{pacbio_lflank}\t{pacbio_seq}\t{pacbio_rflank}\n")
                outfh.write("\n")

    fasta.close()

main('/Users/dashnowh/Downloads/STRchive-disease-loci.v2.2.1.hg38.bed',
     '/Users/dashnowh/Downloads/STRchive-disease-loci.v2.2.1.hg38.TRGT.bed',
     'https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa',
     'output.txt')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import defopt
    defopt.run(main)