import pysam
import re

fasta = pysam.Fastafile('../references/hg38/human_GRCh38_no_alt_analysis_set.fasta')


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


flank = 10
with open('/Users/avvarua/Downloads/STRchive-coordinate_adjustments.txt') as fh:
    for line in fh:
        line = line.strip().split('\t')

        strchive_chrom = line[3]
        strchive_start = int(line[4])
        strchive_end = int(line[5])
        strchive_motifs = [x.strip().strip('"') for x in line[1].split(',')]
        
        pacbio_chrom = line[8]
        pacbio_start = int(line[9])
        pacbio_end = int(line[10])
        pacbio_motifs = [x.strip().strip('"') for x in line[6].split(',')]

        # choose the start of the left flank based on the smallest start coordinate
        lflank_start = strchive_start - flank if strchive_start < pacbio_start else pacbio_start - flank
        
        # choose the end of the right flank based on the largest end coordinate
        rflank_end = strchive_end + flank if strchive_end > pacbio_end else pacbio_end + flank
        
        # fetching the sequences
        # NOTE: the flank sequences are extracted from the start and end of strchive locus
        lflank            = fasta.fetch(strchive_chrom, lflank_start, strchive_start)
        rflank            = fasta.fetch(strchive_chrom, strchive_end, rflank_end)
        strchive_seq = fasta.fetch(strchive_chrom, strchive_start, strchive_end)
        pacbio_seq   = ""

        # if there is coordinate/motif change between strchive and pacbio
        coord_change = False if line[11] == "FALSE" else True
        motif_change = False if line[12] == "FALSE" else True

        # flank sequence initialisers
        strchive_lflank = ""; strchive_rflank = ""
        pacbio_lflank = ""; pacbio_rflank = ""
        
        if coord_change:
            
            if strchive_start > pacbio_start:
                # if strchive start is greater than pacbio start
                diff = strchive_start - pacbio_start
                # pull difference bases into the repeat for pacbio and add gaps in the flank
                pacbio_seq    = lflank[-diff:] + strchive_seq
                pacbio_lflank = lflank[:-diff] + " "*diff
                # add difference gaps in the repeat for pacbio and unchanged flank
                strchive_seq  = " "*diff + strchive_seq
                strchive_lflank = lflank
                
            elif strchive_start < pacbio_start:
                # if pacbio start is greater than strchive start
                diff = pacbio_start - strchive_start
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
                
            if strchive_end < pacbio_end:
                # if strchive end is smaller than pacbio end
                diff = pacbio_end - strchive_end
                # add the difference based from flank to pacbio repeat and gaps in flank
                pacbio_seq = pacbio_seq + rflank[:diff]
                pacbio_rflank = " "*diff + rflank[diff:]
                # add the difference gaps strchive repeat and unchanged flank
                strchive_seq += " "*diff
                strchive_rflank = rflank
            
            elif pacbio_end < strchive_end:
                # if pacbio end is smalled than strchive end
                diff = strchive_end - pacbio_end
                # remove difference bases from pacbio repeat and add gaps and add difference bases to the flank
                pacbio_seq = pacbio_seq[:-diff] + " "*diff
                pacbio_rflank = pacbio_seq[-diff:] + rflank
                # unchanged strchive repeat sequence and add difference gaps to strchive flank
                strchive_rflank = " "*diff + rflank
            else:
                # if the end coordinates are the same
                strchive_rflank = rflank
                pacbio_rflank = rflank

            print(line[0])
            print(strchive_chrom, strchive_start, strchive_end, ','.join(strchive_motifs), 'STRCHIVE', sep='\t')
            print(pacbio_chrom,   pacbio_start,   pacbio_end, ','.join(pacbio_motifs), 'PACBIO', sep='\t')
            strchive_seq = split_repeat_sequence(strchive_motifs, strchive_seq)
            pacbio_seq = split_repeat_sequence(pacbio_motifs, pacbio_seq)
            print(strchive_lflank, strchive_seq, strchive_rflank)
            print(pacbio_lflank, pacbio_seq, pacbio_rflank)
            print()


fasta.close()
