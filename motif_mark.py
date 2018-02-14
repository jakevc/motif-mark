#! /usr/bin/env python3.6

import argparse
# import cairo
import re


def get_args():
    '''Define and return command line arguments'''
    parser = argparse.ArgumentParser(
        prog='motif mark',
        description='Determines positions of transcription factor \
        (rbfox, mbnl) binding motifs accross an exon and flanking \
        intron sequences, producing a to scale SVG representation \
        for each exon')
    parser.add_argument('fasta_file',
                        help='input fasta file (path)',
                        type=str)
    parser.add_argument('motifs_file',
                        type=str,
                        help='Specify the k most common seqs.\
                        Default: 10',
                        default=10,
                        nargs='?')
    return parser.parse_args()

# grab filename
infile = get_args().fasta_file

# grab motifile name
motifile = get_args().motifs_file


def get_seq_regions(fastafile):
    with open(fastafile, 'r') as fh:

        # dict for exon regions
        exon_dict = {}

        for line in fh:
            line = line.strip()

            # retain header as key
            if line[0] == '>':
                header = line
                seq = ''

            # seq as values
            elif line[0] != '>':
                seq += line
                exon_dict[header] = seq

    return exon_dict


def get_motif_list(motif_file):
    with open(motif_file, 'r') as mh:

        # list to hold motifs of interest
        motifs = []

        # loop over each line storing each motif in list
        for line in mh:
            motifs.append(line.strip())

        return motifs


def motif_patterns(motif_list):
    'returns motif regex matches for any IUPAC nucleotide code'

    # list of motif matches
    motif_match_list = []

    # dict of IUPAC nucleotide codes
    IUPAC_dict = {'y': '(C|T|U)',
                  'r': '(A|G)',
                  's': '(C|G)',
                  'w': '(A|T|U)',
                  'k': '(G|T|U)',
                  'm': '(C|A)',
                  'b': '(C|G|T|U)',
                  'd': '(A|G|T|U)',
                  'v': '(A|C|G)',
                  'h': '(A|C|T|U)'
                  }

    # for each motif
    for tif in motif_list:

        # for each nucleotide code
        for key in IUPAC_dict:
            m = re.finditer(r'(?i)'+key, tif)

            # replace the code with a regex match for the possible nucs.
            for i in m:
                tif = tif.replace(i.group(), IUPAC_dict[key])

        # add to motif match list
        motif_match_list.append(tif)
    return motif_match_list


class seq_region:
    '''Class for each exon region'''
    def __init__(self, seq):

        # the sequence for each entry
        self.seq = seq

        # initialize lists for coordinates and seq
        self.coords = []
        self.exon = []

        # loop over seq to fill in lists
        for i, char in enumerate(self.seq):

            # exons are uppercase
            if char.isupper():
                self.coords.append(i)
                self.exon.append(char)

        # return exon sequence // might not need this
        self.exon = ''.join(self.exon)

        # return coordinates of exon
        self.exon_coord = (min(self.coords), max(self.coords))

    def find_motif_coords(self, motif_match):

        self.motif_coords_list = []

        m = re.finditer(r'(?i)'+motif_match, self.seq)
        for i in m:
            self.motif_coords_list.append(str(i.start())+':'+str(i.end()))

        return self.motif_coords_list


def main():
    # read in each fasta entry to a dictonary
    region_dict = get_seq_regions(infile)

    # get motifs
    motif_list = get_motif_list(motifile)

    # get regex matches for each motif
    motif_matches = motif_patterns(motif_list)

    # place to put motif coordinate dicts for each entry
    gene_motif_dict = {}

    # for each fasta entry make seq_region object from the sequence
    for entry in region_dict:

        seq_obj = seq_region(region_dict[entry])

        # place to put coordinates for each motif
        motif_coordinate_dict = {}

        # loop over indexed motif patterns
        for i, pattern in enumerate(motif_matches):

            # determine coordinates of motif occurances in the sequence
            motif_coordinate_dict[motif_list[i]] = \
                                    seq_obj.find_motif_coords(pattern)

            # name them and add data to coordinate dict
            gene_motif_dict[entry] = motif_coordinate_dict

    return gene_motif_dict


if __name__ == main():
    main()

