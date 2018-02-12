#! /usr/bin/env python3.6

import argparse
import cairo
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
                        type=int,
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


def motif_coords(motif, string):
    'returns all case-insensitive matches of the motif in a string'

    motif_coords_list = []

    # find all case-insensitive matches of the motif in a string
    m = re.finditer(r'(?i)'+motif, string)
    for i in m:
        motif_coords_list.append(str(i.start())+':'+str(i.end()))

    return motif_coords_list


class seq_region:
    '''Class for each exon region'''
    def __init__(self, seq):
        self.seq = seq

    def exon(self, string):

        # initialize lists for coordinates and seq
        coords = []
        exon = []

        # loop over seq to fill in lists
        for i, char in enumerate(string):

            # exons are uppercase
            if char.isupper():
                coords.append(i)
                exon.append(char)

        # return exon sequence
        self.exon = exon

        # return coordinates of exon
        self.coords = (min(coords), max(coords))

    def find_motif_coords(self, motif):

        motif_coords_list = []

        m = re.finditer(r'(?i)'+motif, self.seq)
        for i in m:
            motif_coords_list.append(str(i.start())+':'+str(i.end()))

        self.motif_coords_list = motif_coords_list


# read in each fasta entry to a dictonary
region_dct = get_seq_regions(infile)



