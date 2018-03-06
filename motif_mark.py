#! /usr/bin/env python

import argparse
import re
import cairo
import matplotlib.colors as mplc


def get_args():
    '''Define and return command line arguments'''
    parser = argparse.ArgumentParser(
        prog='motif mark',
        description='Determines positions of transcription factor \
        (rbfox, mbnl) binding motifs accross an exon and flanking \
        intron sequences, producing a to scale SVG representation \
        for each exon')
    parser.add_argument('sequences_fasta',
                        help='input fasta file (path)',
                        type=str)
    parser.add_argument('motifs_txt',
                        help='file of possible motifs',
                        type=str)
    parser.add_argument('-o', '--outfile',
                        help='outfile name',
                        type=str,
                        default='motifs.svg',
                        required=False)
    parser.add_argument('-c', '--colors',
                        help='matlab color spec of colors,\
                        one for each motif. Example: rgb, will show\
                        three motifs in red, green, and blue',
                        type=list,
                        required=True)
    return parser.parse_args()


def get_seq_regions(fastafile):
    '''returns dict of each exon region with fasta header as key'''
    with open(fastafile, 'r') as fh:

        # dict for exon regions
        region_dict = {}

        for line in fh:
            line = line.strip()

            # retain header as key
            if line[0] == '>':
                header = line
                seq = ''

            # seq as values
            elif line[0] != '>':
                seq += line
                region_dict[header] = seq

    return region_dict


def get_motif_list(motif_file):
    '''returns list of motifs in the motif file'''
    with open(motif_file, 'r') as mh:

        # list to hold motifs of interest
        motifs = []

        # loop over each line storing each motif in list
        for line in mh:
            motifs.append(line.strip())

        return motifs


def motif_patterns(motif_list):
    '''returns motif regex matches for any IUPAC nucleotide code'''

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
        '''return the start position of each motif'''

        self.motif_coords_list = []

        m = re.finditer(r'(?i)'+motif_match, self.seq)
        for i in m:
            self.motif_coords_list.append(str(i.start()))

        return self.motif_coords_list


def get_svg_arguments():
    '''call main script code'''

    # grab filename
    infile = get_args().sequences_fasta

    # grab motifile name
    motifile = get_args().motifs_txt

    # outfile name
    outfile = get_args().outfile

    # colors
    colors = get_args().colors

    # read in each fasta entry to a dictonary
    region_dict = get_seq_regions(infile)

    # get motifs
    motif_list = get_motif_list(motifile)

    # get regex matches for each motif
    motif_matches = motif_patterns(motif_list)

    exon_coords = []

    # list for all motif coordinates
    motif_coordinate_list = []

    # list of sequence lengths
    lengths = []

    for entry in region_dict:

        # create seq object for each fasta entry
        seq_obj = seq_region(region_dict[entry])

        # add to list of exon corrdinates
        exon_coords.append(seq_obj.exon_coord)

        # add to list of sequence lengths
        lengths.append(len(region_dict[entry]))

        # list of coordinates for each motif
        motif_entry_coords = []

        # loop over indexed motif patterns
        for pattern in motif_matches:

            # tuple of positions for each motif
            motif_entry_coords.append(
                tuple(seq_obj.find_motif_coords(pattern)))

        # list of coordinates for each motif, for each entry
        motif_coordinate_list.append(motif_entry_coords)

    return (motif_coordinate_list,
            exon_coords,
            lengths,
            motif_list,
            outfile,
            colors)


class draw_svg():
    '''draw class for motif mark svg creation'''
    def __init__(self, motifs, lengths, exons):
        self.motifs = motifs
        self.lengths = lengths
        self.exons = exons
        self.outfile = get_svg_arguments()[4]
        self.colors = get_svg_arguments()[5]

    def setup_surface(self):
        '''draws a surface for the number of fasta entries'''

        # dimensions
        n = len(self.exons)
        self.width = 100 * n
        self.height = (100 * n) + 25 * n
        self.svg_len = max(self.lengths)

        # surface
        self.surface = cairo.SVGSurface(self.outfile, 300, self.height)
        self.context = cairo.Context(self.surface)
        self.context.scale(300, 300)

        # centers
        self.centers = []
        top = 1
        self.centers.append(top/(n+1) + 0.25)

        # create the right number of centers for sequences
        while top < n:
            top += 1
            self.centers.append(top/(n+1) + 0.25)

    def draw_seq_regions(self):
        # sequence lines across figure
        for l, c in zip(self.lengths, self.centers):
            x, x2 = 0, l/self.svg_len
            y, y2 = c, c
            self.context.move_to(x, y)
            self.context.set_line_width(0.01)
            self.context.line_to(x2, y2)

        self.context.stroke()

    def draw_exons(self):
        '''draws exon from a string of coordinates "start:stop"'''

        for seqlen, ecenter, exon in zip(
                self.lengths, self.centers, self.exons):

            # set exon color as black with mpl function
            r, g, b, a = mplc.to_rgba('black')
            self.context.set_source_rgba(r, g, b, a)

            # Rectangle(x0, y0, x1, y1)
            self.context.rectangle(exon[0]/self.svg_len,
                                   (ecenter-0.02),
                                   (exon[1]-exon[0])/self.svg_len, 0.04)

        self.context.stroke()

    def draw_motifs(self):
        '''draws motifs of the desired color on the sequence region'''

        for mot, center in zip(self.motifs, self.centers):

            # for each attribute of each seq region
            for m, seqlen, col in zip(mot, self.lengths, self.colors):

                # use mpl color function to get the color you want
                r, g, b, a = mplc.to_rgba(col)

                # set that color as the context
                self.context.set_source_rgba(r, g, b, a)

                for pos in m:

                    motif = int(pos)

                    # scale to svg_len
                    x = motif/self.svg_len

                    # draw around center line
                    y = center + 20/self.svg_len
                    y1 = center - 20/self.svg_len

                    # create line
                    self.context.move_to(x, y)
                    self.context.line_to(x, y1)

                    self.context.set_line_width(0.005)

                    self.context.stroke()

    def draw_legend(self):

        # black legend rectangle
        self.context.set_source_rgba(0, 0, 0, 1)
        self.context.set_line_width(0.004)
        self.context.rectangle(0.004, 0.004, 0.5, 0.25)

        self.context.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,
                                      cairo.FONT_WEIGHT_BOLD)
        self.context.set_font_size(0.05)

        # x, y coordinate of text
        x = (0.5 - 0.004)/4
        y = 0.05
        self.context.move_to(x, y)
        self.context.show_text("Motif Legend")
        self.context.stroke()

        # get motif text from file
        motif_txt = get_svg_arguments()[3]

        # draw each motif with its color in the rectangle
        for i, mtext in enumerate(motif_txt):
            r, b, g, a = mplc.to_rgba(self.colors[i])
            self.context.set_source_rgba(r, b, g, a)

            self.context.move_to(x, 0.10 + 0.05*i)
            self.context.show_text(mtext)
            self.context.stroke()

    def finish_svg(self):
        # start svg
        svg = draw_svg(self.motifs, self.lengths, self.exons)

        # setup svg surface
        svg.setup_surface()

        # draw the scaled regions for each entry
        svg.draw_seq_regions()

        # draw exons
        svg.draw_exons()

        svg.draw_motifs()

        svg.draw_legend()

        # finish the surface!
        svg.surface.finish()


def main():

    # get lists for svg_draw
    motifs, exons, lengths = get_svg_arguments()[0:3]

    # draw the svg
    draw_svg(motifs, lengths, exons).finish_svg()


if __name__ == '__main__':
    main()
