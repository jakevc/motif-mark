
#### Develop a Python script to plot protein binding motifs on an image of
the an exon and flanking introns

**Input file:** fasta file with exons in all caps, and flanking intron sequence in lower-case.

**Output:** To scale SVG image (vector-based) with visaul representation of where the motifs fall on the exons, upstream, downstream, on the exon. Each motif has a different color.


Flow:

###argparse

	* command-line arguments
	* input file
	* motif file
	* RNA, DNA

Input fasta ->

	* parse_fasta entry():
		- returns the current fasta entry as an object with a corresponding coordinate object, including coordinates of exon, as well as the total length.
		- object has gene name attributes
		- parse_fasta entry will be a class with attributes.
		- entry.exon will return coordinates of the exon by searching for capitalized letters (start, end)
		- entry.motifs will return coordinates of the motifs

	* motif_mark():
		- Marks motif coordinates on sequence object.

	* Generate SVG():
		- For each gene in the fasta, write out a exons structure to scale.
		- then draw lines on exon to show where motifs occur
		- blue if color1 if rbfox, color2 if mbnl.
		- label 3' and 5' end based on class attributes













