# motif-mark

TO DO:

 - test on other sequences

MBNL and RBFOX regulate splicing in a similar manner, binding both up and downstream the alternatively spliced exon. The two factors have similar motifs, and are known to regulate each others binding sites. MBNL1 and RBFOX1 work together to regulate splicing events in transcripts that are regulated by both. To visualize the potential interaction of the two transcription factors the program written here generates and SVG image representation of sequence region, with the exon represented as a black box. A colored line is placed on the region of sequence where that transcription factor binding motif is found.

Multiple sequence regions and motifs can be represented in the same image by passing a multi entry sequences.fasta file, and motifs.txt file containing the binding motifs of interest.

The image generated using the provided example files and generated with the following command is shown below:

`./motifmark.py sequence.fasta motifs.txt -c rgb`



![example svg](https://github.com/jakevc/motif-mark/blob/master/example.svg)



