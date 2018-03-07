# WooVariant_1.0
A tool to find all amino acid mutations in a sequencing result of a single gene

This is a new version for WooVariant. Compared to previous versions, this version could help to output handle sequencing result focusing on only one gene.
For example, if your sequencing result is only based on a gene, and you want to find frequence of every amino acid mutation on the gene, you could use this version.
Generally speaking, this program required at least three documents: sequencing result (by sorted BAM format), DNA sequence of gene (by fasta or txt format), and amino acid sequence (by fna or txt)
Notice that DNA sequence should be the complete operon for the amino acid sequence (i.e. begun with ATG and ended with stop codon), and only exon is supported.

Python packages "argparse" and "pysam" will be needed. The usage of this program is as following:

python WooVariant_1_0.py -i file/file.bam -d file/dna.fasta -a file/aa.fna -o file/name

The output format is csv in which both nucleotide and amino acid changes are reported. Two tables in csv format would also be reported, including one with nucleotide changes and their percentages in each position, the other with amino acid changes.
