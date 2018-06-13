import argparse
import pysam


parser = argparse.ArgumentParser(description='WooVariants - A tool to find SNPs and indels in prokaryotic BAM file. '
                                             'Usage: python WooVariants.py -i file/name.sorted.bam -d file/dna.fasta -a'
                                             ' file/aa.faa -o file/name'
                                             'Output: file/name_SNP.vcf file/name_DNA.csv file/name_AA.csv'
                                             'Copyleft: Jacky Woo from ZHK Research team, iSynBio, SIAT.')
parser.add_argument('-a', '--aa', type=str, help='Input amino acid residue sequence of gene in form of faa or txt')
parser.add_argument('-d', '--dna', type=str, help='Input DNA sequence of gene in form of fasta, fna, fa or txt')
parser.add_argument('-i', '--input_file', type=str, help='Input bam file')
parser.add_argument('-o', '--output_file', type=str, help='Output vcf and csv file prefix')

args = parser.parse_args()
codons = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
          "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
keys = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', '*']
template = ['', '']
for count, infile in enumerate(args.dna, args.aa):
    infh = open(infile)
    for line in infh:
        if '>' not in line:
            template[count] += line.split('\n')[0]
    infh.close()
if (len(template[0]) / 3 != len(template[1])) or (len(template[0]) % 3 != 0):
    print 'NOTICE! Your input gene DNA is not matched with your input of amino acid residue'
output_file = args.output_file + '.vcf'
dna_file = args.output_file + '_DNA.csv'
aa_file = args.output_file + '_AA.csv'
outfh = open(output_file, 'w')
outfh.write('##fileformat=VCFv4.1\n##source=WooVariant\n'
            '##FILTER=<ID=IS,Number=1,Type=Integer,Description="Maximum number of reads supporting a mutation>"\n'
            '##FILTER=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">\n'
            '##FILTER=<ID=AF,Number=1,Type=Float,Description="Mutation abundance">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
samfile = pysam.AlignmentFile(args.input_file, 'rb')
dnafh = open(dna_file, 'w')
dnafh.write('Pos\tRef\tA\tC\tG\tT')
aafh = open(aa_file, 'w')
aafh.write('Pos\tRef\t' + '\t'.join(keys))
for line in samfile.text.split('\n'):
    if 'SQ' in line:
        for count, row in enumerate(samfile.pileup(line.split('\t')[1].split(':')[1], 0, int(line.split('\t')[2].split(':')[
                1]), max_depth=2500000)):
            dna_counter = {}
            if count % 3 == 0:
                total = 0
                aa_counter = {}
                for read in row.pileups:
                    if read.indel == 0 and read.query_position:
                        if (read.alignment.query_sequence[read.query_position:(read.query_position + 3)] in codons) and (
                                read.query_position <= (len(read)-3)):
                            aa_position = read.alignment.query_sequence[read.query_position:(read.query_position + 3)]
                            if aa_position != template[0][row.pos: (row.pos + 3)]:
                                if codons[aa_position] not in aa_counter:
                                    aa_counter.update({codons[aa_position]: 0})
                                aa_counter[codons[aa_position]] += 1
                            total += 1
                        if read.alignment.query_sequence[read.query_position] in ['A', 'C', 'G', 'T']:
                            dna_position = read.alignment.query_sequence[read.query_position]
                            if dna_position != template[0][row.pos]:
                                if dna_position not in dna_counter:
                                    dna_counter.update({dna_position: 0})
                            dna_counter[dna_position] += 1
                aafh.write('\n' + str(count / 3 + 1) + '\t' + template[1][count / 3])
                for item in keys:
                    aafh.write('\t')
                    if item in aa_counter:
                        aafh.write(str(float(aa_counter[item]) / float(total)))
                    else:
                        aafh.write('0')
                if dna_counter != {}:
                    for item in dna_counter:
                        outfh.write('T7RNAP\t' + str(row.pos + 1) + '\t.\t' + template[0][row.pos] + '\t' + item +
                            '\t.\tPASS\tIS=' + str(dna_counter[item]) + ';DP=' + str(row.n) + ';AF=%.4f' % (
                                    float(dna_counter[item]) / float(row.n)) + '\n')
                dnafh.write('\n' + str(row.pos + 1) + '\t' + template[0][row.pos])
                for item in ['A', 'C', 'G', 'T']:
                    dnafh.write('\t')
                    if item != template[0][row.pos] and item in dna_counter:
                        dnafh.write(str(float(dna_counter[item]) / float(row.n)))
                    else:
                        dnafh.write('0')
            else:
            # Following codes are the same as previous ones for nucleotide outputs
                for read in row.pileups:
                    if read.indel == 0 and read.query_position:
                        if read.alignment.query_sequence[read.query_position] in ['A', 'C', 'G', 'T']:
                            dna_position = read.alignment.query_sequence[read.query_position]
                            if dna_position != template[0][row.pos]:
                                if dna_position not in dna_counter:
                                    dna_counter.update({dna_position: 0})
                                dna_counter[dna_position] += 1
                if dna_counter != {}:
                    for item in dna_counter:
                        outfh.write('T7RNAP\t' + str(row.pos + 1) + '\t.\t' + template[0][row.pos] + '\t' + item +
                            '\t.\tPASS\tIS=' + str(dna_counter[item]) + ';DP=' + str(row.n) + ';AF=%.4f' % (
                                    float(dna_counter[item]) / float(row.n)) + '\n')
                dnafh.write('\n' + str(row.pos + 1) + '\t' + template[0][row.pos])
                for item in ['A', 'C', 'G', 'T']:
                    dnafh.write('\t')
                    if item != template[0][row.pos] and item in dna_counter:
                        dnafh.write(str(float(dna_counter[item]) / float(row.n)))
                    else:
                        dnafh.write('0')
samfile.close()
outfh.close()
dnafh.close()
aafh.close()
print 'WooVariant finished!'
