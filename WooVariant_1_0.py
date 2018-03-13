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
if len(template[0]) / 3 != len(template[1]):
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
        output = []
        count = 0
        codon = ''
        new_codon = [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}, {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                     {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}]
        for row in samfile.pileup(line.split('\t')[1].split(':')[1], 0, int(line.split('\t')[2].split(':')[
                1]), max_depth=2500000):
            codon += template[0][row.pos]
            counter = {}
            for read in row.pileups:
                if read.indel == 0 and read.query_position is not None:
                    if read.alignment.query_sequence[read.query_position] in ['A', 'C', 'G', 'T']:
                        position = read.alignment.query_sequence[read.query_position]
                        new_codon[count][position] += 1
                        if read.alignment.query_sequence[read.query_position] != template[0][row.pos]:
                            if position not in counter:
                                counter.update({position: 0})
                            counter[position] += 1
            for item in new_codon[count]:
                new_codon[count][item] = float(new_codon[count][item]) / float(row.n)
            if counter != {}:
                for item in counter:
                    outfh.write(line.split('\t')[1].split(':')[1] + '\t' + str(row.pos + 1) + '\t.\t' +
                                template[0][row.pos] + '\t' + item + '\t.\tPASS\tIS=' + str(counter[item]) + ';DP=' +
                                str(row.n) + ';AF=%.4f' % (float(counter[item]) / float(row.n)) + '\n')
            dnafh.write('\n' + str(row.pos + 1) + '\t' + template[0][row.pos])
            for item in ['A', 'C', 'G', 'T']:
                dnafh.write('\t')
                if item != template[0][row.pos] and item in counter:
                    dnafh.write(str(float(counter[item]) / float(row.n)))
            count += 1
            if count == 3:
                output_dict = {}
                for item in keys:
                    output_dict.update({item: 0.0})
                for item in codons:
                    if item != codon:
                        output_dict[codons[item]] += new_codon[0][item[0]] * new_codon[1][item[1]] * new_codon[2][
                            item[2]]
                output += [output_dict]
                count = 0
                codon = ''
                new_codon = [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}, {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
                             {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}]
        count = 0
        for item in output:
            aafh.write('\n' + str(count + 1) + '\t' + template[1][count] + '\t'.join([item[x] for x in keys]))
            count += 1
samfile.close()
outfh.close()
dnafh.close()
aafh.close()
print 'WooVariant finished!'
