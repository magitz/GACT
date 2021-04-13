#!/usr/bin/python

"""
Modified from https://github.com/HuffordLab/GenomeQC
  project to pull out as stand alone script
  
This script calculates all the basic length metrics 
  for the given input genome assembly.

Edited by Matt Gitzendanner
Versions: 
  1.0 - April 6, 2021 - Refactor to write descriptions 
                        and reduce redundancy
"""

from Bio import SeqIO
import sys
import statistics
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Genome assembly stats script.')
parser.add_argument('-i', '--infile', required=True,
                    help='Input genome assembly file (fasta format)')
parser.add_argument('-o','--outfile', help='Output filename')
parser.add_argument('-s', '--size', required=True,
                    help='Genome size estimate in MBp')
parser.add_argument('-n', '--num_longest', default = 10,
                    help='Print the length of the n longest scaffolds. Default=10')
#parser.add_argument('-f', '--format', default=1,
#                    help='Output format: 1: Full text file, 2: Values only text file, 3: STDOUT. Default: 1')

args = parser.parse_args()

inputfile = args.infile
outputfile = args.outfile
estimated_genome_size = float(args.size)
n_longest_to_print = int(args.num_longest)

OUT = open(outputfile, 'w')

records = list(SeqIO.parse(inputfile, "fasta"))

number_of_scaffolds = len(records)

OUT.write(f'Number of scaffolds:\t{number_of_scaffolds}\n')


len_seq = [len(rec) for rec in records]
scaffold_lengths = pd.DataFrame(len_seq)
total_size_scaffolds = sum(len_seq)

OUT.write(f'Total sum of scaffold lengths:\t{total_size_scaffolds}\n')

total_scaffold_length_percentage_genome_size = ((total_size_scaffolds/(estimated_genome_size*1000000))*100)

OUT.write(f'Percent of genome size:\t{total_scaffold_length_percentage_genome_size}\n')
OUT.write(f'Longest scaffold:\t{max(len_seq)}\n')
OUT.write(f'Shortest scaffold:\t{min(len_seq)}\n')
OUT.write('\n')
OUT.write('--------------------------------------------------------------------\n')
OUT.write('\n')

seq_greater_1k = sorted(i for i in len_seq if i>1000)
OUT.write(f'Total no. scaffolds over 1KBp:\t{len(seq_greater_1k)}\n')
OUT.write(f'Sum of scaffold lengths over 1KBp:\t{sum(seq_greater_1k)}\n')
OUT.write(f'Percent genome over 1KBp:\t{(sum(seq_greater_1k)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
seq_greater_10k = sorted(i for i in len_seq if i>10000)
OUT.write(f'Total no. scaffolds over 10KBp:\t{len(seq_greater_10k)}\n')
OUT.write(f'Sum of scaffold lengths over 10KBp:\t{sum(seq_greater_10k)}\n')
OUT.write(f'Percent genome over 10KBp:\t{(sum(seq_greater_10k)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
seq_greater_25k = sorted(i for i in len_seq if i>25000)
OUT.write(f'Total no. scaffolds over 25KBp:\t{len(seq_greater_25k)}\n')
OUT.write(f'Sum of scaffold lengths over 25KBp:\t{sum(seq_greater_25k)}\n')
OUT.write(f'Percent genome over 25KBp:\t{(sum(seq_greater_25k)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
seq_greater_100k = sorted(i for i in len_seq if i>100000)
OUT.write(f'Total no. scaffolds over 100KBp:\t{len(seq_greater_100k)}\n')
OUT.write(f'Sum of scaffold lengths over 100KBp:\t{sum(seq_greater_100k)}\n')
OUT.write(f'Percent genome over 100KBp:\t{(sum(seq_greater_100k)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
seq_greater_1M = sorted(i for i in len_seq if i>1000000)
OUT.write(f'Total no. scaffolds over 1MBp:\t{len(seq_greater_1M)}\n')
OUT.write(f'Sum of scaffold lengths over 1MBp:\t{sum(seq_greater_1M)}\n')
OUT.write(f'Percent genome over MBp:\t{(sum(seq_greater_1M)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
seq_greater_10M = sorted(i for i in len_seq if i>10000000)
OUT.write(f'Total no. scaffolds over 10MBp:\t{len(seq_greater_10M)}\n')
OUT.write(f'Sum of scaffold lengths over 10MBp:\t{sum(seq_greater_10M)}\n')
OUT.write(f'Percent genome over 10MBp:\t{(sum(seq_greater_10M)/(estimated_genome_size*1000000))*100}\n')
OUT.write('\n')
OUT.write('--------------------------------------------------------------------\n')
OUT.write('\n')
#calculates N50 and L50 values

sorted_len = sorted(len_seq, reverse=True)
sum_sorted_length = sum(sorted_len)
half_length = sum_sorted_length/2.0

testSum = 0
N50 = 0
L50 = 0
i = 0

for con in sorted_len:
    testSum += con
    i += 1
    if half_length < testSum:
       N50 = con
       L50 = i
       break

OUT.write(f'N50:\t{N50}\n')
OUT.write(f'L50:\t{L50}\n')


#calculates NG50 and LG50 values
half_genome = (estimated_genome_size*1000000)/2.0

testSumNG50 = 0
NG50 = 0
LG50 = 0
i = 0

for conNG50 in sorted_len:
    testSumNG50 += conNG50
    i += 1
    if  half_genome < testSumNG50:
       NG50 = conNG50
       LG50 = i
       break

OUT.write(f'NG50:\t{NG50}\n')
OUT.write(f'LG50:\t{LG50}\n')
OUT.write('\n')
OUT.write('--------------------------------------------------------------------\n')
OUT.write('\n')
#calculates A,C,G,T,N percentages
counterA = 0
counterC = 0
counterG = 0
counterT = 0
counterN = 0
for record in records:
    counterA += record.seq.count('A') 
    counterC += record.seq.count('C') 
    counterG += record.seq.count('G') 
    counterT += record.seq.count('T') 
    counterN += record.seq.count('N') 

OUT.write(f'%A:\t{((counterA)/total_size_scaffolds)*100}\n')
OUT.write(f'%C:\t{((counterC)/total_size_scaffolds)*100}\n')
OUT.write(f'%G:\t{((counterG)/total_size_scaffolds)*100}\n')
OUT.write(f'%T:\t{((counterT)/total_size_scaffolds)*100}\n')
OUT.write(f'%N:\t{((counterN)/total_size_scaffolds)*100}\n')
OUT.write(f'\n')
OUT.write('--------------------------------------------------------------------\n')
OUT.write('\n')
# Print longest n scaffold lengths
OUT.write(f'Logest {n_longest_to_print} scaffolds:')
OUT.write(f'{scaffold_lengths.sort_values(0,ascending=False).head(n=n_longest_to_print)}')

OUT.close()
