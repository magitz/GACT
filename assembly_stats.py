#!/usr/bin/python

"""
Modified from https://github.com/HuffordLab/GenomeQC
  project to pull out as stand alone script
  
This script calculates many basic length metrics 
  for the given input genome assembly.

Edited by Matt Gitzendanner
Versions: 
  1.0 - April 6, 2021 - Refactor to write descriptions 
                        and reduce redundancy
  2.0 - April 13, 2021 - Refactor to compare multiple input assemblies
          - Break several things into functions
          - Improve stats calculation efficiency

"""
__version__ = '2.0'

from Bio import SeqIO
import sys
import statistics
import numpy as np
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Genome assembly stats script.')
parser.add_argument('-i', '--infile', required=True, nargs='+'
                    help='Input genome assembly file(s) (fasta format). Add multiple assemblies for comparison among them: e.g. -i method1.fa method2.fa method3.fa')
parser.add_argument('-o','--outfile', help='Output filename')
parser.add_argument('-s', '--size', required=True,
                    help='Genome size estimate in MBp')
parser.add_argument('-n', '--num_longest', default = 10,
                    help='Print the length of the n longest scaffolds. Default=10')
#parser.add_argument('-f', '--format', default=1,
#                    help='Output format: 1: Full text file, 2: Values only text file, 3: STDOUT. Default: 1')


def get_genome_stats(infile):
  """Calculate the statistics for an input fasta formatted genome file.

  infile should be a fasta formatted file.
  Returns a dictionary of values for that genome.
  """

  records = list(SeqIO.parse(infile, "fasta"))
  genome_dict={}
  
  genome_dict[number_of_scaffolds] = len(records)
  genome_dict[len_seq] = [len(rec) for rec in records]
  genome_dict[scaffold_lengths] = pd.DataFrame(len_seq)
  genome_dict[total_size_scaffolds] = sum(len_seq)

  genome_dict[total_scaffold_length_percentage_genome_size] = ((total_size_scaffolds/(estimated_genome_size*1000000))*100)
    
  genome_dict[seq_greater_1k] = scaffold_lengths[scaffold_lengths > 1000].count() 
  genome_dict[seq_greater_10k] = scaffold_lengths[scaffold_lengths > 10000].count()
  genome_dict[seq_greater_25k] = scaffold_lengths[scaffold_lengths > 25000].count()
  genome_dict[seq_greater_100k] = scaffold_lengths[scaffold_lengths > 100000].count()
  genome_dict[seq_greater_1M] = scaffold_lengths[scaffold_lengths > 1000000].count()
  genome_dict[seq_greater_10M] = scaffold_lengths[scaffold_lengths > 1000000].count()

  genome_dict[sorted_len] = sorted(len_seq, reverse=True)
  genome_dict[sum_sorted_length] = sum(sorted_len)
  genome_dict[half_length] = sum_sorted_length/2.0

  #calculates N50 and L50 values
  testSum = 0
  genome_dict[N50] = 0
  genome_dict[L50] = 0
  i = 0
  for con in genome_dict[sorted_len]:
      testSum += con
      i += 1
      if half_length < testSum:
        genome_dict[N50] = con
        genome_dict[L50] = i
        break

  #calculates NG50 and LG50 values
  half_genome = (estimated_genome_size*1000000)/2.0
  testSumNG50 = 0
  genome_dict[NG50] = 0
  genome_dict[LG50] = 0
  i = 0
  for conNG50 in sorted_len:
      testSumNG50 += conNG50
      i += 1
      if  half_genome < testSumNG50:
        genome_dict[NG50] = conNG50
        [LG50] = i
        break

  #calculates A,C,G,T,N percentages: No real need for both A&T and G&C? Comment out for now.
  genome_dict[counterA] = 0
  #genome_dict[counterC] = 0
  genome_dict[counterG] = 0
  #genome_dict[counterT] = 0
  genome_dict[counterN] = 0
  for record in records:
      genome_dict[counterA] += record.seq.count('A') 
      #genome_dict[counterC] += record.seq.count('C') 
      genome_dict[counterG] += record.seq.count('G') 
      #genome_dict[counterT] += record.seq.count('T') 
      genome_dict[counterN] += record.seq.count('N') 

  return genome_dict


def write_output_stats(genome_stats, outputfile):

  OUT = open(outputfile, 'w')

  OUT.write('Genome:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome])
  OUT.write('\n')

  OUT.write('Number of scaffolds:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][number_of_scaffolds])
  OUT.write('\n')

  OUT.write('Total sum of scaffold lengths:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][total_size_scaffolds])
  OUT.write('\n')

  OUT.write('Percent of genome size:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][total_scaffold_length_percentage_genome_size])
  OUT.write('\n')

  OUT.write('Longest scaffold:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(max(genome_stats[genome][len_seq]))
  OUT.write('\n')


  OUT.write('Shortest scaffold:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(max(genome_stats[genome][len_seq]))
  OUT.write('\n\n')

  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')

  OUT.write('Total no. scaffolds over 1KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_1k]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 1KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_1k]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 1KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_1k])/(estimated_genome_size*1000000))*100})
  OUT.write('\n')

  OUT.write('\n')
   OUT.write('Total no. scaffolds over 10KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_10k]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 10KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_10k]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 10KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_10k])/(estimated_genome_size*1000000))*100})
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 25KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_25k]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 25KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_25k]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 25KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_25k])/(estimated_genome_size*1000000))*100})
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 100KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_100k]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 100KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_100k]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 100KBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_100k])/(estimated_genome_size*1000000))*100})
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 1MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_1M]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 1MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_1M]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 1MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_1M])/(estimated_genome_size*1000000))*100})
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 10MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(len(genome_stats[genome][seq_greater_10M]))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 10MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_10M]))
  OUT.write('\n')

  OUT.write(f'Percent genome over 10MBp:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(sum(genome_stats[genome][seq_greater_10M])/(estimated_genome_size*1000000))*100})
  OUT.write('\n\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')
  
  OUT.write('N50:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][N50])
  OUT.write('\n')

  OUT.write('L50:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][L50])
  OUT.write('\n')

  OUT.write('NG50:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][NG50])
  OUT.write('\n')

  OUT.write('LG50:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write(genome_stats[genome][LG50])
  OUT.write('\n\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')

  OUT.write('%AT:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write((genome_stats[genome][counterA]/total_size_scaffolds)*100)
  OUT.write('\n')

  OUT.write('%GC:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write((genome_stats[genome][counterG]/total_size_scaffolds)*100)
  OUT.write('\n')

  OUT.write('%N:')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write((genome_stats[genome][counterN]/total_size_scaffolds)*100)
  OUT.write('\n')

  OUT.write(f'\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')
  # Print longest n scaffold lengths
  OUT.write(f'Logest {args.num_longest} scaffolds:')
  OUT.write(f'')
  for genome in genome_stats:
    OUT.write('\t')
    OUT.write((genome_stats[genome][scaffold_lengths.sort_values(0,ascending=False).head(n=args.num_longest)])
  OUT.write('\n')

  OUT.close()


def main():

  args = parser.parse_args()
  estimated_genome_size = float(args.size)

  print(f"Got {len(args.infile)} genome(s) to compare.")

      
  # Genome stats are stored in a dict of dicts, one for each genome.
  genomes_dict = {}
  append = 1

  # For each genome get the stats.
  for genome in args.infile:
    genome_name = os.path.split(genome)[1]

    if genome_name in genomes_dict:
      print(f"Multiple genome assemblies found with name {genome_name}, appending numbers to subsequent genomes")
      genome_name = genome_name + "_" + str(append)
      append += 1

    print(f"Getting statistics for {genome_name}.")
    genomes_dict[genome_name] = get_genome_stats(genome)

  print(f"Done getting stats for {len(args.infile)} genomes. \nSummarizing data.")

  write_output_stats(genomes_dict, args.outfile)



if __name__ == __main__:
  main()