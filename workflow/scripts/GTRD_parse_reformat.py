#!/usr/bin/python

'''

GTRD_parse_closest.py

This script takes the raw .fasta created by prvious step in the pipeline. The
file is reformatted to a three column format describing TF : TG : TFBS

'''

#import libraries
import os
import csv
import argparse

from cisRegion import readFastaGenerator, genomeLoad

def main (FASTA, out_file):
    '''
    Extract info and reformat to three column
    '''
    # open output file
    with open(out_file, 'w') as OUTPUT:
        # process load the fasta as a dict and process each header/seq pair
        for header, sequence in genomeLoad(FASTA).items():

            # from header retreive TF ensembl id
            TF_ID = header.split('_')[0].split(';')[0]

            # from header retrieve TG ensembl id
            TG_ID = header.split('_')[1].split(';')[0]

            # from header retrieve TG strand
            TG_ST = header.split('_')[1].split(';')[1].split('(')[1].replace(')', '')

            # merge TG id and strand with correct fromatting
            TG_W_ST = TG_ID + '_(' + TG_ST + ')'

            # write to utput
            OUTPUT.write('\t'.join([TF_ID, TG_W_ST, sequence]) + '\n')


# Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('fasta_in', type=str, help='path to .fasta file containing sequences under binding event peaks')
parser.add_argument('out_file', type=str, help='File path of the output')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.fasta_in, args.out_file)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
