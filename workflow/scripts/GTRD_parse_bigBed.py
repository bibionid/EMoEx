#!/usr/bin/python

'''
GTRD_parse_bigBed.py

This script will take a bigbed file from the GTRD database describing the
binding events for one transcription factor in a target species and extract the
sequences under peaks.

WIP
Need to add in the decoding trackDb.txt file from the GTRD schema, the reference
annotation and an uncluttering function.

'''

#import libraries
import csv
import argparse

import pyBigWig as pbw
import numpy    as np
import pandas   as pd

#Define custom functions
def parseBIGBEDinput (bigbed_infile):
    '''
    method 1
    '''
    bb = pbw.open(bigbed_infile)
    print(bb.chroms())
    print(bb.header())

    for chrom, limit in bb.chroms().items():
        print(chrom, 0, limit)
        for bindingEvent in bb.entries(chrom, 0, limit):
            start = bindingEvent[0]
            stop  = bindingEvent[1]
            name  = bindingEvent[2]
            print(chrom, start, stop, name)



def main (inBIGBED, out_path):
    '''
    main
    '''
    parseBIGBEDinput(inBIGBED)

#Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('BIGBED_in',  type=str, help='path to GTRD binding events in BIGBED format')
parser.add_argument('out_dir', type=str, help='path to directory where the figure should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.BIGBED_in, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
