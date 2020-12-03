#!/usr/bin/python
'''
This script will make a distribution plot from the bedfile genereated by promBED_fromGeneBED_20kbHardCoded.py
'''
import os
import sys
import csv
import time
import argparse
# import statistics

import numpy as np
import pandas  as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parseBEDinput (bed_infile):
    '''
    read annotations from bed file into dictionary

    load dictionary as pandas DataFrame

    return both
    '''
    #create blank dictionary
    bed_dict = {}

    #create blank dataframe
    bed_dataframe = pd.DataFrame(columns = ['gene_name', 'start', 'stop', 'strand'])

    # open file and stroe in blank dict above
    with open(bed_infile) as input_file:
        for line in csv.reader(input_file, delimiter = '\t'):

            #extract gene name
            gene_name = line[3]

            #check if gene is duplicated in the dict
            if gene_name not in bed_dict:
                #if it is not, store the bed line
                bed_dict[gene_name] = line

                bed_dataframe = bed_dataframe.append(pd.Series([gene_name, int(line[1]), int(line[2]), line[5]], index = ['gene_name', 'start', 'stop', 'strand']), ignore_index = True)
            else:
                #if there is a duplication in the infile throw an error
                raise Exception('\n\tAnnotation dictionary construction stopped\n\tGene {} duplicated in the in file'.format(gene_name))

    return((bed_dict, bed_dataframe))


def main (inBED, out_path):
    '''
    poopy
    '''
    #parse the annotations in BED format
    BedAnnotData = parseBEDinput(inBED)

    #extract the dataframe
    BedAnnotData_df = BedAnnotData[1]

    #set plot layout
    plt.figure(figsize = (16, 6))
    plt.tight_layout()

    #calculate length of promoter annotations
    promoterLengths = BedAnnotData_df['stop'] - BedAnnotData_df['start']

    #make the lengths a list
    promoterLengthList = promoterLengths.values.tolist()

    #plot the histogram
    dist = sns.distplot(promoterLengthList, bins = 50, kde = False)
    dist.set(xlabel = "Length (nucleotides)", ylabel = "Count")

    #write figure
    figure   = dist.get_figure()

    figure.savefig(out_path + '/Om_5kbRegions_lengthDist.pdf')

    print('\nAnnotation size distribution plot written to: ' + out_path + '/Om_5kbRegions_lengthDist.pdf\n')

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('BED_in', type=str, help='path to promoter annotation in BED format')
parser.add_argument('out_dir', type=str, help='path to directory where the figure should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.BED_in, args.out_dir)
