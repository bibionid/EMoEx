#!/usr/bin/python

'''
cisRegion_displot.py

This script will make a distribution plot from the bedfile genereated by cisRegion.py

'''

#import libraries
import csv
import argparse

import numpy   as np
import pandas  as pd
import seaborn as sns

#Define custom functions
def parseBEDinput (bed_infile):
    '''
    Read annotations from bed file into dictionary

    Load dictionary as pandas DataFrame

    Return both
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
    Uses the seaborn library to plot the distribution of the cisRegion lengths
    '''
    #parse the annotations in BED format
    BedAnnotData = parseBEDinput(inBED)

    #extract the dataframe
    BedAnnotData_df = BedAnnotData[1]

    #calculate length of promoter annotations
    promoterLengths = BedAnnotData_df['stop'] - BedAnnotData_df['start']

    #plot the histogram
    dist = sns.displot(data = promoterLengths, x = promoterLengths, bins = 50)
    dist.set(xlabel = "Length (nucleotides)", ylabel = "Count")
    dist.fig.set_size_inches(16,6)

    #write figure
    dist.savefig(out_path + '/cisRegions_lengthDist.pdf')

    print('\nAnnotation size distribution plot written to: ' + out_path + '/cisRegions_lengthDist.pdf\n')

#Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('BED_in',  type=str, help='path to promoter annotation in BED format')
parser.add_argument('out_dir', type=str, help='path to directory where the figure should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.BED_in, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
