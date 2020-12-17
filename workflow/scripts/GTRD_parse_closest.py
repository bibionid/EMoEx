#!/usr/bin/python

'''

GTRD_parse_closest.py

This script takes 2 .bed files created by ealier steps in the pipeline. These
descibe the genes closest to each binding event in both up- and down stream.
This script parses them together, retaining the association between each binding
event and the most logical closest gene.

This script returns a modified .bed file containing binding event coordinates
and the id of the gene association.

'''

#import libraries
import os
import csv
import argparse

#Define custom functions
def bed2Dict (BED_in, storage_dict):
    '''
    Store modified .bed closest in memory
    '''
    # iterate through the input .bed file
    with open(BED_in) as closest_annotations:
        for line in csv.reader(closest_annotations, delimiter = '\t'):
            # use the coords of the binding event and the TF id as an identifier
            if '_'.join(line[:5]) not in storage_dict:
                # store the closest gene coords and distance info
                storage_dict['_'.join(line[:5])] = [line[5:]]
            else:
                storage_dict['_'.join(line[:5])].append(line[5:])


def createOutputLine(in_key, in_vals, output):
    '''
    Created the modified .bed output line style for this script
    '''
    #split up the bindivent coords from the in_key
    be_coords = in_key.split('_')

    #create an identifier: TF ensemble - TF gene name _ Target Gene ensembl - Target Gene name
    ident = '_'.join([';'.join(be_coords[3:]), ';'.join([in_vals[3], in_vals[6]])])

    #create the desired format
    #event coords, gene strand, identifier
    output.write('\t'.join(be_coords[:3]) + '\t' + ident + '\t.\t' + in_vals[5] + '\n')


def main (UPSTREAM, DOWNSTREAM, WINDOW, out_file):
    '''
    Collate the flanking genes by binding event coordinate, then print .bed file
    output based on the strand and proximity of flanking genes. Binding events
    with no flanking genes are discarded.
    '''
    # create blank storage dictionary
    bed_closest = {}

    # load annotation of closest genes
    bed2Dict(UPSTREAM, bed_closest)
    bed2Dict(DOWNSTREAM, bed_closest)

    # open the output
    with open(out_file, 'w') as OUTPUT:
        # iterate through the binding event dictionary created above
        for binding_event_coords, closest_genic_annotations in bed_closest.items():
            # skip events that have >1 closest gene in either direction
            if len(closest_genic_annotations) > 2:
                continue

            # skip events that have no genes within 10kb in either direction
            elif closest_genic_annotations[0][0] == '.' and closest_genic_annotations[1][0] == '.':
                continue

            # filter the remaining events
            else:

                # if the flanking genes are on the same strand
                if closest_genic_annotations[0][5] == closest_genic_annotations[1][5]:

                    # if the upstream gene is on the forward strand
                    if closest_genic_annotations[0][5] == '+':

                        # ensure gene downstream event is within range
                        if int(closest_genic_annotations[1][7]) > 2*int(WINDOW):
                            continue
                        else:
                            # print the .bed line
                            createOutputLine(binding_event_coords, closest_genic_annotations[1], OUTPUT)

                    # if upstream gene is on the reverse strand
                    else:

                        # ensure gene upstream event is within range
                        if int(closest_genic_annotations[0][7]) < 0-(2*int(WINDOW)):
                            continue
                        else:
                            # print the .bed line
                            createOutputLine(binding_event_coords, closest_genic_annotations[0], OUTPUT)

                # if flanking genes are on opposite strands
                else:

                    # if there is no gene upstream
                    if closest_genic_annotations[0][6] == '.':

                        # if the downstream gene is forward strand
                        if closest_genic_annotations[1][5] == '+':
                            # print the .bed line
                            createOutputLine(binding_event_coords, closest_genic_annotations[1], OUTPUT)
                        else:
                            continue

                    # if there is no gene downstream
                    elif closest_genic_annotations[1][6] == '.':

                        #if the upstream gene is on the reverse strand
                        if closest_genic_annotations[0][5] == '-':
                            # print the .bed line
                            createOutputLine(binding_event_coords, closest_genic_annotations[0], OUTPUT)
                        else:
                            continue

                    #if there are genes both up- and downstream
                    else:

                        # if the binding event is closer to the upstream gene
                        if abs(int(closest_genic_annotations[0][7])) < abs(int(closest_genic_annotations[1][7])):

                            # if the upstream gene is reverse strand
                            if closest_genic_annotations[0][5] == '-':
                                # print the .bed line
                                createOutputLine(binding_event_coords, closest_genic_annotations[0], OUTPUT)
                            else:
                                # print the .bed line for the downstream gene
                                createOutputLine(binding_event_coords, closest_genic_annotations[1], OUTPUT)

                        # if the binding event is closer to the downstream gene
                        else:

                            #if the downstream gene is on the forward strand
                            if closest_genic_annotations[1][5] == '+':
                                # print the .bed line
                                createOutputLine(binding_event_coords, closest_genic_annotations[1], OUTPUT)
                            else:
                                # # print the .bed line for the upstream gene
                                createOutputLine(binding_event_coords, closest_genic_annotations[1], OUTPUT)

# Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('upstream_in',   type=str, help='path to .bed file describing closest upstream genes to each binding event')
parser.add_argument('downstream_in', type=str, help='path to .bed file describing closest downstream genes to each binding event')
parser.add_argument('cis_window',    type=int, help='integer value of nt window for putative cis interaction')
parser.add_argument('out_file',      type=str, help='File path of the output')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.upstream_in, args.downstream_in, args.cis_window, args.out_file)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
