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


# def ensemblIDextract (ensemblGeneBed):
#     '''
#     use reference annotation to link ensembl ids to gene ids
#     '''
#     #blank dict for storage
#     ensembl_ids = {}
#
#     #store ensembl ids in memory
#     with open(ensemblGeneBed) as reference_bed:
#         for line in csv.reader(reference_bed, delimiter = '\t'):
#             if line[6] not in ensembl_ids:
#                 ensembl_ids[line[6]] = line[3]
#             else:
#                 print('duplicate gene id')
#                 print(line)
#                 sys.exit()
#     return(ensembl_ids)

def main (UPSTREAM, DOWNSTREAM, WINDOW, out_path):
    '''
    main function - to edit
    '''
    # create blank storage dictionary
    bed_closest = {}

    # load annotation of closest genes
    bed2Dict(UPSTREAM, bed_closest)
    bed2Dict(DOWNSTREAM, bed_closest)

    for k, v in bed_closest.items():
        if len(v) > 2:
            print(k, '\tskipping due to duplication')
        elif v[0][0] == '.' and v[1][0] == '.':
            print(k, '\tskipping due to no closest')
        else:
            if v[0][5] == v[1][5]:
                if v[0][5] == '+':
                    if int(v[0][7]) <= 0:
                        if int(v[1][7]) > 2*int(WINDOW):
                            print(k,'no closest within', str(2*int(WINDOW)))
                        else:
                            print(k, '-upstream-', v[1])
                    else:
                        print(k, v, 'samesies')
                else:
                    if int(v[0][7]) >= 0:
                        if int(v[0][7]) < 0-(2*int(WINDOW)):
                            print(k,'no closest within', str(2*int(WINDOW)))
                        else:
                            print(k, v[0], '-downstream-')
                    else:
                        if abs(int(v[0][7])) < int(int(v[1][7])):
                            print(k, v[0], '-downstream-')
                        else:
                            print(k, '-upstream-', v[1])

            else:
                print(k, v, 'differ')

    # load schema dictionary
    # TF_ID = trackDBparser(trackDB)

    # extract gene name using this schema
    # gene_name = TF_ID[os.path.basename(BIGBED).split('_')[1]][0]

    #use reference annotation to retreive ensembl id
    # ensembl_id = ensemblDecoder[gene_name]

    # open bigbed file
    # bb = pbw.open(BIGBED)

    # intergrate this information with coordinates of binding events, write .bed
    # with open(os.path.join(out_path, 'GTRD_BED', os.path.basename(BIGBED).replace('.bb', '.bed')), 'w') as outfile:
    #     for chrom, limit in bb.chroms().items():
    #
    #         for bindingEvent in bb.entries(chrom, 0, limit):
    #             start = bindingEvent[0]
    #             stop  = bindingEvent[1]
    #
    #             # correction to chrom id made below with .replace() - this may need alteration for generalisation
    #             outfile.write('\t'.join([str(x) for x in [chrom.replace('chr', ''), start, stop, ensembl_id, gene_name]]) + '\n')

# Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('upstream_in',   type=str, help='path to .bed file describing closest upstream genes to each binding event')
parser.add_argument('downstream_in', type=str, help='path to .bed file describing closest downstream genes to each binding event')
parser.add_argument('cis_window',    type=int, help='integer value of nt window for putative cis interaction')
parser.add_argument('out_dir',       type=str, help='path to directory where the figure should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.upstream_in, args.downstream_in, args.cis_window, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
