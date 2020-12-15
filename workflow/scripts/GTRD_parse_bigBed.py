#!/usr/bin/python

'''

GTRD_parse_bigBed.py

This script takes a bigbed file from the GTRD database describing the
binding events for one transcription factor in a target species and transforms
to bed format. 

It also uses the custom reference gene annotation bed used elsewhere in EMotEP
to retreive the Ensembl id of the transcription factor in the
bigbed file.

'''

#import libraries
import os
import csv
# Legacy as requires internet connection
#import mygene
import argparse

import pyBigWig as pbw

#Define custom functions
def trackDBparser (trackDB):
    '''
    Retreive TF gene names fro GTRD Database shcema
    '''
    # blank storage dictionary
    tf_indentifiers = {}

    # iterate through the GTRD schema, extract tf indentifiers, store in memory
    with open(trackDB) as tDB:
        for line in tDB:
            # filter short lines
            if len(line.split()) > 1:
                # extract tf identifiers
                if line.split()[1] == 'tf_name':
                    # extract coordinates of binding events
                    for item in line.split()[3:]:
                        tf_ident  = item.split('=')[0]
                        gene_name = item.split('=')[1].split('\\')[-1]
                        if tf_ident not in tf_indentifiers:
                            tf_indentifiers[tf_ident] = [gene_name]
                        else:
                            if gene_name not in tf_indentifiers[tf_ident]:
                                tf_indentifiers[tf_ident].append(gene_name)
    return(tf_indentifiers)

def ensemblIDextract (ensemblGeneBed):
    '''
    use reference annotation to link ensembl ids to gene ids
    '''
    #blank dict for storage
    ensembl_ids = {}

    #store ensembl ids in memory
    with open(ensemblGeneBed) as reference_bed:
        for line in csv.reader(reference_bed, delimiter = '\t'):
            if line[6] not in ensembl_ids:
                ensembl_ids[line[6]] = line[3]
            else:
                print('duplicate gene id')
                print(line)
                sys.exit()
    return(ensembl_ids)

def main (BIGBED, trackDB, BED, out_path):
    '''
    Create a output dir if needed, generate decoding library, exctract gene name
    , use this to extract Ensembl ID. Write to bed file with gene name and
    Ensembl id.
    '''
    # create outdir if needed
    if not os.path.isdir(os.path.join(out_path, 'GTRD_BED')):
        os.mkdir(os.path.join(out_path, 'GTRD_BED'))

    # Legacy as requires internet connection
    # load mygene tool
    #mg = mygene.MyGeneInfo()

    # load annotation to retreive ensemble ids
    ensemblDecoder = ensemblIDextract(BED)

    # load schema dictionary
    TF_ID = trackDBparser(trackDB)

    # extract gene name using this schema
    gene_name = TF_ID[os.path.basename(BIGBED).split('_')[1]][0]

    # Legacy as requires internet connection
    # use mygene tool to extract ensembl id
    #ensembl_id = mg.query(gene_name, fields = 'ensembl.gene', species = 'fruitfly')['hits'][0]['ensembl']['gene']

    #use reference annotation to retreive ensembl id
    ensembl_id = ensemblDecoder[gene_name]

    # open bigbed file
    bb = pbw.open(BIGBED)

    # intergrate this information with coordinates of binding events, write .bed
    with open(os.path.join(out_path, 'GTRD_BED', os.path.basename(BIGBED).replace('.bb', '.bed')), 'w') as outfile:
        for chrom, limit in bb.chroms().items():

            for bindingEvent in bb.entries(chrom, 0, limit):
                start = bindingEvent[0]
                stop  = bindingEvent[1]

                outfile.write('\t'.join([str(x) for x in [chrom, start, stop, ensembl_id, gene_name]]) + '\n')

# Define arguments used in the script
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('BIGBED_in',  type=str, help='path to GTRD binding events in BIGBED format')
parser.add_argument('trackDB_in', type=str, help='path to GTRD textDB file for decoding TF identity')
parser.add_argument('BED_in',     type=str, help='path to reference gene annotation in custom BED format')
parser.add_argument('out_dir',    type=str, help='path to directory where the figure should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.BIGBED_in, args.trackDB_in, args.BED_in, args.out_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
