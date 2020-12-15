#!/usr/bin/python

'''
cisRegion.py

Extract a user defined region upstream of input annotations

Developed on protein coding genes in vertebrate genomes

Accounts for overlapping annotations and those that share start sites

Input: annotation (.bed), genome (.fasta)
Output cisRegion annotation (.bed), list of flagged genes sharing start & strand (.txt)
'''

#import libraries
import os
import sys
import csv
import argparse


def asint(s):
    '''
    returns an integer value
    '''
    try: return int(s), ''
    except ValueError: return sys.maxsize, s

def readFastaGenerator(seqin):
    '''
    creates a generator the reads a fasta file sequence by sequence
    '''
    for line in seqin:
        yield line.strip()[1:], next(seqin).strip()

def genomeLoad(genome_fasta):
    '''
    loads a genome fasta into a dictionary
    '''
    seq_dict = {}
    with open(genome_fasta) as seq_file:
        seq_gen = readFastaGenerator(seq_file)
        for _id, _seq in seq_gen:
            seq_dict[_id] = _seq
    return(seq_dict)

def chromLimitFind (chrom_set, genomeFastaDict):
    '''
    Compiles a dictionary containing the length of each scaffold
    '''
    chrLim_dict = {}
    for chrom in chrom_set:
        try:
            chrLim_dict[chrom] = len(genomeFastaDict[chrom])
        except KeyError :
            print ('Scaffold', chrom, 'not in fasta - omitting')
            pass

    return(chrLim_dict)


def scaffoldLister(BED_IN):
    '''
    Creates a set of scaffold ids - this is redundant to the above (change in
    ver2.0)
    '''
    chr_list = []
    with open(BED_IN) as GENEBED_archive:
        for ge in GENEBED_archive:
            if str(ge)[0] != '#':
                chr_list.append(ge.split('\t')[0])
    chr_set = set(chr_list)
    return(chr_set)

def outputFormat(output, scaffold, start, stop, gene_id, gene_strand):
    '''
    write coordinates to output in bed format
    '''
    print('Keeping: ' + '\t'.join([str(scaffold), str(start),str(stop), str(gene_id),'1', str(gene_strand)]))
    output.write('\t'.join([str(scaffold), str(start),str(stop), str(gene_id),'1', str(gene_strand) + '\n']))

def hugeCisRegionCallingBehemoth(SCAF_LIST, BED_IN, OUTPUT_DIR, SCAFF_LIMS, CIS_WINDOW):
    '''
    This highly repetitive function works to call cis regions upstream of
    annotated genes. It excludes cis regions that overlap coding sequence. It
    will also split shared cisRegions 50/50 between genes on opposite strands.
    It also excluded genes on the same strand which share start site. These
    putative annotation errors are reported in a log file.

    This needs to be heavily unpicked and edited but is currently functional,
    so maybe in ver2.0...

    This script is so complex as it was developed on experimental annotations of
    short read mammalian genomes. These annotations contained many questionable
    gene models which overlapped each other and sometimes could be multiply
    nested inside each other. This script should be able to handle such
    scenarios, but also still contains a kill switch if the annotation is too
    messy.
    '''
    remove_list       = []
    same_strand_start = []
    filename = os.path.basename(BED_IN)

    for chrom in sorted(SCAF_LIST):
        print('\n\nWorking with chromosome ' + chrom + '\n')
        with open(BED_IN) as GENEBED_archive:
            gene_dict = {}
            for gene in csv.reader(GENEBED_archive, delimiter = '\t'):
                if gene[0][0] != '#':
                    if gene[0] == chrom:
                        try:
                            if gene[3] in gene_dict:
                                print('\n\t!!!\tBreak during gene library construction due to duplicate genes in annotation')
                                sys.exit()
                            else:
                                gene_dict[gene[3]] = [gene[0], gene[1], gene[2], gene[5], gene[3]]

                        except:
                            print('\n\n--- Error during gene library construction ---\n\n')
                            print('\n'+ str(gene) +'\n')
                            sys.exit()

            sorted_gene_list = [(k, gene_dict[k]) for k in sorted(gene_dict, key=asint)]

            chrom_dict = {}
            for gen in sorted_gene_list:
                if gen[1][1] not in chrom_dict:
                    chrom_dict[gen[1][1]] = [gen]

                else:
                    if chrom_dict[gen[1][1]][0][1][3] == gen[1][3]:

                        print('Recording gene ' + gen[1][4] + ' with promoter overlap the same strand')
                        same_strand_start.append(chrom_dict[gen[1][1]][0][1][4])
                        same_strand_start.append(gen[1][4])

                        chrom_dict[gen[1][1]].append(gen)

                    else:
                        chrom_dict[gen[1][1]].append(gen)

            sorted_chrm_list = [(k, chrom_dict[k]) for k in sorted(chrom_dict, key=asint)]

            count           = 0
            switch          = 0
            stored_start    = 0
            stored_stop     = 0
            stored_chr      = '@'
            stored_strand   = '+'
            strand_skip     = '@'
            # with open(os.path.join(OUTPUT_DIR, filename[:-3] + str(CIS_WINDOW) + 'nt_cisRegions.stranded.bed'), 'w') as out_file:
            with open(os.path.join(OUTPUT_DIR, 'cisRegions.bed'), 'a') as out_file:
                for gene in sorted_chrm_list:

                    #this removes genes that start at 0 on a scaffold
                    if int(gene[0]) == 0 and gene[1][0][1][3] == '+':
                        remove_list.append(gene[1][0][1][4])
                        print('\n\n\nRemoving gene with no possible promoter\n\n\n\n')
                        stored_start = 1
                        count        = count + 1
                        continue

                    if switch == 0:
                        stored_stop = gene[1][0][1][2]
                        switch      = 1

                    if len(gene[1]) > 1:
                        print ('\n\n\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' + str(len(gene[1])) + ' DUPLICATES!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n')
                        count = count + 1
                        stored_dup_start = '@'
                        stored_dup_stop  = '@'
                        stored_dup_prom  = '@'
                        for d_c, dup in enumerate(gene[1]):
                            print('--- On gene =\t' + dup[0] + ' ' + dup[1][3])
                            print('pum')
                            if d_c > 0:

                                if dup[1][4] in same_strand_start:
                                    print('\t' + str(dup[1][4]) + ' is a start site same strand duplicate - mirroring promoter\n')
                                    if stored_dup_prom == '@':
                                        print('\nNo stored promoter to mirror - annotation not recorded\n')
                                        remove_list.append(dup[1][4])
                                    else:
                                        dup[1].append(stored_dup_prom)

                                    if dup[1][3] == '+':
                                        if int(dup[1][2]) > int(stored_stop) and int(dup[1][1]) < int(stored_start):
                                            print('\tRemoving duplicate startsite gene 2 as nested within previous gene\n')
                                            remove_list.append(dup[1][4])
                                        elif int(dup[1][2]) > int(stored_stop):
                                            stored_stop = int(dup[1][2])

                                    elif dup[1][3] == '-':
                                        if int(dup[1][2]) < int(stored_start):
                                            print('Gene overlapping - same start site on the negative strand - annotation not recorded\n')
                                            remove_list.append(dup[1][4])

                                        elif int(dup[1][2]) > int(stored_start):
                                            print('Previous gene overlapping - same start site on the negative strand - annotation not recorded\n')
                                            remove_list.append(stored_dup_name)
                                            remove_list.append(dup[1][4])

                                            try:
                                                dup[1][5] = int(dup[1][2]) + CIS_WINDOW
                                            except(IndexError):
                                                dup[1].append(int(dup[1][2]) + CIS_WINDOW)
                                            stored_start = int(dup[1][2])

                                        else:
                                            print(str(dup[1][4]) + ' is a start site same strand duplicate - mirroring promoter\n')
                                            dup[1].append(stored_dup_prom)


                                    if d_c == len(gene[1]) - 1:
                                        stored_strand = dup[1][3]
                                else:
                                    if dup[1][3] == '-':

                                        if int(dup[1][2]) >= int(stored_dup_stop):
                                            dup[1].append(int(dup[1][2]) + CIS_WINDOW)
                                            print('Gene ' + dup[1][4] + ',\tstop:\t' + str(dup[1][2]) + ',\treplaces with ' + dup[1][2] + ' + ' + str(CIS_WINDOW) + '\n')
                                            stored_start = int(dup[1][2])
                                            stroed_stop  = int(dup[1][1])
                                        else:
                                            stored_start = int(stored_dup_stop)
                                            stroed_stop  = int(dup[1][1])
                                            remove_list.append(dup[1][4])
                                            print('Promoter overlapping CDS - promoter not annotated')
                                    elif dup[1][3] == '+':
                                        if stored_strand == '-':
                                            dist = int((int(dup[1][1]) - int(stored_start)) / 2)

                                            if dist <= 1:
                                                print('\nOverlapping gene - promoter not annotated')
                                                remove_list.append(dup[1][4])
                                                print('\nGene prior to overlapping pair also promoter not annotated due to overlap')
                                                remove_list.append(sorted_chrm_list[count - 2][1][0][1][4])
                                                continue

                                            elif dist < CIS_WINDOW :
                                                print('Duplicate start site gene - calculating mid distance for promoter')
                                                print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - ' + str(int((int(dup[1][1])-int(stored_start))/2)))
                                                dup[1].append(int(int(dup[1][1]) - int((int(dup[1][1])-int(stored_start))/2)))

                                                print('Adjusting previous promoter stop to: ' + str(stored_start) + ' + ' +  str(int((int(dup[1][1])-int(stored_start))/2)))
                                                sorted_chrm_list[count-2][1][0][1][5] = int(int(stored_start) + int((int(dup[1][1])-int(stored_start))/2))


                                            else:
                                                print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                                dup[1].append(int(dup[1][1]) - CIS_WINDOW)

                                        else:
                                            if int(dup[1][1]) - CIS_WINDOW >= stored_stop:
                                                dup[1].append(int(dup[1][1]) - CIS_WINDOW)
                                                print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                                stored_start = int(dup[1][1])
                                                stroed_stop  = int(dup[1][2])
                                            else:
                                                dup[1].append(int(stored_stop))
                                                print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + str(stored_stop) + '\n')
                                                stored_start = int(dup[1][1])
                                                stroed_stop  = int(dup[1][2])

                                        if stored_strand == dup[1][3] and int(dup[1][1]) >= int(stored_start):
                                            print('Removing promoter annotation as overlapping previous gene.')
                                            remove_list.append(dup[1][4])

                                        if int(dup[1][2]) > int(stored_dup_start):
                                            print('Promoter of gen prior to ' + str(dup[1][4]) + ' overlapping CDS. Promoter not annotated.\n\n\n')
                                            remove_list.append(stored_dup_name)

                                        if stored_dup_strand == '-' and int(stored_dup_start) > int(stored_start):
                                            stored_start = stored_dup_start
                                            if int(stored_dup_stop) > int(stored_stop):
                                                stored_stop = stored_dup_stop

                            else:
                                stored_dup_name   = dup[1][4]
                                stored_dup_strand = dup[1][3]

                                if dup[1][0] == stored_chr:

                                    if dup[1][3] == stored_strand:

                                        if dup[1][3] == '-':
                                            if int(dup[1][1]) >= int(stored_start) + CIS_WINDOW:
                                                print('\nAdding preliminary promoter annotation to dup gene 1\n')
                                                dup[1].append(int(dup[1][2]) + CIS_WINDOW)
                                                stored_dup_start = int(dup[1][2])
                                                stored_dup_stop  = int(dup[1][1])
                                            else:
                                                print('\nGenes closer than ' + str(CIS_WINDOW) + ' - editting previous promoter \n') ####CHONDO

                                                print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4].strip() + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\treplaces with ' +  str(dup[1][1]))

                                                try:
                                                    if int(dup[1][1]) >= int(sorted_chrm_list[count-2][1][0][1][1]) and int(dup[1][2]) <= int(sorted_chrm_list[count-2][1][0][1][5]):
                                                        print('\nPrevious gene as duplicates nested within it - removing annotations')
                                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4].strip())
                                                        for dup in gene[1]:
                                                            print(dup[0])
                                                            remove_list.append(dup[0])
                                                except IndexError:
                                                    if int(dup[1][1]) >= int(sorted_chrm_list[count-2][1][0][1][1]):
                                                        print('\nPrevious gene as duplicates overlapping it - removing annotations')
                                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4].strip())
                                                        for dup in gene[1]:
                                                            print(dup[0])
                                                            remove_list.append(dup[0])

                                                for loc in range(len(sorted_chrm_list[count-2][1])):
                                                    try:
                                                        sorted_chrm_list[count-2][1][loc][1][5] = int(dup[1][1])
                                                    except IndexError:
                                                        pass

                                                print('\nAdding preliminary promoter annotation to dup gene 1\n')
                                                dup[1].append(int(dup[1][2]) + CIS_WINDOW)
                                                stored_dup_start = int(dup[1][2])
                                                stored_dup_stop  = int(dup[1][1])
                                                stored_dup_prom  = int(dup[1][2]) + CIS_WINDOW
                                                print('\n')

                                        elif dup[1][3] == '+' and int(dup[1][1]) - CIS_WINDOW <= 0:
                                            dup[1].append(0)
                                            print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with 0')
                                            print('\n')
                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                            stored_dup_prom  = int(dup[1][5])

                                        else:
                                            if int(dup[1][1]) >= int(stored_stop) + CIS_WINDOW:
                                                print(stored_strand)
                                                print(stored_start)
                                                print(stored_stop)
                                                print(stored_start - 5000)
                                                print(stored_stop - 5000)

                                                print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                                dup[1].append(int(dup[1][1]) - CIS_WINDOW)
                                                print('\n')
                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])
                                                stored_dup_prom  = int(dup[1][5])

                                            else:
                                                if int(dup[1][1]) > int(stored_start) and int(dup[1][2]) < int(stored_stop):
                                                    print('Duplicate start site nested within previous gene - promoter not annotated')
                                                    remove_list.append(dup[1][4])
                                                    stored_dup_start = int(dup[1][1])
                                                    stored_dup_stop  = int(dup[1][2])

                                                elif int(dup[1][1]) < int(stored_start):
                                                    print('Duplicate start site overlaps previous gene on same strand - promoter not annotated')
                                                    remove_list.append(dup[1][4])
                                                    stored_dup_start = int(dup[1][1])
                                                    stored_dup_stop  = int(dup[1][2])

                                                else:
                                                    if int(dup[1][1]) - int(stored_stop) > 0:
                                                        print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - ' + str(int(dup[1][1]) - int(stored_stop)))
                                                        dup[1].append(int(stored_stop))
                                                        print('\n')
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])
                                                        stored_dup_prom  = int(dup[1][5])

                                                    else:
                                                        print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\toverlaps previous stop: ' + str(stored_stop) + '. No promoter annotated.\n')
                                                        remove_list.append(dup[1][4])
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])

                                    elif dup[1][3] != stored_strand:
                                        print('pum')
                                        if dup[1][3] == '+':
                                            if int(dup[1][1]) < int(stored_start) and int(dup[1][2]) > int(stored_stop):
                                                print('\tRemoving duplicate start site gene 1 as overlapping/nested within previous gene\n')
                                                remove_list.append(dup[1][4])
                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])

                                                print('\tAlso removing previous gene on -ve strand due to overlap\n\n')
                                                remove_list.append(sorted_chrm_list[count-2][1][0][1][4])

                                            else:
                                                if (stored_start + CIS_WINDOW) > int(dup[1][1]) and stored_start < int(dup[1][2]):
                                                    dist = int(dup[1][1]) - stored_start
                                                    print('Stored start:\t' + str(stored_start))
                                                    print('Querey start:\t' + str(dup[1][1]))
                                                    print('The distance between the genes is:\t' + str(dist))

                                                    if dist <= 1:
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])

                                                        print('\nOverlapping gene - promoter not annotated')
                                                        remove_list.append(dup[1][4])
                                                        print('\nGene prior to overlapping pair also promoter not annotated due to overlap')
                                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                                        continue

                                                    else:
                                                        print('Allowed promoter region is\t' + str(int(dist/2)))
                                                        try:
                                                            print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4]) + ',\tstop:\t' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with ' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                            sorted_chrm_list[count-2][1][0][1][5] = int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2)
                                                        except IndexError:
                                                            print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4]) + ',\tstop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\treplaces with ' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                            sorted_chrm_list[count-2][1][0][1].append(int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2))

                                                        print('Gene ' + dup[1][4].strip() + ',\tstart:\t' + str(dup[1][1]) + '\treplaces with ' + str(dup[1][1]) + ' - ' + str(int(dist/2)))
                                                        dup[1].append(int(dup[1][1]) - int(dist/2))
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])
                                                        stored_dup_prom  = int(dup[1][5])
                                                        print('\n')

                                                elif int(dup[1][1]) - CIS_WINDOW <= 0:
                                                    dup[1].append(0)
                                                    print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with 0')
                                                    stored_dup_start = int(dup[1][1])
                                                    stored_dup_stop  = int(dup[1][2])
                                                    stored_dup_prom  = int(dup[1][5])
                                                    print('\n')

                                                else:
                                                    if int(dup[1][1]) - CIS_WINDOW >= int(stored_stop):
                                                        print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                                        gene[1][0][1].append(int(dup[1][1]) - CIS_WINDOW)
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])
                                                        stored_dup_prom  = int(dup[1][5])

                                                    else:
                                                        print('\nGenes closer than ' + str(CIS_WINDOW) + ' - taking shorter promoter\n')
                                                        print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + str(stored_stop))
                                                        dup[1].append(int(stored_stop))
                                                        stored_dup_start = int(dup[1][1])
                                                        stored_dup_stop  = int(dup[1][2])
                                                        stored_dup_prom  = int(dup[1][5])

                                        else:
                                            print('pum')
                                            print(stored_start)
                                            print(stored_stop)
                                            print(dup[1])
                                            dup[1].append(int(dup[1][2]) + CIS_WINDOW)

                                            if int(dup[1][1]) >= int(stored_stop) and int(dup[1][2]) <= int(stored_start):
                                                print('\nDuplicate is nested within earlier gene - removing annotation\n')
                                                remove_list.append(dup[0])
                                            else:
                                                print('Gene ' + dup[1][4] + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][2] + ' + ' + str(CIS_WINDOW) + '\n')
                                                print('\n')
                                            stored_dup_start = int(dup[1][2])
                                            stored_dup_stop  = int(dup[1][1])
                                            stored_dup_prom  = int(dup[1][5])


                                elif dup[1][0] != stored_chr:
                                    print('First of two start site duplicates is at the start of a new scaffold...\n')
                                    if dup[1][3] == '+':
                                        if int(dup[1][1]) - CIS_WINDOW <= 0:
                                            dup[1].append(0)
                                            print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with 0')
                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                            stored_dup_prom  = int(dup[1][5])

                                            print('\n')
                                        else:
                                            dup[1].append(int(dup[1][1]) - CIS_WINDOW)
                                            print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with ' + str(dup[1][1]) + ' - ' + str(CIS_WINDOW) + '\n')
                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                            stored_dup_prom  = int(dup[1][5])

                    else:
                        print('--- On gene =\t' + gene[1][0][0] + ' ' + gene[1][0][1][3])

                        count = count + 1

                        if gene[1][0][1][0] == stored_chr:

                            if gene[1][0][1][3] == stored_strand:

                                if gene[1][0][1][3] == '-':

                                    if int(gene[1][0][1][2]) == int(stored_start) and int(gene[1][0][1][1]) == int(stored_stop):

                                        print('\n' + str(gene))
                                        print(stored_start)
                                        print('\n\n\nExact duplicate gene discovered on opposite strand - manually consult annotation...\n\n\n')
                                        sys.exit()

                                    print('Gene ' + gene[1][0][1][4].strip() + ' stop:\t' + str(gene[1][0][1][2]) + '\treplaces with ' + gene[1][0][1][2] + ' + ' + str(CIS_WINDOW) + '\n')
                                    gene[1][0][1].append(int(gene[1][0][1][2]) + CIS_WINDOW)

                                    strand_skip = '@'

                                    if int(gene[1][0][1][1]) >= int(stored_start) + CIS_WINDOW:

                                        stored_start = int(gene[1][0][1][2])
                                        stored_stop  = int(gene[1][0][1][1])

                                    else:

                                        print('\nGenes closer than ' + str(CIS_WINDOW) + ' - editting previous promoter\n')

                                        if len(sorted_chrm_list[count-2][1]) == 1:

                                            if int(stored_start) > int(gene[1][0][1][2]) > int(stored_stop):
                                                print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                                remove_list.append(gene[1][0][1][4])

                                            if int(gene[1][0][1][2]) >= int(stored_start) > int(gene[1][0][1][1]):
                                                print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                                remove_list.append(gene[1][0][1][4])
                                                print('Previous gene nested or overlapping previous gene on same strand - promoter annotation removed\n')
                                                remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                                stored_start = int(gene[1][0][1][2])
                                                stored_stop  = int(gene[1][0][1][1])

                                            else:

                                                if gene[1][0][1][4] not in remove_list:

                                                    if sorted_chrm_list[count-2][1][0][0] not in remove_list:

                                                        print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])

                                                        if len(sorted_chrm_list[count-2][1][0][1]) == 6:
                                                            sorted_chrm_list[count-2][1][0][1][5] = int(gene[1][0][1][1])
                                                        else:
                                                            sorted_chrm_list[count-2][1][0][1].append(gene[1][0][1][1])

                                                        if sorted_chrm_list[count-2][1][0][1][4] not in remove_list:
                                                            stored_start = int(gene[1][0][1][2])
                                                            stored_stop  = int(gene[1][0][1][1])

                                                    else:
                                                        print('Gene to be edited is to be removed - skipping')
                                                        stored_start = int(gene[1][0][1][2])
                                                        stored_stop  = int(gene[1][0][1][1])

                                        elif len(sorted_chrm_list[count-2][1]) == 2:
                                            if sorted_chrm_list[count-2][1][1][1][3] == gene[1][0][1][3]:
                                                print('\tGene ' + sorted_chrm_list[count-2][1][1][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][1][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])

                                                try:
                                                    sorted_chrm_list[count-2][1][1][1][5] = int(gene[1][0][1][1])
                                                except(IndexError):
                                                    print('\n\t!!! Missing annotation appended...\n')
                                                    sorted_chrm_list[count-2][1][1][1].append(int(gene[1][0][1][1]))
                                            else:
                                                print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])

                                                try:
                                                    sorted_chrm_list[count-2][1][0][1][5] = int(gene[1][0][1][1])
                                                except(IndexError):
                                                    print('\n\t!!! Missing annotation appended...\n')
                                                    sorted_chrm_list[count-2][1][0][1].append(int(gene[1][0][1][1]))
                                            stored_start = int(gene[1][0][1][2])
                                            stored_stop  = int(gene[1][0][1][1])

                                    print('\n')

                                elif gene[1][0][1][3] == '+':

                                    if stored_start == 0:
                                        gene[1][0][1].append(0)
                                        print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with 0')
                                        stored_start  = int(gene[1][0][1][1])
                                        stored_stop   = int(gene[1][0][1][2])
                                        stored_strand = gene[1][0][1][3]
                                        print('\n')

                                    else:

                                        if int(gene[1][0][1][1]) - CIS_WINDOW >= int(stored_stop):
                                            print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with ' + gene[1][0][1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                            gene[1][0][1].append(int(gene[1][0][1][1]) - CIS_WINDOW)
                                            stored_start = int(gene[1][0][1][1])
                                            stored_stop  = int(gene[1][0][1][2])
                                            strand_skip = '@'

                                        elif int(stored_stop) > int(gene[1][0][1][1]) > int(stored_start):
                                            print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                            remove_list.append(gene[1][0][1][4])

                                        else:

                                            print('\nGenes closer than ' + str(CIS_WINDOW) + ' - taking shorter promoter\n')

                                            print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with ' + gene[1][0][1][1] + ' - ' + str(int(gene[1][0][1][1]) - int(stored_stop)))

                                            gene[1][0][1].append(int(stored_stop))
                                            stored_start = int(gene[1][0][1][1])
                                            stored_stop  = int(gene[1][0][1][2])
                                            strand_skip = '@'

                                        print('\n')

                                    stored_strand = gene[1][0][1][3]

                            elif gene[1][0][1][3] != stored_strand:
                                print ('Current strand: '+ gene[1][0][1][3] + ' | Stored strand: ' + stored_strand)

                                if gene[1][0][1][3] == '+':

                                    if int(stored_start + CIS_WINDOW) > int(gene[1][0][1][1]) and stored_start < int(gene[1][0][1][2]):
                                        dist = int(gene[1][0][1][1]) - int(stored_start)
                                        print('Stored start:\t' + str(stored_start))
                                        print('Querey start:\t' + str(gene[1][0][1][1]))
                                        print('The distance between the genes is:\t' + str(dist))

                                        if dist <= 1:
                                            print('\nOverlapping gene - promoter not annotated')

                                            remove_list.append(gene[1][0][1][4])
                                            print('Also removing previous promoter annotation due to overlap')
                                            remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                            stored_start = int(gene[1][0][1][1])
                                            stored_stop  = int(gene[1][0][1][2])
                                            print('\n')
                                            strand_skip  = '@'

                                        else:

                                            print('Allowed promoter region is\t' + str(int(dist/2)))

                                            if len(sorted_chrm_list[count-2][1]) == 1:

                                                if strand_skip == '@':
                                                    print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                    print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                                    gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))

                                                    sorted_chrm_list[count-2][1][0][1][5] = int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2)
                                                    stored_start = int(gene[1][0][1][1])
                                                    stored_stop  = int(gene[1][0][1][2])
                                                    print('\n')
                                                else:
                                                    print(strand_skip)
                                                    try:
                                                        print('PREVIOUS Gene ' + str(sorted_chrm_list[count-nesting_counter][1][0][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-nesting_counter][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-nesting_counter][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                        print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                                        gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))

                                                        sorted_chrm_list[count-nesting_counter][1][0][1][5] = int(sorted_chrm_list[count-nesting_counter][1][0][1][2]) + int(dist/2)
                                                        stored_start = int(gene[1][0][1][1])
                                                        stored_stop  = int(gene[1][0][1][2])
                                                        print('\n')
                                                    except IndexError:
                                                        print(str(sorted_chrm_list[count-4][1][0][1][4].strip()))
                                                        print(str(sorted_chrm_list[count-3][1][0][1][4].strip()))
                                                        print(str(sorted_chrm_list[count-2][1][0][1][4].strip()))
                                                        print(count)
                                                        print(nesting_counter)
                                                        print(str(sorted_chrm_list[count-nesting_counter][1][0][1][4].strip()))
                                                        print('fooked')
                                                        # print('PREVIOUS Gene ' + str(sorted_chrm_list[count-4][1][0][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-4][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-3][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                        # print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                                        # gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))
                                                        #
                                                        # sorted_chrm_list[count-4][1][0][1][5] = int(sorted_chrm_list[count-4][1][0][1][2]) + int(dist/2)
                                                        # stored_start = int(gene[1][0][1][1])
                                                        # stored_stop  = int(gene[1][0][1][2])
                                                        # print('\n')

                                            elif len(sorted_chrm_list[count-2][1]) > 1:

                                                if dup[1][3] != sorted_chrm_list[count-2][1][0][1][3]:
                                                    try:
                                                        print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4]) + '  stop:  ' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                    except(IndexError):
                                                        print('\n\t!!! Missing annotation appended to dup 1...\n')
                                                        sorted_chrm_list[count-2][1][0][1].append(int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2))
                                                else:

                                                    try:
                                                        print('Gene ' + str(sorted_chrm_list[count-2][1][-1][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-2][1][-1][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][-1][1][2]) + ' + ' + str(int(dist/2)))

                                                        for loc in range(len(sorted_chrm_list[count-2][1])):
                                                            sorted_chrm_list[count-2][1][loc][1][5] = int(sorted_chrm_list[count-2][1][-1][1][2]) + int(dist/2)

                                                    except(IndexError):
                                                        print('\n\t!!! Missing annotation appended to dup 2...\n')
                                                        sorted_chrm_list[count-2][1][-1][1].append(int(sorted_chrm_list[count-2][1][-1][1][2]) + int(dist/2))

                                                print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                                gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))

                                                try:
                                                    print(sorted_chrm_list[count-2][1][1][1][5])
                                                except(KeyError):
                                                    print('None stored')

                                                stored_start = int(gene[1][0][1][1])
                                                stored_stop  = int(gene[1][0][1][2])
                                                print('\n')
                                            strand_skip  = '@'

                                    elif int(gene[1][0][1][1]) - CIS_WINDOW <= 0:

                                        if int(gene[1][0][1][1]) - int(stored_stop) > 0:
                                            gene[1][0][1].append(0)
                                            print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with 0\n')
                                            stored_start = int(gene[1][0][1][1])
                                            print('\n')

                                        else:
                                            remove_list.append(gene[1][0][1][4])
                                            print('Gene overlapping previous on negative strand - promoter not annotated\n')

                                        strand_skip  = '@'

                                    elif int(stored_start + CIS_WINDOW) > int(gene[1][0][1][1]) and int(stored_start) >= int(gene[1][0][1][2]):
                                        print('Gene nested within another - promoter not annotated')
                                        remove_list.append(gene[1][0][1][4])
                                        if strand_skip == '@':
                                            nesting_counter = 3
                                        else:
                                            nesting_counter = nesting_counter + 1
                                        print(nesting_counter)
                                        strand_skip = stored_strand
                                        print('\n')

                                    else:
                                        gene[1][0][1].append(int(gene[1][0][1][1]) - CIS_WINDOW)
                                        print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with ' + gene[1][0][1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                        stored_start = int(gene[1][0][1][1])
                                        stored_stop  = int(gene[1][0][1][2])
                                        strand_skip  = '@'
                                        print('\n')

                                else:
                                    gene[1][0][1].append(int(gene[1][0][1][2]) + CIS_WINDOW)
                                    print('Gene ' + gene[1][0][1][4].strip() + ' stop:\t' + str(gene[1][0][1][2]) + '\treplaces with ' + gene[1][0][1][2] + ' + ' + str(CIS_WINDOW) + '\n')

                                    if strand_skip == gene[1][0][1][3] and int(stored_start) + CIS_WINDOW > int(gene[1][0][1][1]):

                                        print('\nGene prior to nested genes closer than ' + str(CIS_WINDOW) + ' - editting previous promoter\n')

                                        if len(sorted_chrm_list[count-2][1]) == 1:

                                            if int(stored_start) > int(gene[1][0][1][2]) > int(stored_stop):
                                                print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                                remove_list.append(gene[1][0][1][4])

                                            if int(gene[1][0][1][2]) > int(stored_start) > int(gene[1][0][1][1]):
                                                print('Previous gene nested or overlapping previous gene on same strand - promoter annotation removed\n')
                                                remove_list.append(sorted_chrm_list[count-2][1][0][1][4])

                                            else:

                                                if sorted_chrm_list[count-2][1][0][1][4] in remove_list:
                                                    print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-3][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                                    strand_skip = '@'

                                                elif sorted_chrm_list[count-3][1][0][1][4] not in remove_list:
                                                    print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-3][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                                    sorted_chrm_list[count-3][1][0][1][5] = int(gene[1][0][1][1])
                                                    strand_skip = '@'

                                                elif sorted_chrm_list[count-4][1][0][1][4] not in remove_list:
                                                    print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-4][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                                    sorted_chrm_list[count-4][1][0][1][5] = int(gene[1][0][1][1])
                                                    strand_skip = '@'

                                                else:
                                                    print('\n\n\nToooo much nesting - this script cannot handle the level of gene overlap in the annotation provided\n\n\n')
                                                    sys.exit()

                                        elif len(sorted_chrm_list[count-2][1]) == 2:
                                            if sorted_chrm_list[count-2][1][1][1][3] == gene[1][0][1][3]:

                                                print('\tGene ' + sorted_chrm_list[count-2][1][1][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][1][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])

                                                try:
                                                    sorted_chrm_list[count-2][1][1][1][5] = int(gene[1][0][1][1])
                                                except(IndexError):
                                                    print('\n\t!!! Missing annotation appended...\n')
                                                    sorted_chrm_list[count-2][1][1][1].append(int(gene[1][0][1][1]))
                                            else:
                                                print('\tGene ' + sorted_chrm_list[count-(2+len(sorted_chrm_list[count-2][1]))][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-(2+len(sorted_chrm_list[count-2][1]))][1][0][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])

                                                try:
                                                    sorted_chrm_list[count-(2+len(sorted_chrm_list[count-2][1]))][1][0][1][5] = int(gene[1][0][1][1])
                                                except(IndexError):
                                                    print('\n\t!!! Missing annotation appended...\n')
                                                    sorted_chrm_list[count-(2+len(sorted_chrm_list[count-2][1]))][1][0][1].append(int(gene[1][0][1][1]))

                                    print('\n')

                                    stored_start = int(gene[1][0][1][2])
                                    stored_stop  = int(gene[1][0][1][1])

                                if strand_skip == '@':
                                    stored_strand = gene[1][0][1][3]
                                else:
                                    stored_strand = strand_skip

                        elif gene[1][0][1][0] != stored_chr:

                            if gene[1][0][1][3] == '-':
                                print('Gene ' + gene[1][0][1][4] + ' stop:\t' + str(gene[1][0][1][2]) + '\treplaces with ' + gene[1][0][1][2] + ' + ' + str(CIS_WINDOW) + '\n')
                                gene[1][0][1].append(int(gene[1][0][1][2])+ CIS_WINDOW)
                                stored_start  = int(gene[1][0][1][2])
                                stored_stop   = int(gene[1][0][1][1])
                                stored_strand = gene[1][0][1][3]

                            elif gene[1][0][1][3] == '+' and int(gene[1][0][1][1]) - CIS_WINDOW <= 0:
                                gene[1][0][1].append(0)
                                print('Gene ' + gene[1][0][1][4] + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with 0')
                                stored_start = int(gene[1][0][1][1])
                                stored_stop  = int(gene[1][0][1][2])

                            else:
                                print('Gene ' + gene[1][0][1][4] + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with '+ gene[1][0][1][1] + ' - ' + str(CIS_WINDOW) + '\n')
                                gene[1][0][1].append(int(gene[1][0][1][1]) - CIS_WINDOW)
                                stored_start = int(gene[1][0][1][1])
                                stored_stop  = int(gene[1][0][1][2])

                            stored_chr = gene[1][0][1][0]

                for stored_cisRegion_coords in sorted_chrm_list:

                    if len(stored_cisRegion_coords[1]) > 1:

                        for duplicate_cisRegion_coords in stored_cisRegion_coords[1]:

                            if duplicate_cisRegion_coords[1][4] not in remove_list:

                                # print('Keeping: ' + str(duplicate_cisRegion_coords))

                                if duplicate_cisRegion_coords[1][3] == '+':

                                    if int(duplicate_cisRegion_coords[1][1]) - int(duplicate_cisRegion_coords[1][5]) != 0:
                                        #out_file.write(duplicate_cisRegion_coords[1][0] + '\t' + str(int(duplicate_cisRegion_coords[1][5])) + '\t' + str(duplicate_cisRegion_coords[1][1]) + '\t' + str(duplicate_cisRegion_coords[0]) + '\t1\t' + str(duplicate_cisRegion_coords[1][3]) + '\n')
                                        outputFormat(out_file, duplicate_cisRegion_coords[1][0], duplicate_cisRegion_coords[1][5], duplicate_cisRegion_coords[1][1], duplicate_cisRegion_coords[0], duplicate_cisRegion_coords[1][3])

                                    else:

                                        print('--- ' + str(duplicate_cisRegion_coords[0]) + ' removed as ZERO LENGTH' + '\t' + duplicate_cisRegion_coords[1][3] + '\t' + duplicate_cisRegion_coords[1][1] + '\t' + duplicate_cisRegion_coords[1][2])
                                        remove_list.append(duplicate_cisRegion_coords[0])

                                elif duplicate_cisRegion_coords[1][3] == '-':

                                    if int(duplicate_cisRegion_coords[1][5]) > int(SCAFF_LIMS[duplicate_cisRegion_coords[1][0]]):

                                        if int(SCAFF_LIMS[duplicate_cisRegion_coords[1][0]]) - int(duplicate_cisRegion_coords[1][2]) != 0:

                                            print('\n\n\n END OF SCAFF REACHED  -- CLIPPING STOP ANNOT FROM ' + str(duplicate_cisRegion_coords[1][5]) + ' to ' + str(SCAFF_LIMS[duplicate_cisRegion_coords[1][0]]) + '\n\n\n')

                                            # out_file.write(duplicate_cisRegion_coords[1][0] + '\t' + str(int(duplicate_cisRegion_coords[1][2])) + '\t' + str(SCAFF_LIMS[duplicate_cisRegion_coords[1][0]]) + '\t' + str(duplicate_cisRegion_coords[0]) + '\t1\t' + str(duplicate_cisRegion_coords[1][3]) + '\n')
                                            outputFormat(out_file, duplicate_cisRegion_coords[1][0], duplicate_cisRegion_coords[1][2], str(SCAFF_LIMS[duplicate_cisRegion_coords[1][0]]), duplicate_cisRegion_coords[0], duplicate_cisRegion_coords[1][3])
                                        else:
                                            print('--- ' + str(duplicate_cisRegion_coords[0]) + ' removed as ZERO LENGTH AND END OF SCAFF' + '\t' + duplicate_cisRegion_coords[1][3] + '\t' + duplicate_cisRegion_coords[1][1] + '\t' + duplicate_cisRegion_coords[1][2])
                                            remove_list.append(duplicate_cisRegion_coords[0])

                                    elif int(duplicate_cisRegion_coords[1][5]) - int(duplicate_cisRegion_coords[1][2]) != 0:
                                        # out_file.write(duplicate_cisRegion_coords[1][0] + '\t' + str(int(duplicate_cisRegion_coords[1][2])) + '\t' + str(duplicate_cisRegion_coords[1][5]) + '\t' + str(duplicate_cisRegion_coords[0]) + '\t1\t' + str(duplicate_cisRegion_coords[1][3]) + '\n')
                                        outputFormat(out_file, duplicate_cisRegion_coords[1][0], duplicate_cisRegion_coords[1][2], duplicate_cisRegion_coords[1][5], duplicate_cisRegion_coords[0], duplicate_cisRegion_coords[1][3])

                                    else:
                                        print('--- ' + str(duplicate_cisRegion_coords[0]) + ' removed as ZERO LENGTH' + '\t' + duplicate_cisRegion_coords[1][3] + '\t' + duplicate_cisRegion_coords[1][1] + '\t' + duplicate_cisRegion_coords[1][2])
                                        remove_list.append(duplicate_cisRegion_coords[0])

                            else:
                                print('--- ' + str(duplicate_cisRegion_coords[1][4]) + ' removed' + '\t' + duplicate_cisRegion_coords[1][3] + '\t' + duplicate_cisRegion_coords[1][1] + '\t' + duplicate_cisRegion_coords[1][2])
                    else:

                        if stored_cisRegion_coords[1][0][1][4] not in remove_list:

                            # print('Keeping: ' + str(stored_cisRegion_coords))

                            if stored_cisRegion_coords[1][0][1][3] == '+':

                                if int(stored_cisRegion_coords[1][0][1][1]) - int(stored_cisRegion_coords[1][0][1][5]) != 0:
                                    # out_file.write(str(stored_cisRegion_coords[1][0][1][0]) + '\t' + str(stored_cisRegion_coords[1][0][1][5]) + '\t' + str(stored_cisRegion_coords[1][0][1][1]) + '\t' + str(stored_cisRegion_coords[1][0][1][4].strip()) + '\t1\t' + str(stored_cisRegion_coords[1][0][1][3].strip()) + '\n')
                                    outputFormat(out_file, stored_cisRegion_coords[1][0][1][0], stored_cisRegion_coords[1][0][1][5], stored_cisRegion_coords[1][0][1][1], stored_cisRegion_coords[1][0][1][4].strip(), stored_cisRegion_coords[1][0][1][3].strip())

                                else:

                                    print('--- ' + str(stored_cisRegion_coords[1][0][1][4]) + ' removed as ZERO LENGTH' + '\t' + str(stored_cisRegion_coords[1][0][1][3]) + '\t' + str(stored_cisRegion_coords[1][0][1][1]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]))


                            elif stored_cisRegion_coords[1][0][1][3] == '-':

                                if int(stored_cisRegion_coords[1][0][1][5]) > int(SCAFF_LIMS[stored_cisRegion_coords[1][0][1][0]]):

                                    if int(stored_cisRegion_coords[1][0][1][2]) - int(SCAFF_LIMS[stored_cisRegion_coords[1][0][1][0]]) != 0:
                                        print('\n\n\n END OF SCAFF REACHED  -- CLIPPING STOP ANNOT FROM ' + str(stored_cisRegion_coords[1][0][1][5]) + ' to ' + str(SCAFF_LIMS[stored_cisRegion_coords[1][0][1][0]]) + '\n\n\n')
                                        # out_file.write(str(stored_cisRegion_coords[1][0][1][0]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]) + '\t' + str(SCAFF_LIMS[stored_cisRegion_coords[1][0][1][0]]) + '\t' + str(stored_cisRegion_coords[1][0][1][4].strip()) + '\t1\t' + str(stored_cisRegion_coords[1][0][1][3].strip()) + '\n')
                                        outputFormat(out_file, stored_cisRegion_coords[1][0][1][0], stored_cisRegion_coords[1][0][1][2], SCAFF_LIMS[stored_cisRegion_coords[1][0][1][0]], stored_cisRegion_coords[1][0][1][4].strip(), stored_cisRegion_coords[1][0][1][3].strip())

                                    else:
                                        print('--- ' + str(stored_cisRegion_coords[1][0][1][4]) + ' removed as ZERO LENGTH AT END OF SCAFF' + '\t' + str(stored_cisRegion_coords[1][0][1][3]) + '\t' + str(stored_cisRegion_coords[1][0][1][1]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]))
                                        remove_list.append(stored_cisRegion_coords[1][0][1][4])

                                elif int(stored_cisRegion_coords[1][0][1][2]) - stored_cisRegion_coords[1][0][1][5] != 0:
                                    # out_file.write(str(stored_cisRegion_coords[1][0][1][0]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]) + '\t' + str(stored_cisRegion_coords[1][0][1][5]) + '\t' + str(stored_cisRegion_coords[1][0][1][4].strip()) + '\t1\t' + str(stored_cisRegion_coords[1][0][1][3].strip()) + '\n')
                                    outputFormat(out_file, stored_cisRegion_coords[1][0][1][0], stored_cisRegion_coords[1][0][1][2], stored_cisRegion_coords[1][0][1][5], stored_cisRegion_coords[1][0][1][4].strip(), stored_cisRegion_coords[1][0][1][3].strip())

                                else:

                                    print('--- ' + str(stored_cisRegion_coords[1][0][1][4]) + ' removed as ZERO LENGTH' + '\t' + str(stored_cisRegion_coords[1][0][1][3]) + '\t' + str(stored_cisRegion_coords[1][0][1][1]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]))
                                    remove_list.append(stored_cisRegion_coords[1][0][1][4])
                        else:
                            print('--- ' + str(stored_cisRegion_coords[1][0][1][4]) + ' removed' + '\t' + str(stored_cisRegion_coords[1][0][1][3]) + '\t' + str(stored_cisRegion_coords[1][0][1][1]) + '\t' + str(stored_cisRegion_coords[1][0][1][2]))

    # with open(os.path.join(OUTPUT_DIR, filename[:-3] + 'same_strand+start.out'), 'w') as out2:
    with open(os.path.join(OUTPUT_DIR, 'same_strand+start.out'), 'w') as out2:
        for out in same_strand_start:
            if out not in remove_list:
                out2.write(str(out) + '\n')
    # print('\n\n\nOutputs:\n\t' + os.path.join(OUTPUT_DIR, filename[:-3] + str(CIS_WINDOW) + 'nt_cisRegions.stranded.bed') + '\n\t' + os.path.join(OUTPUT_DIR, filename[:-3] + 'same_strand+start.out') + '\n\n')
    print('\n\n\nOutputs:\n\t' + os.path.join(OUTPUT_DIR, 'cisRegions.bed') + '\n\t' + os.path.join(OUTPUT_DIR, 'same_strand+start.out') + '\n\n')


def main(GENE_BED, GENOME_FASTA, ntWINDOW, OUTPUT_DIR):
    print('\n\ncisRegion.py\n\n')
    print('Loading annotations from:\t' + GENE_BED)
    print('Loading sequences from:\t'   + GENOME_FASTA)
    print('Extracting cisRegions of:\t' + GENOME_FASTA + 'nt')
    print('Writing cis annotations to:\t' + OUTPUT_DIR + '\n')

    # create a list of scaffold IDs
    scaffold_set = scaffoldLister(GENE_BED)

    # create a dictionary of the lengths of all scaffolds
    scaff_limits = chromLimitFind(scaffold_set, genomeLoad(GENOME_FASTA))

    #call cis regions - to be updated
    hugeCisRegionCallingBehemoth(scaffold_set, GENE_BED, OUTPUT_DIR, scaff_limits, ntWINDOW)

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('BED_in',   type=str, help='path to cisRegion annotation in .BED format')
parser.add_argument('FASTA_in', type=str, help='path to genome sequences in .FASTA format')
parser.add_argument('WINDOW',   type=int, help='integer value describing nucleotide length of max cis region to be extracted')
parser.add_argument('OUT_dir',  type=str, help='path to directory where the output should be written')

if __name__ == '__main__':

    args = parser.parse_args()

    main(args.BED_in, args.FASTA_in, args.WINDOW, args.OUT_dir)

__author__ = "Will Nash"
__copyright__ = "Copyright 2020, The Earlham Institute"
__credits__ = ["Will Nash", "Wilfried Haerty"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Will Nash"
__email__ = "will.nash@earlham.ac.uk"
__status__ = "Testing"
