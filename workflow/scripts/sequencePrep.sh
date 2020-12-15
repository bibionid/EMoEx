#!/bin/bash -e
##
## sequencePrep.sh
##
## Converts multi-line fasta to single line fasta, then simplifies the header to
## a single term by splitting on the spcace character. Single line version is
## then remove.
##
##
## Will Nash
## Github: https://github.com/bibionid/
## will.nash@earlham.ac.uk

# Recieve genome from first commandline argument
INPUT_GENOME=$1

# Use awk to convert fasta from multi to single line
awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${INPUT_GENOME} > ${INPUT_GENOME}.singleLine

# Use awk to simplify the headers in the single line fasta
awk -F" " '/^>/{print $1; next}{print}' < ${INPUT_GENOME}.singleLine > ${INPUT_GENOME}.singleLine.simpleHeader

# Remove the intermediary file
rm ${INPUT_GENOME}.singleLine
