#!/bin/bash -e
##
## annotPrep.sh
##
## extracts 'gene' lines from whole genome annotations in gff/gff3 format. These
## lines are then searched for 'biotype=protein' annotations. The format is them
## converted to a modified .bed through a series of awk and sed calls.
##
## The output of this script is 6 column .bed format with a 7th column describing
## describing gene name.
##
## Will Nash
## Github: https://github.com/bibionid/
## will.nash@earlham.ac.uk


INPUT_ANNOTATION=$1

if [[ "${INPUT_ANNOTATION}" == *.gff ]]
then

  printf "\ngff to bed\n\n"
  awk '$3 == "gene"' inputs/ceratitis/GCF_000347755.3_Ccap_2.1_genomic.gff | grep 'biotype=protein_coding' | awk -F";" 'OFS="\t" {print $1,$3}' | sed 's/ID=gene-//g' | sed 's/Name=//g' | awk 'OFS="\t" {print $1,$4,$5,$9,$6,$7,$10}' > ${INPUT_ANNOTATION}.protein_coding.bed

elif [[ "${INPUT_ANNOTATION}" == *.gff3 ]]
then

  printf "\ngff3 to bed\n\n"
  awk '$3 == "gene"' inputs/ceratitis/GCF_000347755.3_Ccap_2.1_genomic.gff | grep 'biotype=protein_coding' | awk -F";" 'OFS="\t" {print $1,$2}' | sed 's/ID=gene://g' | sed 's/Name=//g' | awk 'OFS="\t" {print $1,$4,$5,$9,$6,$7,$10}' > ${INPUT_ANNOTATION}.protein_coding.bed

else

  printf "\nthis parser was not designed for this format\n\n"

fi
