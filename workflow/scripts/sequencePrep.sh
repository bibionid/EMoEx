#!/bin/bash -e

INPUT_GENOME=$1

awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ${INPUT_GENOME} > ${INPUT_GENOME}.singleLine

awk -F" " '/^>/{print $1; next}{print}' < ${INPUT_GENOME}.singleLine > ${INPUT_GENOME}.singleLine.simpleHeader
