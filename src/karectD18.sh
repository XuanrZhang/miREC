#!/bin/bash

input_dir=./Data/Verify_by963miRXploreData/D18_filtered_NN_removed_raw_fq/
files=$(ls ./Data/Verify_by963miRXploreData/D18_filtered_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./karect -correct -threads=12 -matchtype=hamming -celltype=haploid -inputfile=$input_dir$file 
    #./miREC.sh -f $input_dir$file -s 8 -e 15 -t 20 -o $Output_dir$file
done
