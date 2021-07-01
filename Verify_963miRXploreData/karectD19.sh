#!/bin/bash

input_dir=./Data/963miRXplore_data/D19_NN_removed_raw_fq/
files=$(ls ./Data/963miRXplore_data/D19_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./Verify_963miRXploreData/karect -correct -threads=12 -matchtype=hamming -celltype=haploid -inputfile=$input_dir$file 
    # move the corrected files to data directory
    mv ./karect_* ./Data/963miRXplore_data/corrected/D19_karect/
    rm *.txt ./temp_res*
done
