#!/bin/bash

input_dir=./Data/963miRXplore_data/D18_filtered_NN_removed_raw_fq/
files=$(ls ./Data/963miRXplore_data/D18_filtered_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./Verify_963miRXploreData/karect -correct -threads=12 -matchtype=hamming -celltype=haploid -inputfile=$input_dir$file 
    mv ./karect_* ./Data/963miRXplore_data/corrected/D18_karect/
    rm *.txt ./temp_res*
done
