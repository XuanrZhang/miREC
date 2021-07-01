#!/bin/bash

Output_dir=./Data/963miRXplore_data/corrected/D18_miREC_8_20/
input_dir=./Data/963miRXplore_data/D18_filtered_NN_removed_raw_fq/
files=$(ls ./Data/963miRXplore_data/D18_filtered_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./miREC.sh -f $input_dir$file -s 8 -e 20 -t 20 -o $Output_dir$file
done