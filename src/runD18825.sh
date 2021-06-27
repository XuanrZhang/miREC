#!/bin/bash

Output_dir=./Data/synthetic_963miRNAs_Reads/corrected/D18_miREC_8_25/
input_dir=./Data/synthetic_963miRNAs_Reads/filtered_NN_removed_raw_fq/
files=$(ls ./Data/synthetic_963miRNAs_Reads/filtered_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./miREC.sh -f $input_dir$file -s 8 -e 25 -t 20 -o $Output_dir$file
done