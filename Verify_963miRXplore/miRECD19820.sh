#!/bin/bash

Output_dir=./Data/Verify_by963miRXploreData/corrected/D19_miREC_8_20/
input_dir=./Data/Verify_by963miRXploreData/D19_NN_removed_raw_fq/
files=$(ls ./Data/Verify_by963miRXploreData/D19_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./miREC.sh -f $input_dir$file -s 8 -e 20 -t 20 -o $Output_dir$file
done