#!/bin/bash

Output_dir=./Data/963miRXplore_data/corrected/D19_miREC/
input_dir=./Data/963miRXplore_data/D19_NN_removed_raw_fq/
files=$(ls ./Data/963miRXplore_data/D19_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./miREC.sh -f $input_dir$file -s 8 -e 25 -t 20 -o $Output_dir$file
    rm ./*.freq ./*.txt ./input.fq
done