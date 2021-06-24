#!/bin/bash

Output_dir=/home/pping/Data/miREC/data/input/corrected/D19_miREC_8_25/
input_dir=/home/pping/Data/miREC/data/input/D19_NN_removed_raw_fq/
files=$(ls /home/pping/Data/miREC/data/input/D19_NN_removed_raw_fq/)
for file in $files
do
    echo $input_dir$file
    ./miREC.sh -f $input_dir$file -s 8 -e 25 -t 20 -o $Output_dir$file
done