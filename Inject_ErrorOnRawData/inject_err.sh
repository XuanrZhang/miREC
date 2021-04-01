#!/bin/bash

S=1;
O=output.fq;

while getopts f:s:o: op
do 
    case $op in
        f)
            echo "Input sequence file name is: $OPTARG"
            F=$OPTARG;;
        s)
            echo "seed number for error generation is: $OPTARG"
            S=$OPTARG;;
        o)
            echo "Output file name with errors is : $OPTARG"
            O=$OPTARG;;
        \?)
            echo "Usage: args [-f] [-s] [-d] [-o] [-g] [-t]"
            echo "-f means Input sequence fastq file name "
            echo "-s means seed number for error generation"
            echo "-o means Output file name with errors (default: simulated.fa) "
            exit 1;;
    esac
done


awk '{if((NR%2)==1)print $1;else print $0}' $F > ./input.fq;
awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' input.fq | awk '{print $1 " " $(NF-2) " " $NF}' > ./id_read.txt;
awk '{print $2}' ./id_read.txt |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor.txt;

./induce_err -s $S -c expreLevel_cor.txt -f id_read.txt -o $O;
