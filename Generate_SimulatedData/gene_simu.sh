#!/bin/bash

D=distrubution.txt;
F=mature.fa;
S=1;
O=simulated.fa;
G=truth.fa;
T=0;

while getopts d:f:s:o:g:t op
do 
    case $op in
        d)
            echo "Input distrubution file name is: $OPTARG"
            D=$OPTARG;;
        f)
            echo "Input sequence template file name is: $OPTARG"
            F=$OPTARG;;
        s)
            echo "seed number for error generation is: $OPTARG"
            S=$OPTARG;;
        o)
            echo "Output file name with errors is : $OPTARG"
            O=$OPTARG;;
        g)
            echo "Output file name is : $OPTARG"
            G=$OPTARG;;
        t) 
            echo "error type (mix or subs only)"
            T=1;;
        \?)
            echo "Usage: args [-f] [-s] [-d] [-o] [-g] [-t]"
            echo "-f means Input distrubution file name "
            echo "-d means Input sequence template file name "
            echo "-s means seed number for error generation"
            echo "-o means Output file name with errors (default: simulated.fa) "
            echo "-g means Output file name (default: truth.fa)"
            echo "-t means subs error only"
            exit 1;;
    esac
done


if [ $T -eq 1 ]
then
    echo "start generate dataset with subs error only ----";
    echo "./gene_simu_sub -s ${S} -c ${D} -f ${F} -o ${O} -g ${G}";
    ./gene_simu_sub -s ${S} -c ${D} -f ${F} -o ${O} -g ${G};
else
    echo "start generate dataset with mixed error  ----";
    echo "./gene_simu_mixerr -s ${S} -c ${D} -f ${F} -o ${O} -g ${G}";
    ./gene_simu_mixerr -s ${S} -c ${D} -f ${F} -o ${O} -g ${G};
fi

