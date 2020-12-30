#!/bin/bash

R=0;
T=5;
S=15;
E=20;

while getopts f:t:s:e:o op
do 
    case  $op in
        f)
            echo "Input file name is: $OPTARG"
            F=$OPTARG;;

        t)
            echo "Threshod value is: $OPTARG"
            T=$OPTARG;;
        s)
            echo "K_1 value is: $OPTARG"
            S=$OPTARG;;
        e)
            echo "K_end value is: $OPTARG"
            E=$OPTARG;;
        o)
            echo "correct subs error only"
            R=1;;
        \?)
            echo "Usage: args [-f] [-s] [-e] [-t] [-o]"
            echo "-f means Input file name "
            echo "-t means Threshod value"
            echo "-s means k_1 value"
            echo "-e means k_end value"
            echo "-o means run_type is subs error only"
            exit 1;;
    esac
done


echo "$F $T $E $S $R";

if [ $R -eq 1 ]
then
    echo "running subs error correction only";
else
    echo "running mix error correction";
fi



-----
	wc -l $2;
	for i in $(seq $3 $4)
	do
	echo “$i”;
	done

