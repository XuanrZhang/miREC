#!/bin/bash

if [ "$1" == "mix" ]
then
	wc -l $2;
	for i in $(seq $3 $4)
	do
	echo “$i”;
	done
else
	wc $2;
fi
