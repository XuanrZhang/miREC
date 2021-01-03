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
    awk '{if((NR%2)==1)print $1;else print $0}' ${F} > input.fq
    awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' input.fq | awk '{print $1 " " $(NF-2) " " $NF}' > ./id_read.txt;
    awk '{print $2}' ./id_read.txt |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor.txt   
    cp ./id_read.txt ./ID_read_quality_cor.txt
    cp ./ID_read_quality_cor.txt ./ID_read_quality_input.txt
    cp input.fq ./correct_read.fastq
    
    for i in $(seq $S $E )
    do
        #----recurring prepare ----: finish error correction, recount 'kmer frequency' and 'read frequency' and 'id_read'
        #recount 'kmer frequency', then create kmer.freq (e.g. 5mer.freq)
        ./kmc -k${i} -fq -ci1 ./correct_read.fastq tmp${i} ./
        ./kmc_dump tmp${i} tmpkc${i}
        sort -nk2 -r tmpkc${i} > ./${i}mer.freq
        rm tmp*

        #recount 'read frequency', then create read_expresslevel data ([read_freq] [read])
        awk '{print $2}' ID_read_quality_cor.txt |sort |uniq -c| sort -r -nk1 > expreLevel_cor.txt
        echo "----------------------${i} mer frequency preparation ready";

        #error correction
        #echo "./miREC_fq -k ${i} -m /home/xuanzhan/Data/miRNA/simu/${i}mer.freq -l expreLevel_cor.txt -f ID_read_quality_input.txt >> miREC_subindel${file_id[${j}]}.log;"

        ./miREC_fq -k ${i} -m ./${i}mer.freq -l expreLevel_cor.txt -f ID_read_quality_input.txt >> miREC_subindel${file_id[${j}]}.log;

        cp ID_read_quality_cor.txt ID_read_quality_input.txt 
        #cp correct_read.fa correct_read_cp.fa
    done
    rm *.freq *.txt

else
    echo "running mix error correction";
    awk '{if((NR%2)==1)print $1;else print $0}' ${F} > input.fq
    awk '{if(NR%4!=0)ORS=" ";else ORS="\n"}1' input.fq | awk '{print $1 " " $(NF-2) " " $NF}' > ./id_read.txt;
    awk '{print $2}' ./id_read.txt |sort |uniq -c| sort -r -nk1 > ./expreLevel_cor.txt   
    cp ./id_read.txt ./ID_read_quality_cor.txt
    cp ./ID_read_quality_cor.txt ./ID_read_quality_input.txt
    cp input.fq ./correct_read.fastq
    
    for i in $(seq $S $E )
    do
        #----recurring prepare ----: finish error correction, recount 'kmer frequency' and 'read frequency' and 'id_read'
        #recount 'kmer frequency', then create kmer.freq (e.g. 5mer.freq)
        ./kmc -k${i} -fq -ci1 ./correct_read.fastq tmp${i} ./
        ./kmc_dump tmp${i} tmpkc${i}
        sort -nk2 -r tmpkc${i} > ./${i}mer.freq
        rm tmp*

        #recount '(k-1)mer frequency', then create kmer.freq (e.g. 5mer.freq)
        tmpm=$(($i-1))
        ./kmc -k${tmpm} -fq -ci1 ./correct_read.fastq tmp${tmpm} ./
        ./kmc_dump tmp${tmpm} tmpkc${tmpm}
        sort -nk2 -r tmpkc${tmpm} > ./${tmpm}mer.freq
        rm tmp*

        #recount '(k+1)mer frequency', then create kmer.freq (e.g. 5mer.freq)
        tmpa=$(($i+1))
        ./kmc -k${tmpa} -fq -ci1 ./correct_read.fastq tmp${tmpa} ./
        ./kmc_dump tmp${tmpa} tmpkc${tmpa}
        sort -nk2 -r tmpkc${tmpa} > ./${tmpa}mer.freq
        rm tmp*

        #recount 'read frequency', then create read_expresslevel data ([read_freq] [read])
        awk '{print $2}' ID_read_quality_cor.txt |sort |uniq -c| sort -r -nk1 > expreLevel_cor.txt
        echo "----------------------${i} mer frequency preparation ready";

        #error correction
        #echo "./miREC_mix_fq -k ${i} -m /home/xuanzhan/Data/miRNA/simu/${i}mer.freq -s /home/xuanzhan/Data/miRNA/simu/${tmpm}mer.freq -b /home/xuanzhan/Data/miRNA/simu/${tmpa}mer.freq -l expreLevel_cor.txt -f ID_read_quality_input.txt >> miREC_mix${file_id[${j}]}.log;"

        # time ./miREC_fq_update -k ${i} -m /home/xuanzhan/Data/miRNA/simu/${i}mer.freq -l expreLevel_cor.txt -f ID_read_quality_input.txt >> miREC_mix${file_id[${j}]}.log;
        ./miREC_mix_fq -k ${i} -m ${i}mer.freq -s ${tmpm}mer.freq -b ${tmpa}mer.freq -l expreLevel_cor.txt -f ID_read_quality_input.txt >> miREC_mix${file_id[${j}]}.log;

        echo "----------------------${i} mer correction finished";

        cp ID_read_quality_cor.txt ID_read_quality_input.txt 
        #cp correct_read.fa correct_read_cp.fa
    done
    rm *.freq *.txt
    
fi
