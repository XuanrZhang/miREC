CC= g++
CFLAG = -O3 -Wall -std=c++0x

test: miREC.c


$(CC) -fopenmp  miREC_mix_fq.c -o miREC_mix_fq $(CFLAG)
$(CC) -fopenmp  miREC_fq.c -o miREC_fq $(CFLAG)
# $(CC) -fopenmp miREC.c -o miREC $(CFLAG)

# $(CC) find_common.c -o find_common $(CFLAG)
# $(CC) miREC_subindel.c -o miREC_subindel $(CFLAG)
# $(CC) find_isomir.c -o find_isomir $(CFLAG)
# $(CC) simu_tofastq.c -o simu_tofastq $(CFLAG)
# $(CC) -fopenmp  copy_change.c -o copy_change $(CFLAG)
# $(CC) -fopenmp  miREC_mix.c -o miREC_mix $(CFLAG)

# $(CC) miREC_indel.c -o miREC_indel $(CFLAG)
# $(CC) -fopenmp  miREC_indel.c -o miREC_indel $(CFLAG)
# $(CC) -fopenmp  miREC_indel_fq.c -o miREC_indel_fq $(CFLAG)
# $(CC) miREC_fq_update.c -o miREC_fq_update $(CFLAG)
# $(CC) miREC.c -o miREC $(CFLAG)
# $(CC)  miREC_update.c -o miREC_update $(CFLAG)
# $(CC) gene_simu_mixerr.c -o gene_simu_mixerr $(CFLAG)
# $(CC) test_h.c -o test_h $(CFLAG)
# $(CC) miREC_fq.c -o miREC_fq $(CFLAG)
# $(CC) miREC_in.c -o miREC_in $(CFLAG)
# $(CC) miREC_infq.c -o miREC_infq $(CFLAG)
# $(CC) miREC_de.c -o miREC_de $(CFLAG)
# $(CC) miREC_defq.c -o miREC_defq $(CFLAG)
# $(CC) gene_simu_sub.c -o gene_simu_sub $(CFLAG)
# $(CC) gene_simu.c -o gene_simu $(CFLAG)
# $(CC) gene_simu_sub.c -o gene_simu_sub $(CFLAG)
# $(CC) gene_simu_del.c -o gene_simu_del $(CFLAG)
# $(CC) gene_simu_ins.c -o gene_simu_ins $(CFLAG)
# $(CC) check.c -o check $(CFLAG)
# $(CC) check_fa.c -o check_fa $(CFLAG)

