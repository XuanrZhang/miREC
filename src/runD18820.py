'''
Author: Pengyao PING
Date: 2021-06-13 12:52:48
LastEditors: Pengyao PING
LastEditTime: 2021-06-13 13:00:04
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
import time
from multiprocessing import Process
from Utils import Utils

def main():
    start = time.time()
    miRNA_file = "./data/input/963/GSE139936_180719_GEO_miRNAs.txt"
    #output_dir = "../data/output/"
    input_dir = "/home/pping/Data/miREC/data/input/"
    output_dir = "/home/pping/Data/miREC/data/output/countD18_8_20/"

    U = Utils()

    #remove NN 
    
    #out_dir1 = input_dir + "NN_removed_raw_fq/"
    #in_dir1 = input_dir + "raw_fq/"
    #U.mutliprocess_remove_NN(out_dir1, in_dir1)
    #U.mutliprocess_count(miRNA_file, output_dir + "count/NN_removed_raw_fq/", out_dir1)

    #filter_NN_removed_fq_dir = input_dir + 'filtered_NN_removed_raw_fq/'
    #U.mutliprocess_filter(miRNA_file, out_dir1, filter_NN_removed_fq_dir)
    #U.mutliprocess_count(miRNA_file, output_dir + "count/filtered_NN_removed_raw_fq/", filter_NN_removed_fq_dir)

    #corrected_fq_dir = input_dir + 'corrected/filtered_NN_removed_raw_fq/'
    #U.multiprocess_miREC_call(filter_NN_removed_fq_dir, corrected_fq_dir, 8, 15, 10)

    #filter_NN_removed_fq_dir2 = input_dir + 'filtered_NN_removed_raw_fq2/'
    #inter_filter_NN_removed_fq_dir = input_dir + 'inter_filtered_NN_removed_raw_fq/'
    #U.multiprocess_forward_reverse_inter(filter_NN_removed_fq_dir2, inter_filter_NN_removed_fq_dir)

    input_dir1 = input_dir + 'filtered_NN_removed_raw_fq/'
    input_dir2 = input_dir + 'corrected/D18_miREC_8_20/'
    input_dir3 = input_dir + 'corrected/filtered_NN_removed_raw_fq_karect/'
    U.mutliprocess_multi_count(miRNA_file, output_dir, input_dir1, input_dir2, input_dir3)
    #asm_fq_dir = input_dir + 'asm_fq/'
    #U.multiprocess_pear_call(filter_NN_removed_fq_dir, asm_fq_dir, threads=10, n=0)

    '''
    input_dir = './data/input/'
    output_dir = "./data/filtered_input/"
    U.mutliprocess_filter(miRNA_file, input_dir, output_dir)

    #input_dir = "./data/filtered_input/"
    #output_dir = "./data/output/"
    #U.mutliprocess_count(miRNA_file, output_dir, input_dir)
    '''
    '''
    #ff = input_dir + '180719Ded_D18-6962_1_sequence.10bp.fasta'
    #fr = input_dir + '180719Ded_D18-6962_2_sequence.10bp.fasta'
    ff = './mit_data/remote/dataset2a/180719Ded_D18-6962_1_sequence.3clipped.fq'
    fr = './mit_data/remote/dataset2a/180719Ded_D18-6962_2_sequence.rc3clipped.fq'
    fout = "./data/mid_input/asm_input/D18.6962"
    threads = 20
    n = 0
    U.pear_call(ff, fr, fout, threads, n)
    '''
    end = time.time()
    print("Time:{}".format(end - start))
    return

if __name__ == '__main__':
    main()