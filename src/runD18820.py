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
    miRNA_file = "./Data/synthetic_963miRNAs_Reads/963miRNAs/GSE139936_180719_GEO_miRNAs.txt"
    input_dir = "./Data/synthetic_963miRNAs_Reads/"
    output_dir = "./Data/synthetic_963miRNAs_Reads/output/countD18_8_20/"

    U = Utils()

    #remove NN 
    out_dir1 = input_dir + "NN_removed_raw_fq/"
    in_dir1 = input_dir + "raw_fq/"
    U.mutliprocess_remove_NN(out_dir1, in_dir1)
    #U.mutliprocess_count(miRNA_file, output_dir + "count/NN_removed_raw_fq/", out_dir1)

    # filter noise sequences for datasets included in the directory of NN_removed_raw_fq
    filter_NN_removed_fq_dir = input_dir + 'filtered_NN_removed_raw_fq/'
    U.mutliprocess_filter(miRNA_file, out_dir1, filter_NN_removed_fq_dir)
    # The following code only count copy numbers for the data without correction
    #U.mutliprocess_count(miRNA_file, output_dir + "count/filtered_NN_removed_raw_fq/", filter_NN_removed_fq_dir)

    input_dir1 = input_dir + 'filtered_NN_removed_raw_fq/'
    input_dir2 = input_dir + 'corrected/D18_miREC_8_20/'
    input_dir3 = input_dir + 'corrected/D18_Karect/'
    U.mutliprocess_multi_count(miRNA_file, output_dir, input_dir1, input_dir2, input_dir3)

    end = time.time()
    print("Time:{}".format(end - start))
    return

if __name__ == '__main__':
    main()