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
    input_dir = "./data/input/"
    output_dir = "./data/output/corrected"

    U = Utils()

    filter_NN_removed_fq_dir = input_dir + 'filtered_NN_removed_raw_fq2/'
    U.multiprocess_miREC_call(filter_NN_removed_fq_dir, output_dir, 8, 15, 8)
    end = time.time()
    print("Time:{}".format(end - start))
    return

if __name__ == '__main__':
    main()