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
    #output_dir = "../data/output/"
    input_dir = "./Data/synthetic_963miRNAs_Reads/"
    output_dir = "./Data/synthetic_963miRNAs_Reads/output/countD18_8_25/"

    U = Utils()

    input_dir1 = input_dir + 'filtered_NN_removed_raw_fq/'
    input_dir2 = input_dir + 'corrected/D18_miREC_8_25/'
    input_dir3 = input_dir + 'corrected/D18_Karect/'
    U.mutliprocess_multi_count(miRNA_file, output_dir, input_dir1, input_dir2, input_dir3)


    end = time.time()
    print("Time:{}".format(end - start))
    return

if __name__ == '__main__':
    main()