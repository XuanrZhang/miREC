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
import subprocess

def main():
    start = time.time()

    ###############################################################################################
    # To create some directories for saving dataset
    try:
        subprocess.call(["./Verify_963miRXploreData/mkdir.sh"], shell=False)
    except:
        print('You may need using "chmod +x ./Verify_963miRXploreData/mkdir.sh" to make it has permission to run')

    ###############################################################################################
    miRNA_file = "./Data/963miRXplore_data/963miRNAs/GSE139936_180719_GEO_miRNAs.txt"
    input_dir = "./Data/963miRXplore_data/"
    output_dir = "./Data/963miRXplore_data/output/countD19/"

    U = Utils()

    input_dir1 = input_dir + 'D19_NN_removed_raw_fq/'
    input_dir2 = input_dir + 'corrected/D19_miREC/'
    input_dir3 = input_dir + 'corrected/D19_karect/'
    
    ###############################################################################################
    # Do errorection by miREC
    try:
        subprocess.call(["./Verify_963miRXploreData/miRECD19.sh"], shell=False)
    except:
        print('Please make sure miREC work well')

    ###############################################################################################
    # Do errorection by Karect
    try:
        subprocess.call(["./Verify_963miRXploreData/karectD19.sh"], shell=False)
    except:
        print('Please make sure karect work well')
    
    ##############################################################################################
    U.mutliprocess_multi_count(miRNA_file, output_dir, input_dir1, input_dir2, input_dir3)

    end = time.time()
    print("Time:{}".format(end - start))
    return

if __name__ == '__main__':
    main()