'''
Author: Pengyao PING
Date: 2021-06-10 14:23:32
LastEditors: Pengyao PING
LastEditTime: 2021-06-13 11:54:06
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from collections import Counter
import editdistance
import os
import multiprocessing
import subprocess

class Utils:
    def __init__(self):
        pass

    def selection(self, fin, fout):
        df = pd.read_table(fin, sep = '\t', header=None)
        df_select = df[df[2] <5]
        df_select.to_csv(fout, sep='\t', header=None, index=None)
        return df_select

    def add_seqs(self, fin1, fout1, fin2, fout2):
        mirnas_df = self.selection(fin1, fout2)
        record_iterator = SeqIO.parse(fin2, "fa")
        records = []
        for i in mirnas_df.index:
            item = mirnas_df.loc[i]
            #print(len(item[1]))
            #print(len(item[0]))
            for j in range(1,6):
                qua = {}
                qua['phred_quality'] = [40] * len(str(item[1]))
                iddes = str(item[0] + "addition.{}".format(j))
                seqq = Seq(str(item[1]))
                temp_rec = SeqRecord(seqq, id=iddes, description=iddes, letter_annotations=qua)
                records.append(temp_rec)
        for item in record_iterator:
            qual = {}
            val = item.letter_annotations['phred_quality']
            qual['phred_quality'] = val
            tmp_seq = str(item.seq)
            tmp_id = item.id
            tmp_ds = item.description
            #print(type(Seq(tmp_seq)))
            tmp_rec = SeqRecord(Seq(tmp_seq), id=tmp_id, description=tmp_ds, letter_annotations=qual)
            records.append(tmp_rec)
        SeqIO.write(records, fout2, "fastq")
        return

    # remove the last two bases for each sequence in fastq file
    def remove_NN(self, filename, output_fastq):
        record_iterator = SeqIO.parse(filename, "fastq")
        records = []
        for item in record_iterator:
            qual = {}
            val = item.letter_annotations['phred_quality'][:-2]
            qual['phred_quality'] = val
            tmp_seq = str(item.seq)[:-2]
            tmp_id = item.id
            tmp_ds = item.description
            tmp_rec = SeqRecord(Seq(tmp_seq), id=tmp_id, description=tmp_ds, letter_annotations=qual)
            records.append(tmp_rec)
        #SeqIO.write(records, output_fasta, "fasta")
        SeqIO.write(records, output_fastq, "fastq")
        return

    def mutliprocess_remove_NN(self, output_dir, input_dir):
        file_names = []
        for fn in os.listdir(input_dir):
            if fn.find("clipped.fq") > 0:
                print(fn)
                file_names.append(fn)
                
        n = len(file_names)
        print("File Number:{}".format(n))
        procs = []
        n_cpu = multiprocessing.cpu_count()
        print("n_cpu:{}".format(n_cpu))
        #chunk_size = int(n/n_cpu)
        #print("chunk_size:{}".format(chunk_size))
        if n_cpu > n:
            print("The CPU/core's number is greater than that of files")
            for filename in file_names:
                print(filename.strip(".fq") + "_NN_removed.fq")
                procs.append(multiprocessing.Process(target=self.remove_NN, args=(input_dir+filename, output_dir+filename.strip(".fq") + "_NN_removed.fq")))
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
        else:
            print("The CPU/core's number is less than that of files")
            times = int(n/n_cpu)
            remain = n % n_cpu
            start = 0
            end = n_cpu
            while times:
                print("Multiprocess {}".format(times))
                print("Start Index:{}".format(start))
                procs1 = []
                print("End Index:{}".format(end))
                for filename1 in file_names[start:end]:
                    procs1.append(multiprocessing.Process(target=self.remove_NN, args=(input_dir+filename, output_dir+filename.strip(".fq") + "_NN_removed.fq")))
                for proc1 in procs1:
                    proc1.start()
                for proc1 in procs1:
                    proc1.join()
                start = end
                end = end + n_cpu
                times = times-1
            if remain:
                print("Multiprocess {}".format(times))
                procs2 = []
                print("Start Index:{}".format(end - n_cpu))
                print("End Index:{}".format(end - n_cpu + remain))
                for filename1 in file_names[end : (end + remain)]:
                    procs2.append(multiprocessing.Process(target=self.remove_NN, args=(input_dir+filename, output_dir+filename.strip(".fq") + "_NN_removed.fq")))
                for proc2 in procs2:
                    proc2.start()
                for proc2 in procs2:
                    proc2.join()
        return

    def count1(self, filename1, filename2, filename3, filename4):
        miRNA_seq = {}
        with open(filename1, 'r') as f1:
            #next(f3)
            for line in f1.readlines()[3:]:
                val = line.split('\t')
                tmp = str(val[1]).replace('U', 'T')
                print(val[0])
                count1 = 0
                count2 = 0
                #for item1, item2 in zip(record_iterator1, record_iterator2):
                '''
                record_iterator1 = SeqIO.parse(filename2, "fastq")
                for item1 in record_iterator1:
                    seq1 = item1.seq
                    if str(seq1) == tmp:
                        #print(seq1)
                        #print(tmp)
                        count1 += 1
                    else:
                        continue
                    #if str(item2.seq) == val[1]:
                    #    count2 += 1
                '''
                with open(filename2, 'r') as f2:
                    for item1 in f2.readlines():
                        #print(item2)
                        #print(item2)
                        tmp_seq1 = item1.strip('\n')
                        if tmp_seq1 == tmp:
                            count1 += 1
                        else:
                            continue
                        #print(tmp)
                with open(filename3, 'r') as f3:
                    for item2 in f3.readlines():
                        #print(item2)
                        #print(item2)
                        tmp_seq = item2.strip('\n')
                        if tmp_seq == tmp:
                            count2 += 1
                        else:
                            continue
                        #print(tmp)
                miRNA_seq[val[0]] = [tmp, count1, count2]
        #js = json.dumps(miRNA_seq)
        with open(filename4, 'w') as f4:
            for mirna in miRNA_seq:
                value = miRNA_seq[mirna]
                f4.write(mirna)
                f4.write('\t')
                f4.write(value[0])
                f4.write('\t')
                f4.write(str(value[1]))
                f4.write('\t')
                f4.write(str(value[2]))
                f4.write('\n')
        return
        
    def file_iter(self, file_name):
        Fn=open(file_name)
        fileiter=(linE.strip() for linE in Fn)
        for Line in fileiter:
            yield(Line)
            
    def counting(self, miRNA_seqs, filename, output_dir):
        '''
        if ".fq" or ".fastq" in filename:
            record_iterator = SeqIO.parse(filename, "fastq")
        elif ".fa" or ".fasta" in filename:
            record_iterator = SeqIO.parse(filename, "fasta")
        seqs = []
        for item in record_iterator:
            seqs.append(str(item.seq))
        '''
        seqs = []
        if ".fa" or "fasta" in filename:
            fa = self.file_iter(filename)
            for i in fa:
                if not ">" in i:
                    seqs.append(i)
                else:
                    continue
        elif ".fq" or ".fastq" in filename:
            fq = self.file_iter(filename)
            for i in fq:
                if not "@" and not "+" and not "E" and not "/" in i:
                    seqs.append(i)
                else:
                    continue
        
        seqs_count = Counter(seqs)
        sqs_keys = seqs_count.keys()
        miRNA_seq_count = {}

        for mirna_id, raw_seq in miRNA_seqs:
            if raw_seq in sqs_keys:
                miRNA_seqs[mirna_id, raw_seq] = seqs_count[raw_seq]
            else:
                miRNA_seqs[mirna_id, raw_seq] = 0
        fout = output_dir + filename.split('/')[-1] + ".count.txt"
        with open(fout, 'w') as fo:
            for mirna, mirna_seq in miRNA_seqs:
                value = miRNA_seqs[mirna, mirna_seq]
                fo.write(mirna)
                fo.write('\t')
                fo.write(mirna_seq)
                fo.write('\t')
                fo.write("{:d}".format(value))
                fo.write('\n')
        return

    def mutliprocess_count(self, miRNA_file, output_dir, input_dir):
        miRNA_seq = {}
        with open(miRNA_file, 'r') as f1:
            for line in f1.readlines()[3:]:
                val = line.split('\t')
                tmp = str(val[1]).replace('U', 'T')
                miRNA_seq[val[0], tmp] = 0

        file_names = []
        for fn in os.listdir(input_dir):
            #if ".fa" or ".fq" or ".fasta" or ".fastq" in fn:
            if fn.find(".fa") > 0:
                print(fn)
                file_names.append(fn)
            elif fn.find(".fq") > 0:
                print(fn)
                file_names.append(fn)
                
        n = len(file_names)
        print("File Number:{}".format(n))
        procs = []
        n_cpu = multiprocessing.cpu_count()
        print("n_cpu:{}".format(n_cpu))
        #chunk_size = int(n/n_cpu)
        #print("chunk_size:{}".format(chunk_size))
        if n_cpu > n:
            print("The CPU/core's number is greater than that of files")
            for filename in file_names:
                procs.append(multiprocessing.Process(target=self.counting, args=(miRNA_seq, input_dir + filename, output_dir)))
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
        else:
            print("The CPU/core's number is less than that of files")
            times = int(n/n_cpu)
            remain = n % n_cpu
            start = 0
            end = n_cpu
            while times:
                print("Multiprocess {}".format(times))
                print("Start Index:{}".format(start))
                procs1 = []
                print("End Index:{}".format(end))
                for filename1 in file_names[start:end]:
                    procs1.append(multiprocessing.Process(target=self.counting, args=(miRNA_seq, input_dir + filename1, output_dir)))
                for proc1 in procs1:
                    proc1.start()
                for proc1 in procs1:
                    proc1.join()
                start = end
                end = end + n_cpu
                times = times-1
            if remain:
                print("Multiprocess {}".format(times))
                procs2 = []
                print("Start Index:{}".format(end - n_cpu))
                print("End Index:{}".format(end - n_cpu + remain))
                for filename1 in file_names[end : (end + remain)]:
                    procs2.append(multiprocessing.Process(target=self.counting, args=(miRNA_seq, input_dir + filename1, output_dir)))
                for proc2 in procs2:
                    proc2.start()
                for proc2 in procs2:
                    proc2.join()
        return

    def mutliprocess_multi_count(self, miRNA_file, output_dir, input_dir1, input_dir2, input_dir3):
        miRNA_seq = {}
        with open(miRNA_file, 'r') as f1:
            for line in f1.readlines()[3:]:
                val = line.split('\t')
                tmp = str(val[1]).replace('U', 'T')
                miRNA_seq[val[0], tmp] = 0

        file_names = []
        for fn in os.listdir(input_dir1):
            #if ".fa" or ".fq" or ".fasta" or ".fastq" in fn:
            if fn.find(".fa") > 0:
                print(fn)
                file_names.append(fn)
            elif fn.find(".fq") > 0:
                print(fn)
                file_names.append(fn)
                
        n = len(file_names)
        print("File Number:{}".format(n))
        procs = []
        n_cpu = multiprocessing.cpu_count()
        print("n_cpu:{}".format(n_cpu))
        #chunk_size = int(n/n_cpu)
        #print("chunk_size:{}".format(chunk_size))
        if n_cpu > n:
            print("The CPU/core's number is greater than that of files")
            for filename in file_names:
                procs.append(multiprocessing.Process(target=self.multi_count, args=(miRNA_seq, input_dir1 + filename, input_dir2 + filename, input_dir3 + filename, output_dir+filename + '.count.txt', output_dir+filename + '.stats.txt', output_dir+filename + '.freqnochangereads.txt')))
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
        else:
            print("The CPU/core's number is less than that of files")
            times = int(n/n_cpu)
            remain = n % n_cpu
            start = 0
            end = n_cpu
            while times:
                print("Multiprocess {}".format(times))
                print("Start Index:{}".format(start))
                procs1 = []
                print("End Index:{}".format(end))
                for filename1 in file_names[start:end]:
                    procs1.append(multiprocessing.Process(target=self.multi_count, args=(miRNA_seq, input_dir1 + filename, input_dir2 + filename, input_dir3 + filename, output_dir+filename + '.count.txt', output_dir+filename + '.stats.txt', output_dir+filename + '.freqnochangereads.txt')))
                for proc1 in procs1:
                    proc1.start()
                for proc1 in procs1:
                    proc1.join()
                start = end
                end = end + n_cpu
                times = times-1
            if remain:
                print("Multiprocess {}".format(times))
                procs2 = []
                print("Start Index:{}".format(end - n_cpu))
                print("End Index:{}".format(end - n_cpu + remain))
                for filename1 in file_names[end : (end + remain)]:
                    procs2.append(multiprocessing.Process(target=self.multi_count, args=(miRNA_seq, input_dir1 + filename, input_dir2 + filename, input_dir3 + filename, output_dir+filename + '.count.txt', output_dir+filename + '.stats.txt', output_dir+filename + '.freqnochangereads.txt')))
                for proc2 in procs2:
                    proc2.start()
                for proc2 in procs2:
                    proc2.join()
        return

    def minmum_edit_distance(self, tar_seq, control_set):
        cur_dis = 100000
        for seq in control_set:
            real_dis = editdistance.eval(tar_seq, seq)
            if cur_dis >= real_dis:
                cur_dis = real_dis
            else:
                continue
        return cur_dis

    def distinct_freqNoChange(self, tar_seq_set, freq_count_dict1, freq_count_dict2, control_set, fout):
        distinct_freq_noChange = {}
        dis0_num = 0
        dis1_num = 0
        dis2_num = 0
        dis3p_num = 0
        dis_num_percent = {}
        for seq in tar_seq_set:
            if freq_count_dict1[seq] == freq_count_dict2[seq]:
                min_edit_dis = self.minmum_edit_distance(seq, control_set)
                distinct_freq_noChange[seq] = [freq_count_dict1[seq], freq_count_dict2[seq], min_edit_dis]
                if min_edit_dis == 0:
                    dis0_num += 1
                elif min_edit_dis == 1:
                    dis1_num += 1
                elif min_edit_dis == 2:
                    dis2_num += 1
                elif min_edit_dis >=3:
                    dis3p_num += 1
        total = dis0_num + dis1_num + dis2_num + dis3p_num
     
        dis_num_percent[0] = [total, dis0_num, (dis0_num*1.0)/(total*1.0)]
        dis_num_percent[1] = [total, dis1_num, (dis1_num*1.0)/(total*1.0)]
        dis_num_percent[2] = [total, dis2_num, (dis2_num*1.0)/(total*1.0)]
        dis_num_percent[3] = [total, dis3p_num, (dis3p_num*1.0)/(total*1.0)]

        with open(fout, 'w') as fo:
            fo.write("Sequence")
            fo.write('\t')
            fo.write("distinct_reads_frequence_noChange_BeforeCorrection")
            fo.write('\t')
            fo.write("distinct_reads_frequence_noChange_afterCorrection")
            fo.write('\t')
            fo.write("min_edit_dis")
            fo.write('\n')
            for seq in distinct_freq_noChange:
                value = distinct_freq_noChange[seq]
                fo.write(str(seq))
                fo.write('\t')
                fo.write(str(value[0]))
                fo.write('\t')
                fo.write(str(value[1]))
                fo.write('\t')
                fo.write(str(value[2]))
                fo.write('\n')

        return dis_num_percent

    def multi_count(self, miRNA_seq, filename1, filename2, filename3, filename4, filename5, filename6):
        record_iterator1 = SeqIO.parse(filename1, "fastq")
        record_iterator2 = SeqIO.parse(filename2, "fastq")
        record_iterator3 = SeqIO.parse(filename3, "fastq")
        record_seqs1 = []
        record_seqs2 = []
        record_seqs3 = []
        errors_sum12 = 0
        errors_sum13 = 0
        for item1, item2, item3 in zip(record_iterator1, record_iterator2, record_iterator3):
            seq1 = item1.seq
            seq2 = item2.seq
            seq3 = item3.seq
            record_seqs1.append(seq1)
            record_seqs2.append(seq2)
            record_seqs3.append(seq3)
            seq12_dis = editdistance.eval(seq1, seq2)
            seq13_dis = editdistance.eval(seq1, seq3)
            errors_sum12 = errors_sum12 + seq12_dis
            errors_sum13 = errors_sum13 + seq13_dis

        seqs_963 = []
        for miRNA, seq in miRNA_seq:
            seqs_963.append(seq)

        seqs_963_set = set(seqs_963)
        record_seqs1_set = set(record_seqs1)
        record_seqs2_set = set(record_seqs2)
        record_seqs3_set = set(record_seqs3)

        noCorrection_963 = record_seqs1_set - seqs_963_set
        miRECCorrection_963 = record_seqs2_set - seqs_963_set
        karectCorrection_963 = record_seqs3_set - seqs_963_set

        noCorrection_963_ssd = record_seqs1_set ^ seqs_963_set
        miRECCorrection_963_ssd = record_seqs2_set ^ seqs_963_set
        karectCorrection_963_ssd = record_seqs3_set ^ seqs_963_set
        
        ###
        noCorrection_miRECCorrection_inter = record_seqs1_set & record_seqs2_set
        noCorrection_karectCorrection_inter = record_seqs1_set & record_seqs3_set

        miRECCorrection_noCorrection_miRECCorrection_inter = record_seqs2_set - noCorrection_miRECCorrection_inter
        karectCorrection_noCorrection_karectCorrection_inter = record_seqs3_set - noCorrection_karectCorrection_inter

        miRECincreased_seq_set = (miRECCorrection_963 - record_seqs1_set)
        karectincreased_seq_set = (karectCorrection_963 - record_seqs1_set)

        miRECnotcorrection_set = noCorrection_963 & miRECCorrection_963
        karectnotcorrection_set = noCorrection_963 & karectCorrection_963

        record_seqs1_dict = Counter(record_seqs1)
        record_seqs2_dict = Counter(record_seqs2)
        record_seqs3_dict = Counter(record_seqs3)

        dis_num_percent_miREC = self.distinct_freqNoChange(noCorrection_miRECCorrection_inter, record_seqs1_dict, record_seqs2_dict, seqs_963_set, filename6)
        dis_num_percent_karect = self.distinct_freqNoChange(noCorrection_karectCorrection_inter, record_seqs1_dict, record_seqs3_dict, seqs_963_set, filename6)

        miRNA_seq_count = {}
        count1_sum = 0
        count2_sum = 0
        count3_sum = 0
        for miRNA, seq in miRNA_seq:
            if seq in record_seqs1_dict:
                count1 = record_seqs1_dict[seq]
                count1_sum = count1_sum + count1
            else:
                count1 = 0
            if seq in record_seqs2_dict:
                count2 = record_seqs2_dict[seq]
                count2_sum = count2_sum + count2
            else:
                count2 = 0
            if seq in record_seqs2_dict:
                count3 = record_seqs3_dict[seq]
                count3_sum = count3_sum + count3
            else:
                count3 = 0
            miRNA_seq_count[miRNA, seq] = [count1, count2, count3]

        with open(filename4, 'w') as fo:
            fo.write("miRNA_ID")
            fo.write('\t')
            fo.write("Sequence")
            fo.write('\t')
            fo.write("CopyNum_BeforeCorrection")
            fo.write('\t')
            fo.write("CopyNum_miREC")
            fo.write('\t')
            fo.write("CopyNum_Karect")
            fo.write('\n')
            for mirna, seq in miRNA_seq_count:
                value = miRNA_seq_count[mirna, seq]
                fo.write(mirna)
                fo.write('\t')
                fo.write(seq)
                fo.write('\t')
                fo.write(str(value[0]))
                fo.write('\t')
                fo.write(str(value[1]))
                fo.write('\t')
                fo.write(str(value[2]))
                fo.write('\n')

        with open(filename5, 'w') as f5:
            print('Total reads count of the dataset is: {}'.format(len(record_seqs1)),file=f5)
            print("miREC:", file=f5)
            print('The number of errors corrected by miREC is: {}'.format(errors_sum12),file=f5)
            print('The number of distinct reads was decreased from {0} to {1} after the correction by miREC'.format(len(record_seqs1_set), len(record_seqs2_set)), file=f5)
            print("The number of distinct sequences in the intersection of the pre-correction and miREC correction sets is: {}".format(len(noCorrection_miRECCorrection_inter)),file=f5)
            print("The number of distinct sequences in the miREC correction set but not in the intersection of the pre-correction and miREC correction sets is: {}".format(len(miRECCorrection_noCorrection_miRECCorrection_inter)), file=f5)
            print("The number of increased new sequences by miREC correction is: {}".format(len(miRECincreased_seq_set)), file=f5)
            #print("The number of sequences that are not corrected by miREC is: {}".format(len(miRECnotcorrection_set)), file=f5)
            print("The read counts of the 963 miRNAs increased about {:.4%} on average (from total {} to {}) after miREC correction".format((count2_sum-count1_sum)/count1_sum, count1_sum, count2_sum), file=f5)            
            print("The number of total distinct reads whose frequences are not changed is：{}".format(dis_num_percent_miREC[0][0]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 0 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_miREC[0][1], dis_num_percent_miREC[0][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 1 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_miREC[1][1], dis_num_percent_miREC[1][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 2 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_miREC[2][1], dis_num_percent_miREC[2][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance >=3 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_miREC[3][1], dis_num_percent_miREC[3][2]), file=f5)
            print("\n", file=f5)
            print("Karect:", file=f5)
            print('The number of errors corrected by Karect is: {}'.format(errors_sum13),file=f5)
            print('The number of distinct reads was decreased from {0} to {1} after the correction by Karect'.format(len(record_seqs1_set), len(record_seqs3_set)), file=f5)
            print("The number of distinct sequences in the intersection of the pre-correction and Karect correction sets is: {}".format(len(noCorrection_karectCorrection_inter)), file=f5)
            print("The number of distinct sequences in the Karect correction set but not in the intersection of the pre-correction and Karect correction sets is: {}".format(len(karectCorrection_noCorrection_karectCorrection_inter)), file=f5)
            print("The number of increased new sequences by Karect correction is: {}".format(len(karectincreased_seq_set)), file=f5)
            #print("The number of sequences that are not corrected by Karect is: {}".format(len(karectnotcorrection_set)), file=f5)
            print("The read counts of the 963 miRNAs increased about {:.4%} on average (from total {} to {})after Karect correction".format((count3_sum-count1_sum)/count1_sum, count1_sum, count3_sum), file=f5)
            print("The number of total distinct reads whose frequences are not changed is: {}".format(dis_num_percent_karect[0][0]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 0 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_karect[0][1], dis_num_percent_karect[0][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 1 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_karect[1][1], dis_num_percent_karect[1][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance 2 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_karect[2][1], dis_num_percent_karect[2][2]), file=f5)
            print("The number and percentage of distinct reads that have a minimum distance >=3 with the set of the 963 miRNAs is: {} and {:.4%}".format(dis_num_percent_karect[3][1], dis_num_percent_karect[3][2]), file=f5)

        return

    def edit_distance(self, fin1, fin2, fout):
        record_iterator1 = SeqIO.parse(fin1, "fastq")
        seqs = []
        for item1 in record_iterator1:
            seqs.append(str(item1.seq))
        #seqs_unique = list(set(seqs))
        
        result = Counter(seqs)
        infrequent_dict = dict((k, v) for k, v in result.items() if v < 5)
        frequent_dict = dict((k, v) for k, v in result.items() if v >= 5)

        infre_control_edit_distance = {}
        control_list = []
        with open(fin2, 'r') as f2:
            for line in f2.readlines()[3:]:
                val = line.split('\t')
                seq22 = str(val[1]).replace('U', 'T')
                control_list.append(seq22)

        control_set = list(set(control_list))
        '''
        control_count = Counter(control_list)
        for k,v in control_count.items():
            if v==2:
                print(k.replace('T', 'U'), v)
        '''
        for seq11 in list(infrequent_dict.keys()):
        #for seq11 in list(infrequent_dict.keys())[:10]:
            cur_dis2 = 1000000
            for tmp_control_seq in control_set:
                real_dis = editdistance.eval(seq11, tmp_control_seq)
                if cur_dis2 >= real_dis:
                    cur_dis2 = real_dis
                    cur_seq2 = tmp_control_seq
                else:
                    continue
            infre_control_edit_distance[seq11, cur_seq2] = cur_dis2    

        infre_control_edit_distance1 = dict((k, v) for k, v in infre_control_edit_distance.items() if v == 1)
        
        frequent_minus_control_dict = dict((k, v) for k, v in frequent_dict.items() if k not in list(set(control_set)))

        selectInfre_fre_edit_distance = {}
        for infre_seq, control_seq in infre_control_edit_distance1:
            cur_dis1 = 1000000
            for fmc_seq in frequent_minus_control_dict:
                real_dis = editdistance.eval(infre_seq, fmc_seq)
                if cur_dis1 >= real_dis:
                    cur_dis1 = real_dis
                    cur_seq1 = fmc_seq
                else:
                    continue
            selectInfre_fre_edit_distance[infre_seq, control_seq, cur_seq1] = [infre_control_edit_distance1[infre_seq, control_seq], cur_dis1]

        with open(fout, 'w') as ff:
            ff.write("infre_seq")
            ff.write('\t')
            ff.write("control_seq")
            ff.write('\t')
            ff.write("fre-control_seq")
            ff.write('\t')
            ff.write("control_edit_distance")
            ff.write('\t')
            ff.write("fre-control_edit_distance")
            ff.write('\n')
            for s1, s2, s3 in selectInfre_fre_edit_distance:
                value1, value2 = selectInfre_fre_edit_distance[s1, s2, s3]
                ff.write(s1)
                ff.write('\t')
                ff.write(s2)
                ff.write('\t')
                ff.write(s3)
                ff.write('\t')
                ff.write(str(value1))
                ff.write('\t')
                ff.write(str(value2))
                ff.write('\n')
        return

    def relevant_seqs(self, fin1, fin2):
        if fin1.find('fq') > 0:
            record_iterator1 = SeqIO.parse(fin1, "fastq")
        elif fin1.find('fa') > 0:
            record_iterator1 = SeqIO.parse(fin1, "fasta")
        seqs = []
        for item1 in record_iterator1:
            seqs.append(str(item1.seq))
        #seqs_unique = list(set(seqs))
        
        result = Counter(seqs)
        control_list = []
        with open(fin2, 'r') as f2:
            for line in f2.readlines()[3:]:
                val = line.split('\t')
                seq22 = str(val[1]).replace('U', 'T')
                control_list.append(seq22)

        control_set = set(control_list) 
        otherthan963_dict = dict((k, v) for k, v in result.items() if k not in control_set)

        other963_infrequent_set = set(k for k, v in otherthan963_dict.items() if v < 5) 

        selection_seqs = other963_infrequent_set | control_set
        return selection_seqs

    def filter(self, miRNA963, fin, max_len, min_len, fout):
        selection_seqs = self.relevant_seqs(fin, miRNA963)
        print('Selection_seqs Number:{}'.format(len(selection_seqs)))
        #print(selection_seqs[0])
        if fin.find('.fq' or 'fastq') > 0:
            record_iterator_fq = SeqIO.parse(fin, "fastq")
            records_fq = []
            #num = 0
            for item in record_iterator_fq:
                qual = {}
                qual['phred_quality'] = item.letter_annotations['phred_quality']
                tmp_seq = str(item.seq)
                tmp_id = item.id
                tmp_len = len(tmp_seq)
                #if (min_len - 2) < tmp_len and tmp_len <= (max_len + 2) and tmp_seq in selection_seqs:
                
                if tmp_seq in selection_seqs:
                    if (min_len - 2) < tmp_len and tmp_len <= (max_len + 2):
                        #print(num)
                        tmp_rec = SeqRecord(Seq(tmp_seq), id=tmp_id, description=tmp_id, letter_annotations=qual)
                        records_fq.append(tmp_rec)
                        #num += 1
                    else:
                        continue
                else:
                    continue
            SeqIO.write(records_fq, fout, "fastq")
        elif fin.find('.fa' or 'fasta') > 0: 
            record_iterator_fa = SeqIO.parse(fin, "fasta")
            records_fa = []
            num = 0
            for item in record_iterator_fa:
                tmp_seq = str(item.seq)
                tmp_id = item.id
                tmp_len = len(tmp_seq)
                if tmp_seq in selection_seqs:
                    if (min_len - 2) < tmp_len and tmp_len <= (max_len + 2):
                        print(num)
                        tmp_rec = SeqRecord(Seq(tmp_seq), id=tmp_id, description=tmp_id)
                        records_fa.append(tmp_rec)
                        num += 1
                    else:
                        continue
                else:
                    continue
            SeqIO.write(records_fa, fout, "fasta")
        return

    def mutliprocess_filter(self, miRNA_file, input_dir, output_dir):
        miRNA_seq = {}
        with open(miRNA_file, 'r') as f1:
            #next(f3)
            for line in f1.readlines()[3:]:
                val = line.split('\t')
                tmp = str(val[1]).replace('U', 'T')
                print(val[0])
                miRNA_seq[val[0], tmp] = len(tmp)
        values = miRNA_seq.values()
        max_len = max(values)
        min_len = min(values)

        file_names = []
    
        for fn in os.listdir(input_dir):
            #if ".fa" or ".fq" or ".fasta" or ".fastq" in fn:
            if fn.find(".fq") > 0:
                print(fn)
                file_names.append(fn)
            elif fn.find(".fa") > 0:
                print(fn)
                file_names.append(fn)
            else:
                continue

        n = len(file_names)
        print("File Number:{}".format(n))
        procs = []
        n_cpu = multiprocessing.cpu_count()
        print("n_cpu:{}".format(n_cpu))
        #chunk_size = int(n/n_cpu)
        #print("chunk_size:{}".format(chunk_size))
        if n_cpu > n:
            print("The CPU/core's number is greater than that of files")
            for filename in file_names:
                procs.append(multiprocessing.Process(target=self.filter, args= (miRNA_file,input_dir + filename, max_len, min_len, output_dir + filename)))
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
        else:
            print("The CPU/core's number is less than that of files")
            times = int(n/n_cpu)
            remain = n % n_cpu
            start = 0
            end = n_cpu
            while times:
                print("Multiprocess {}".format(times))
                print("Start Index:{}".format(start))
                procs1 = []
                print("End Index:{}".format(end))
                for filename1 in file_names[start:end]:
                    procs1.append(multiprocessing.Process(target=self.filter, args= (miRNA_file, input_dir + filename1, max_len, min_len, output_dir + filename1)))
                for proc1 in procs1:
                    proc1.start()
                for proc1 in procs1:
                    proc1.join()
                start = end
                end = end + n_cpu
                times = times-1
            if remain:
                print("Multiprocess {}".format(times))
                procs2 = []
                print("Start Index:{}".format(end - n_cpu))
                print("End Index:{}".format(end - n_cpu + remain))
                for filename1 in file_names[end : (end + remain)]:
                    procs2.append(multiprocessing.Process(target=self.filter, args= (miRNA_file, input_dir + filename1, max_len, min_len, output_dir + filename1)))
                for proc2 in procs2:
                    proc2.start()
                for proc2 in procs2:
                    proc2.join()
        return
    
    def multiprocess_forward_reverse_inter(self, input_dir, output_dir):
        filenames = []
        for fn in os.listdir(input_dir):
            #print(fn)
            if fn.find(".fastq") > 0:
                print(fn)
                filenames.append(fn)
            else:
                continue
        print(len(filenames))
        group_filenames = []
        num = 0
        while len(filenames):
            print(len(filenames))
            file1 = filenames[0]
            #print(file1)
            file1_split = file1.split('_')
            filenames.remove(file1)
            for item in filenames:
                #print(item)
                if file1_split[1] in item:
                    if file1_split[2] == '2':
                        group_filenames.append([item, file1])
                    elif file1_split[2] == '1':
                        group_filenames.append([file1, item])
                    print(file1, item)
                    filenames.remove(item)
                    num = num + 1
                    break
                else:
                    continue 
        
        n = len(group_filenames)
        print("File Number:{}".format(n))
        procs = []
        n_cpu = multiprocessing.cpu_count()
        print("n_cpu:{}".format(n_cpu))
        #chunk_size = int(n/n_cpu)
        #print("chunk_size:{}".format(chunk_size))
        if n_cpu > n:
            print("The CPU/core's number is greater than that of files")
            for item in group_filenames:
                f = item[0]
                r = item[1]
                print(input_dir+f, input_dir+r)
                procs.append(multiprocessing.Process(target=self.forward_reverse_inter, args= (input_dir+f, input_dir+r, output_dir+f.split('.fastq')[0] +'_intersection.fastq', output_dir+r.split('.fastq')[0] +'_intersection.fastq')))
            for proc in procs:
                proc.start()
            for proc in procs:
                proc.join()
        else:
            print("The CPU/core's number is less than that of files")
            times = int(n/n_cpu)
            remain = n % n_cpu
            start = 0
            end = n_cpu
            while times:
                print("Multiprocess {}".format(times))
                print("Start Index:{}".format(start))
                procs1 = []
                print("End Index:{}".format(end))
                for item in group_filenames[start:end]:
                    f = item[0]
                    r = item[1]
                    procs1.append(multiprocessing.Process(target=self.forward_reverse_inter, args= (input_dir+f, input_dir+r, output_dir+f.split('.fastq')[0] +'_intersection.fastq', output_dir+r.split('.fastq')[0] +'_intersection.fastq')))
                for proc1 in procs1:
                    proc1.start()
                for proc1 in procs1:
                    proc1.join()
                start = end
                end = end + n_cpu
                times = times-1
            if remain:
                print("Multiprocess {}".format(times))
                procs2 = []
                print("Start Index:{}".format(end - n_cpu))
                print("End Index:{}".format(end - n_cpu + remain))
                for item in group_filenames[end : (end + remain)]:
                    f = item[0]
                    r = item[1]
                    procs2.append(multiprocessing.Process(target=self.forward_reverse_inter, args= (input_dir+f, input_dir+r, output_dir+f.split('.fastq')[0] +'_intersection.fastq', output_dir+r.split('.fastq')[0] +'_intersection.fastq')))
                for proc2 in procs2:
                    proc2.start()
                for proc2 in procs2:
                    proc2.join()
        return

    def forward_reverse_inter(self, forward_input, reverse_input, forward_output, reverse_output):
        #print(forward_input)
        #print(forward_input.find('.fq' or '.fastq'))
        if '.fq' or '.fastq' in forward_input:
            record_iterator_forward = SeqIO.parse(forward_input, "fastq")
        else:
            print("Forward {} may not be a fatsq file.".format(forward_input))

        if '.fq' or '.fastq' in reverse_input:
            record_iterator_reverse = SeqIO.parse(reverse_input, "fastq")
        else:
            print("Reveser {} may not be a fatsq file.".format(reverse_input))

        forward_id = []
        reverse_id = []
        for forward_item, reverse_item in zip(record_iterator_forward, record_iterator_reverse):
            forward_tmp_id = (forward_item.id).split('#')[0]
            reverse_tmp_id = (reverse_item.id).split('#')[0]
            forward_id.append(forward_tmp_id)
            reverse_id.append(reverse_tmp_id)
        
        forward_reverse_id_inter = set(forward_id) & set(reverse_id)

        if '.fq' or '.fastq' in forward_input:
            record_iterator_forward = SeqIO.parse(forward_input, "fastq")
        else:
            print("Forward {} may not be a fatsq file.".format(forward_input))

        if '.fq' or '.fastq' in reverse_input:
            record_iterator_reverse = SeqIO.parse(reverse_input, "fastq")
        else:
            print("Reveser {} may not be a fatsq file.".format(reverse_input))
        count1 = 0
        count2 = 0
        forward_records = []
        reverse_records = []
        for forward_item, reverse_item in zip(record_iterator_forward, record_iterator_reverse):
            qual_forward = {}
            qual_forward['phred_quality'] = forward_item.letter_annotations['phred_quality']
            forward_tmp_seq = str(forward_item.seq)
            forward_tmp_id = forward_item.id

            qual_reverse = {}
            qual_reverse['phred_quality'] = reverse_item.letter_annotations['phred_quality']
            reverse_tmp_seq = str(reverse_item.seq)
            reverse_tmp_id = reverse_item.id

            forward_tmp_id_part = forward_tmp_id.split('#')[0]
            #print(forward_tmp_id_part)
            reverse_tmp_id_part = reverse_tmp_id.split('#')[0]
            if forward_tmp_id_part in forward_reverse_id_inter:
                #print(count1)
                forward_tmp_rec = SeqRecord(Seq(forward_tmp_seq), id=forward_tmp_id, description=forward_tmp_id, letter_annotations=qual_forward)
                forward_records.append(forward_tmp_rec)
                count1 = count1+1
            if reverse_tmp_id_part in forward_reverse_id_inter:
                #print(count2)
                reverse_tmp_rec = SeqRecord(Seq(reverse_tmp_seq), id=reverse_tmp_id, description=reverse_tmp_id, letter_annotations=qual_reverse)
                reverse_records.append(reverse_tmp_rec)
                count2 = count2+1
        SeqIO.write(forward_records, forward_output, "fastq")
        SeqIO.write(reverse_records, reverse_output, "fastq")
        return
