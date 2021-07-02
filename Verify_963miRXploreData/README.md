# Verification results on the sequencing reads of the 963 miRXplore Universal Reference miRNAs (pure control and spike-in)

The code in the **Verify_963miRXploreData** directory is used to test **miREC** in comparison with **Karect**(https://github.com/aminallam/karect) on sequencing reads datasets of **963 miRXplore Universal Reference miRNAs** (three replicates) and their spike-in at eukaryotic cells.

The test was to **verify**：
- whether our detected erroneous reads can be all corrected into one of the 963 miRNAs, and
- whether any new sequences are generated after the correction step.

An ideal performance should be: the corrected reads are all the sequencing reads of the 963 miRNAs, and previously non-existing reads are never not created.

# Dependencies
## To support **miREC**
Please follow the instructions of README.md in the miREC directory.
## To support this analysis
 - biopython
 - editdistance
 - pandas
 - python 3.9
 - multiprocessing
 - collections
 - os
 - subprocess

# Download & Usage
## **Step 1**
```
git clone https://github.com/XuanrZhang/miREC
cd miREC
chmod +x ./Verify_963miRXploreData/miRECD19.sh
chmod +x ./Verify_963miRXploreData/miRECD18.sh
chmod +x ./Verify_963miRXploreData/karectD19.sh
chmod +x ./Verify_963miRXploreData/karectD18.sh
chmod +x ./Verify_963miRXploreData/karect
chmod +x ./Verify_963miRXploreData/mkdir.sh
```
## **Step 2**
Download the dataset (https://drive.google.com/drive/folders/1YDPxrH_B-StPKLgYnkcDXhtCnq5vatt3?usp=sharing) and put them into the relevant folders(like the following tree)  
```
├── Data
│   ├── 963miRXplore_data
│   │   ├── 963miRNAs
│   │   │   └── GSE139936_180719_GEO_miRNAs.txt
│   │   ├── D18_raw_fq
│   │   │   ├── 180719Ded_D18-6962_1_sequence.3clipped.fq
│   │   │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped.fq
│   │   │   ├── 180719Ded_D18-6963_1_sequence.3clipped.fq
│   │   │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped.fq
│   │   │   ├── 180719Ded_D18-6964_1_sequence.3clipped.fq
│   │   │   └── 180719Ded_D18-6964_2_sequence.rc3clipped.fq
│   │   ├── D19_NN_removed_raw_fq
│   │   │   └── D19-10246.assembled.2NN.fastq
```
## **Step 2**
- *Notes*: please make sure you have installed the dependencies in your environment.
### **For the dataset D19-10246**
```
python ./Verify_963miRXploreData/runD19.py 
```
- For runD19.py, it runs with the parameters of kmers [8, 25] on the datasets included in the folders of ./Data/963miRXplore_data/D19_NN_removed_raw_fq/, ./Data/963miRXplore_data/corrected/D19_karect/, ./Data/963miRXplore_data/corrected/D19_miREC/ and ./Data/963miRXplore_data/963miRNAs/.
### **For the datasets of D18-6962_1, D18-6962_2, D18-6963_1, D18-6963_2, D18-6964_1 and D18-6964_2**
```
python ./Verify_963miRXploreData/runD18.py 
```
- For runD18.py, it not only runs the same test with the same parameters like runD19.py on the datasets included in the folders of ./Data/963miRXplore_data/D18_filtered_NN_removed_raw_fq/, ./Data/963miRXplore_data/corrected/D18_karect/, ./Data/963miRXplore_data/corrected/D18_miREC/ and ./Data/963miRXplore_data/963miRNAs/, but also runs some other steps for preprocessing the datasets D18. These preprocessing steps are unnecessary for the dataset D19-10246.
- runD18.py runs about an hour to be finished on a computational node with >=6 cores for 6 .fastq files.

## You also can download the datasets from NCBI
1. GEO accession GSE139936.GSM4149813 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400510)
- 180719Ded_D18-6962_1_sequence.fastq 
- 180719Ded_D18-6962_2_sequence.fastq 
2. GEO accession GSE139936.GSM4149814 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400511)
- 180719Ded_D18-6963_1_sequence.fastq
- 180719Ded_D18-6963_2_sequence.fastq
3. GEO accession GSE139936.GSM4149815 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400512)
- 180719Ded_D18-6964_1_sequence.fastq
- 180719Ded_D18-6964_2_sequence.fastq
4. GEO accession GSE159434 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159434)
- D1910246.assembled.2NN.fastq 

# Output
The output will be saved to the directory of ./Data/963miRXplore_data/output/. For each dataset (e.g. D1910246.assembled.2NN.fastq), it will produce three output (e.g. D1910246.assembled.2NN.fastq.count.txt, D1910246.assembled.2NN.fastq.frenochangereads.txt and D1910246.assembled.2NN.fastq.stats.txt)
- *.count.txt saves the copy numbers before and after correction by **miREC** and **Karect**.
- *.frenochangereads.txt records the sequences whose frequencies are unchanged between before and after correction and the editing distance with the 963 miRNAs.
- *.stats.txt has detailed descriptions for itself.
# Question
If you run with any bugs, please email to pengyao.ping@student.uts.edu.au
