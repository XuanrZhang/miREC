# miREC
An error correction tool for miRNA reads at single-base resolution, which proposed a novel 3-layer lattice structure combining kmers, (k-1)mers and (k+1)mers to first time solve the problem of correcting erroneous bases in miRNA sequencing data. 

# src
The code in the src directory is used to test miREC in comparison with Karect on sequencing reads datasets of 963 miRXplore Universal Reference miRNAs (three replicates) and their spike-in at eukaryotic cells.
The test was to verify
- whether our detected erroneous reads can be all corrected into one of the 963 miRNAs, and
- whether any new sequences are generated after the correction step.

An ideal performance should be: the corrected reads are all the sequencing reads of the 963 miRNAs, and previously non-existing reads are never not created.

# Dependencies
## To support miREC
Please follow the instructions of README.md in the root directory
## To support this analysis
 - biopython
 - editdistance
 - pandas
 - python 3
 - multiprocessing
 - collections
 - os
 - subprocess

# Download & Usage
## 1. If you do not want to do error correction, only to perform and get the performance comparison result
- Step1 
```
git clone https://github.com/Jappy0/miREC
``` 
- Step2 Download the dataset (https://drive.google.com/drive/folders/1YDPxrH_B-StPKLgYnkcDXhtCnq5vatt3?usp=sharing) and put them into the relevant folders(like the following tree)  
```
│   └── synthetic_963miRNAs_Reads
│       ├── 963miRNAs
│       │   └── GSE139936_180719_GEO_miRNAs.txt
│       ├── corrected
│       │   ├── D18_karect
│       │   │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── D18_miREC_8_20
│       │   │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── D18_miREC_8_25
│       │   │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── D19_Karect
│       │   │   └── D19-10246.assembled.2NN.fastq
│       │   ├── D19_miREC_8_20
│       │   │   └── D19-10246.assembled.2NN.fastq
│       │   └── D19_miREC_8_25
│       │       └── D19-10246.assembled.2NN.fastq
│       ├── D19_NN_removed_raw_fq
│       │   └── D19-10246.assembled.2NN.fastq
│       ├── filtered_NN_removed_raw_fq
│       │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       ├── NN_removed_raw_fq
│       │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
```
- Step 3 
```
python ./src/runD18825.py 

```
For runD18825.py, it uses parameters of kmers [8, 25] and the datasets included in the folders of filtered_NN_removed_raw_fq, corrected/D18_miREC_8_25, corrected/D18_karect and 963miRNAs

Notes: please make sure you have installed the dependencies in your environment.

### You also can download the datasets from NCBI
1. GEO accession GSE139936.GSM4149813 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400510)
- 180719Ded_D18-6962_1_sequence.fastq 
- 180719Ded_D18-6962_2_sequence.fastq 
2. GEO accession GSE139936.GSM4149814 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400511)
- 180719Ded_D186963_1_sequence.fastq
- 180719Ded_D18-6963_2_sequence.fastq
3. GEO accession GSE139936.GSM4149815 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10400512)
- 180719Ded_D18-6964_1_sequence.fastq
- 180719Ded_D18-6964_2_sequence.fastq
4. GEO accession GSE159434 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159434)
- D1910246.assembled.2NN.fastq 

## 2. If you want to do error correction by yourself
### 2.1 Do correction using miREC 
- First, please follow the instructions of README.md in the root directory to make sure miREC.sh works well
- Then, run the following script file.
```
./src/runD18825.sh 
```
### 2.2 Do error correction using Karect
Please see the README.md of https://github.com/aminallam/karect, after correction, you should move or copy the corrected results to the folder like ./Data/synthetic_963miRNAs_Reads/corrected/D18_Karect/. And you should remove the file names prefix of "Karect_". 
Notes: For using the files conveniently, I used the same file names to save the sequences after correction.
# Output
The output will be saved to the directory of ./Data/synthetic_963miRNAs_Reads/output/ . For each dataset (e.g. D1910246.assembled.2NN.fastq), it will produce three output (e.g. D1910246.assembled.2NN.fastq.count.txt, D1910246.assembled.2NN.fastq.frenochangereads.txt and D1910246.assembled.2NN.fastq.stats.txt)
- *.count.txt saves the copy numbers before and after correction by miREC and Karect.
- *.frenochangereads.txt records the sequences whose frequencies are unchanged between before and after correction and the editing distance with the 963 miRNAs.
- *.stats.txt has detailed descriptions for itself.
