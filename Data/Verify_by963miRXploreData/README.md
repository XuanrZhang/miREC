# Download the dataset
## 1. The datasets can be downloaded from Google Drive (https://drive.google.com/drive/folders/1YDPxrH_B-StPKLgYnkcDXhtCnq5vatt3?usp=sharing). 

For reproducing the test process without any modification, you can put the datasets into the relevant folders(like the following tree)  
```
│   └── Verify_by963miRXploreData
│       ├── 963miRNAs
│       │   └── GSE139936_180719_GEO_miRNAs.txt
│       ├── corrected
│       │   ├── D18_Karect
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
│       ├── D18_filtered_NN_removed_raw_fq
│       │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       ├── D18_NN_removed_raw_fq
│       │   ├── 180719Ded_D18-6962_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_1_sequence.3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped_NN_removed.fq
│       │   ├── 180719Ded_D18-6964_1_sequence.3clipped_NN_removed.fq
│       │   └── 180719Ded_D18-6964_2_sequence.rc3clipped_NN_removed.fq
│       ├── D18_raw_fq
│       │   ├── 180719Ded_D18-6962_1_sequence.3clipped.fq
│       │   ├── 180719Ded_D18-6962_2_sequence.rc3clipped.fq
│       │   ├── 180719Ded_D18-6963_1_sequence.3clipped.fq
│       │   ├── 180719Ded_D18-6963_2_sequence.rc3clipped.fq
│       │   ├── 180719Ded_D18-6964_1_sequence.3clipped.fq
│       │   └── 180719Ded_D18-6964_2_sequence.rc3clipped.fq
│       ├── D19_NN_removed_raw_fq
│       │   └── D19-10246.assembled.2NN.fastq

```

### 2. You also can download the datasets from NCBI
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