# miREC
An error correction tool for miRNA reads at single-base resolution, which proposed a novel 3-layer lattice structure combining kmers, (k-1)mers and (k+1)mers to first time solve the problem of correcting erroneous bases in miRNA sequencing data. 
The novelty of our method is the use of a 3-layer (k-1)mer-kmer-(k+1)mer lattice structure to maintain the frequency differences of the kmers.
The method is particularly effective for the accurate correction of indel errors. 

Our miREC is a user-friendly tool and suits for different research needs. It provides sev- eral parameters for users to specify according to specific tasks. Three most useful settings are: the error types, the frequency threshold τ, and the kmer range [k_1,k_end]. Our miREC has two running modes. One is for the substitution error correction only, the other is Unlike the exist- ing methods which use a fixed threshold, miREC provides options to obtain a good threshold to identify erroneous kmers, thereby avoiding over-correction or insufficient performance when datasets change.  

## Dependancies
KMC3 tool, a kmer counter, is used to obtain kmer frequences. Here is the instruction of KMC3 (http://sun.aei.polsl.pl/REFRESH/kmc).

## Download & Usage

	git clone https://github.com/XuanrZhang/miREC
	cd InsEC	
	chmod +x run.sh
	pip install pandas (optional command, if you don't install pandas before)
	
	./run.sh [sequence.fa] [reads.fa]
	e.g 
	./run.sh ref.fa reads.fa
  
## Data format
Input: A read dataset in .fasta\fastq format and an nucleotide sequence of interests in .fasta format

	- reads.fa(\fq) : store the whole read data needed to be corrected.
	- ref.fa : store the nucleotide sequence of interests (e.g. a sequence of genes or a sequence of coding region).

Output: A corrected read dataset only related to the given nucleotide sequence and An assembled contig based on corrected data

	- corrected.fa :store corrected read data.
	- contigs.fa :store the updated nucleotide sequence of the gene (or the genome region of interests).

	
## Citation
Please cite the work "Aberration-corrected ultrafine analysis of miRNA reads at single-base resolution: better data makes better conclusion."

## Citation
If any bugs during your running, please email to xuan.zhang-5@student.uts.edu.au
