# Simulated data generation

This simulated data generation tool has been implemented as a software prototype. It provides induced-error profile, groud-truth datasets and  error-injected datasets.

Aim to generate datasets with a close nature to wet-lab miRNA sequencing reads, we have two considerations in the process. One is to computationally replicate lab-verified miRNA sequences as templates to form the basic sequences of the simulated datasets, then we duplicate these basic sequences such that the copy counts of them follow a real distribution from a wet-lab dataset of miRNA sequencing reads. In fact, we replicated the mature miRNA sequences in miRBase as the templates, and made the copy count distribution of these template sequences to follow the distribution drawn from a typical miRNA dataset.

- The miRNA sequence file from miRBase is uploaed and named "mature.fa".

- A real distribution file from a wet-lab dataset (Accession ID:SRR866573) is uploaded and named "distrubution.txt"

For users' convenience, users can use their own sequence templates and copy number distribution as well.

# Make and Usage

cd miREC/Generate_SimulatedData/

make

chmod +x gene_simu.sh


## example

Input files: 

- "-d" : distrubution file: ID and sequence (default "distrubution.txt")
- "-f" : sequence templates: fasta file (default "mature.fa")

Output files:

-  "-o" : simulated dataset, reads with errors (default "simulated.fa")
-  "-g" : the groud-truth of simulated dataset  (default "truth.fa")
- error profiles: (default "err.txt")

Optional parameter:

- "-t" : generate datasets with subs errors only (don't use -t, generate datasets with mixed errors)
- "-s" : random seed number (setting different seed number to obtain different datasets)

**Usage**

./gene_simu.sh -s 1 -d ./distrubution.txt -f ./mature.fa -t (generate datasets with subs errors only)

./gene_simu.sh -s 1 -d ./distrubution.txt -f ./mature.fa (generate datasets with mixed errors only)



