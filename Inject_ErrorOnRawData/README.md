# Inject known errors into wet-lab sequence datasets

To rigorously evaluate the error correction performance, we also randomly and purposely inject a small number of errors into these wet-lab datasets, rather than the simulated datasets, to see whether our algorithm can detect and correct these errors with ground truth, together with other errors without ground truth. Only when all of these artificial errors in the real- life miRNA sequencing reads can be detected and corrected, the corrections on the other bases (without ground truth) can be highly trustable. This small number of artificial errors constitutes only 0.5% of total corrections in each dataset to avoid changing the original nature of the data. 

Users can provide their own sequence datasets (fastq files) to randomly injected some errors and record its profile for futher evaluation.

# Make and Usage

cd miREC/Inject_ErrorOnRawData/

make

chmod +x inject_error.sh


## example

Input files: 

- "-f" : sequence dataset provided (default "distrubution.txt")

Output files:

-  "-o" : simulated dataset, reads with errors (default "simulated.fa")
-  "-g" : the groud-truth of simulated dataset  (default "truth.fa")
- error profiles: (default "err.txt")
