
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

f = "./data/input/180719Ded_D18-6962_1_sequence.10bp.fa"
record_iterator1 = SeqIO.parse(f, "fasta")
for item1 in record_iterator1:
    seq1 = item1.seq
    seq_id = item1.id
    print("test")

    #if str(item2.seq) == val[1]:
    #    count2 += 1