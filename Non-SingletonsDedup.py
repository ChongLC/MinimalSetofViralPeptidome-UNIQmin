from Bio import SeqIO
import pandas as pd

fileA = list(SeqIO.parse("cdhitParamyxoviridae","fasta"))

#frequency count
kmers = pd.read_csv("Output_kmers_Paramyxoviridae.txt", header=None)
kmers.columns = ['kmer']
kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

#nonSingleOccuranceNonamer
more1 = kmers['freq']!=1
more1.head()
kmer_more1 = kmers[more1]
more1List = kmer_more1['kmer']
more1List.to_csv("seqmore1List.txt", index = False, header = False)

lines_seen = set()
outfile = open("nr_more1List.txt","w")
for line in open("seqmore1List.txt","r"):
	if line not in lines_seen:
		outfile.write(line)
		lines_seen.add(line)
outfile.close()

#Remove/delete sequences by ID from multifasta
header_set = set(line.strip() for line in open("seqfileZ.txt"))
remainingSeq = open("remainingSeq.fasta","w")

for seq_record in fileA:
	try:
		header_set.remove(seq_record.name)
	except KeyError:
		remainingSeq.write(seq_record.format("fasta"))
remainingSeq.close()

#Extract A Group Of Fasta Sequences From A File
fasta_file = "cdhitParamyxoviridae" # Input fasta file
wanted_file = "seqfileZ.txt" # Input interesting sequence IDs, one per line
result_file = "result_file.fasta" # Output fasta file

wanted = set()
with open(wanted_file) as f:
	for line in f:
		line = line.strip()
		if line != "":
			wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
	for seq in fasta_sequences:
		if seq.id in wanted:
			SeqIO.write([seq], f, "fasta")
