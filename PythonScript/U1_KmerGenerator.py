from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import math

fileA = list(SeqIO.parse("inputfile.fas","fasta"))
file_id = "Output_kmers.txt"

def generate_kmers(start, end):
	for record in fileA[start:end]:
		nr_sequence = record.seq
		seq_len = len(nr_sequence)
		kmer = 9
		count = 0
		temp = []
		for seq in list(range(seq_len-(kmer-1))):
			count += 1
			my_kmer = (nr_sequence[seq:seq+kmer])
			temp.append(str(my_kmer))
		with open(file_id, 'a') as f:
			f.writelines("%s\n" % kmer for kmer in temp)

if __name__ == '__main__':
  n = len(fileA)
  pool = ProcessPoolExecutor(14)
  futures = []
  perCPUSize = math.ceil(n/14)
  for i in range(0,14):
  	futures.append(pool.submit(generate_kmers, i * perCPUSize, (i+1) * perCPUSize))
