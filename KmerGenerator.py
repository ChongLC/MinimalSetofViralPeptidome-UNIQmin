from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import math

fileA = list(SeqIO.parse("cdhitParamyxoviridae","fasta"))
file_id = "Output_kmers_Paramyxoviridae.txt"

def generate_kmers(start, end):
	for record in fileA[start:end]:
		#file_id = "Fflaviviridae/Output_kmers_All-" + str(start) + "-" +  str(end) + ".txt"
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
  pool = ProcessPoolExecutor(12)
  futures = []
  perCPUSize = math.ceil(n/12)
  for i in range(0,12):
  	futures.append(pool.submit(generate_kmers, i * perCPUSize, (i+1) * perCPUSize))
