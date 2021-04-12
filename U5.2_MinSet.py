import sys
from Bio import SeqIO

fasta_file =  sys.argv[1] # Input fasta file
wanted_file = "fileZ.txt" # Input interesting sequence IDs, one per line
result_file = "FileZ.fasta" # Output fasta file

wanted = set()
with open (wanted_file) as f: 
  for line in f: 
    line = line.strip()
    if line != "":
      wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open (result_file, "w") as f: 
  for seq in fasta_sequences: 
    if seq.id in wanted: 
      SeqIO.write([seq], f, "fasta")
