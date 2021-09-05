from Bio import SeqIO

fileA = list(SeqIO.parse("cdhitCoronaviridae","fasta"))
header_set = set(line.strip() for line in open("seqfileZ.txt"))
remainingSeq = open("remainingSeq.fasta","w")

for seq_record in fileA:
	try:
		header_set.remove(seq_record.name)
	except KeyError:
		remainingSeq.write(seq_record.format("fasta"))
remainingSeq.close()

fasta_file = "cdhitCoronaviridae" 
wanted_file = "seqfileZ.txt" 
result_file = "result_file.fasta" 

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
