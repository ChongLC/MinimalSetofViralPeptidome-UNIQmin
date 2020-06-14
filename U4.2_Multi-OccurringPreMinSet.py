from Bio import SeqIO
import ahocorasick
import logging 

def load_data(fasta_file, kmer_file):
	logging.info("Loading fasta file")
	fasta_list = list(SeqIO.parse(fasta_file,"fasta"))
	logging.info("Loading kmer list")
	kmer_list = [line.rstrip('\n') for line in open (kmer_file)]
	return fasta_list, kmer_list

def find_match(line, A):
	found_kmers = []
	for end_index, kmer in A.iter(line):
		found_kmers.append(kmer)
	return found_kmers

def setup_automaton(kmer_list):
	logging.info("Setting up kmer lookup")
	auto = ahocorasick.Automaton()
	for seq in kmer_list:
		auto.add_word(seq, seq)
	auto.make_automaton()
	logging.info("Completed set-up of kmer lookup")
	return auto 

def match_kmers(fasta_list, kmer_auto):
	logging.info("Writing output")
	with open(output_file, "w") as f:
		for record in fasta_list:
			match = find_match(str(record.seq), kmer_auto)
			if match:
				line = str(match) + "\n"
				f.write(line)
	logging.info("Completed")

if __name__ == '__main__':
	fasta_file = "result_file.fasta"
	kmer_file = "nr_more1List.txt"
	output_file = "matchKmer4CleanKmer.txt"

	logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
	
	fasta_list, kmer_list = load_data(fasta_file, kmer_file)
	kmer_auto = setup_automaton(kmer_list)
	match_kmers(fasta_list, kmer_auto)
