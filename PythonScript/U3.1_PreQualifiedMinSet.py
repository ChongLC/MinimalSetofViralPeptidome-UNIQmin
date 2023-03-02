import argparse
import logging
import ahocorasick as ahc
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(prog = 'U3.1_PreQualifiedMinSet.py',
                    description = 'Step 3.1 of UNIQmin'
                    )
    parser.add_argument('-i', '--input', dest="input", help='Path of the input file (in FASTA format)')
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')

    return parser.parse_args()


class PreQualifiedMinSet:

    def load_data(self, fasta_file, kmer_file):
        logging.info("Loading fasta file to determine pre-qualified minimal set")
        fasta_list = list(SeqIO.parse(fasta_file, "fasta"))
        logging.info("Loading kmer list")
        kmer_list = [line.rstrip('\n') for line in open(kmer_file)]
        return fasta_list, kmer_list

    def __find_match(self, line, A):
        found_kmers = []
        for end_index, kmer in A.iter(line):
            found_kmers.append(kmer)
        return found_kmers

    def setup_automaton(self, kmer_list):
        logging.info("Setting up kmer lookup")
        auto = ahc.Automaton()
        for seq in kmer_list:
            auto.add_word(seq, seq)
        auto.make_automaton()
        logging.info("Completed set-up of kmer lookup")
        return auto

    def match_kmers(self, output_file, fasta_list, kmer_auto):
        logging.info("Writing output (pre-qualified set)")
        with open(output_file, "w") as f:
            for record in fasta_list:
                match = self.__find_match(str(record.seq), kmer_auto)
                if match:
                    line = record.id + "\n"
                    f.write(line)
        logging.info("Completed determining pre-qualified minimal set")


if __name__ == '__main__':
    args = get_args()

    fasta_file = args.input
    kmer_file = args.output + "/seqSingleList.txt"
    output_file = args.output + "/seqfileZ.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    preQualified = PreQualifiedMinSet()

    fasta_list, kmer_list = preQualified.load_data(fasta_file, kmer_file)
    kmer_auto = preQualified.setup_automaton(kmer_list)
    preQualified.match_kmers(output_file, fasta_list, kmer_auto)