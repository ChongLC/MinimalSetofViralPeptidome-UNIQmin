import argparse
import logging
import ahocorasick as ahc
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(prog = 'U4.2_Multi-OccurringPreMinSet.py',
                    description = 'Step 4.2 of UNIQmin'
                    )
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')
    
    return parser.parse_args()


class MultiOccurringPreMinSet:

    def load_data_multi(self, fasta_file, kmer_file):
        logging.info("Loading fasta file for minimal set from multi-occuring kmer list")
        fasta_list = list(SeqIO.parse(fasta_file, "fasta"))
        logging.info("Loading multi-occuring kmer list")
        kmer_list = [line.rstrip('\n') for line in open(kmer_file)]
        return fasta_list, kmer_list

    def __find_match_multi(self, line, A):
        found_kmers = []
        for end_index, kmer in A.iter(line):
            found_kmers.append(kmer)
        return found_kmers

    def setup_automaton_multi(self, kmer_list):
        logging.info("Setting up multi-occuring kmer lookup")
        auto = ahc.Automaton()
        for seq in kmer_list:
            auto.add_word(seq, seq)
        auto.make_automaton()
        logging.info("Completed set-up of multi-occuring kmer lookup")
        return auto

    def match_kmers_multi(self, output_file, fasta_list, kmer_auto):
        logging.info("Writing output (pre-qualified-matched kmers)")
        with open(output_file, "w") as f:
            for record in fasta_list:
                match = self.__find_match_multi(str(record.seq), kmer_auto)
                if match:
                    line = str(match) + "\n"
                    f.write(line)
        logging.info("Completed determining the multi-occuring kmers that matched the pre-qualified minimal set")


if __name__ == '__main__':
    args = get_args()

    fasta_file = args.output + "/result_file.fasta"
    kmer_file = args.output + "/nr_more1List.txt"
    output_file = args.output + "/matchKmer4CleanKmer.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    multiOccuring = MultiOccurringPreMinSet()

    fasta_list, kmer_list = multiOccuring.load_data_multi(fasta_file, kmer_file)
    kmer_auto = multiOccuring.setup_automaton_multi(kmer_list)
    multiOccuring.match_kmers_multi(output_file, fasta_list, kmer_auto)