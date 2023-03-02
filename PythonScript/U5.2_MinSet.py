import argparse
import logging
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(prog = 'U5.2_MinSet.py',
                    description = 'Step 5.2 of UNIQmin'
                    )
  
    parser.add_argument('-i', '--input', dest="input", help='Path of the input file (in FASTA format)')
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')
  
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()

    fasta_file = args.input  # Input fasta file
    wanted_file = args.output + "/minimalSet/fileZ.txt"  # Input interesting sequence IDs, one per line
    result_file = args.output + "/minimalSet/fileZ.fasta"  # Output fasta file

    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    with open(result_file, "w") as f:
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")