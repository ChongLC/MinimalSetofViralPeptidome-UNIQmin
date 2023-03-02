import argparse
import logging
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(prog = 'U3.2_UnmatchedSingletons.py',
                    description = 'Step 3.2 of UNIQmin'
                    )
    parser.add_argument('-i', '--input', dest="input", help='Path of the input file (in FASTA format)')
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    fileA = list(SeqIO.parse(args.input, "fasta"))
    header_set = set(line.strip() for line in open(args.output + "/seqfileZ.txt"))
    remainingSeq = open(args.output + "/remainingSeq.fasta", "w")

    for seq_record in fileA:
        try:
            header_set.remove(seq_record.name)
        except KeyError:
            remainingSeq.write(seq_record.format("fasta"))
    remainingSeq.close()

    fasta_file = args.input
    wanted_file = args.output + "/seqfileZ.txt"
    result_file = args.output + "/result_file.fasta"

    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    seq_list = []
    for seq in fasta_sequences:
        if seq.id in wanted:
            seq_list.append(seq)
    with open(result_file, "w") as f:
        SeqIO.write(seq_list, f, "fasta")