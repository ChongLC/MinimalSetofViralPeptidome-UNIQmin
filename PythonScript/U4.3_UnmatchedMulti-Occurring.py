import argparse
import logging
import ast
import itertools


def get_args():
    parser = argparse.ArgumentParser(prog='U4.3_UnmatchedMulti-Occurring.py',
                                    description='Step 4.3 of UNIQmin')
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    lines = []
    with open(f"{args.output}/matchKmer4CleanKmer.txt", "r") as myfile:
        for line in myfile:
            line = ast.literal_eval(line)
            lines.append(line)

    full_list = list(itertools.chain(*lines))
    lines_seen = set()
    clean_lines = []
    with open(f"{args.output}/fullList.txt", "w") as f:
        f.write('\n'.join(i for i in full_list))

    with open(f"{args.output}/fullList.txt", "r") as f:
        for line in f.readlines():
            if line not in lines_seen:
                clean_lines.append(line)
                lines_seen.add(line)

    with open(f"{args.output}/Clean_lines.txt", "w") as nr_lines:
        nr_lines.write(''.join(i for i in clean_lines))

    with open(f"{args.output}/nr_more1List.txt", "r") as a, open(f"{args.output}/Clean_lines.txt", "r") as b, open(
        f"{args.output}/remainingKmer.txt", "w") as result:
        remain_kmer_list = list(set(a) - set(b))
        result_kmer_list = list()
        for i in remain_kmer_list:
            result_kmer_list.append(i)
        result.write(''.join(i for i in result_kmer_list))