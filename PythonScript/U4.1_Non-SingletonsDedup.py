import argparse
import logging


def get_args():
    parser = argparse.ArgumentParser(prog = 'U4.1_Non-SingletonsDedup.py',
                    description = 'Step 4.1 of UNIQmin'
                    )
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    
    lines_seen = set()
    outfile = open(args.output + "/nr_more1List.txt", "w")
    for line in open(args.output + "/seqmore1List.txt", "r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
