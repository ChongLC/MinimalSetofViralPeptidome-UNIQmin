import argparse, sys
import string
import random

#--------------#
# Argument     #
#--------------#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output',dest="output", help='Path of the output file to be created')
    parser.add_argument('-l', '--seqlen', dest="seqlen", help='The length of random protein sequences to be generated')
    parser.add_argument('-n', '--seqnum', dest="seqnum", help='The number of random protein sequences to be generated')

    return parser.parse_args()

if __name__ == '__main__':

    args = get_args()

    seqlen = int(args.seqlen)
    seqnum = int(args.seqnum)
    counter = 0
    aa_fullList = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X']
    random_list = []

    with open(args.output, 'w') as f:
        for i in range(seqnum):
            counter += 1
            randomSeqList = ''.join(random.choices(aa_fullList, weights = [6.64, 4.22, 5.2, 5.04, 2.56, 3.72, 5.28, 6.2, 1.92, 5.55, 9.13, 5.86, 2.33, 4.45, 4.28, 6.7, 7.07, 1.35, 3.99, 7.34, 1.16], k=seqlen))
            random_list.append(randomSeqList)
            i+=1
            f.write(f">Sequence {counter}" + '\n' + randomSeqList + '\n')
