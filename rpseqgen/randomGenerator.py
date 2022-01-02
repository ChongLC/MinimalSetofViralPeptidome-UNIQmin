import argparse, sys
import string
import random

#--------------#
# Argument     #
#--------------#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output',dest="output", help='Path of the output file to be created')
    parser.add_argument('-l', '--seqlen', dest="seqlen", help='The length of protein sequences')
    parser.add_argument('-n', '--seqnum', dest="seqnum", help='The number of protein sequences')

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
            randomSeqList = ''.join(random.choices(aa_fullList, weights = [6.42, 5.51, 5.08, 5.01, 1.97, 3.96, 6.02, 6.76, 2.04, 6.28, 8.61, 5.93, 2.55, 3.77, 4.99, 6.85, 6.68, 1.76, 3.28, 6.34, 0.19], k=seqlen))
            random_list.append(randomSeqList)
            i+=1
            f.write(f">Sequence {counter}" + '\n' + randomSeqList + '\n')
