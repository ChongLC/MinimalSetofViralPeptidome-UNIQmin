import argparse
import logging
import os
import queue
import threading
import ahocorasick as ahc
import pandas as pd
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(prog = 'U5.1_RemainingMinSet.py',
                    description = 'Step 5.1 of UNIQmin'
					)
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')
    parser.add_argument('-k', '--kmer', dest="kmerlength", help='The length of k-mers to be used', default=9, nargs='?')
    parser.add_argument('-t', '--thread', dest="threads", help="The number of threads to be used", default=4, type=int)

    return parser.parse_args()


class RemainingMinSet:

    def make_automaton(self, kmer_list):
        A = ahc.Automaton()
        for kmer in kmer_list:
            A.add_word(kmer, kmer)
        A.make_automaton()
        return A

    def find_matching(self, line, A):
        found_kmers = []
        for end_index, kmer in A.iter(str(line)):
            found_kmers.append(kmer)
        return found_kmers


def proccess_match(remain_seq_queue: queue.Queue, remain_Seq_copy, remaining, A, lines: list):
    while True:
        try:
            remain_Seq = remain_seq_queue.get_nowait()
            x = remain_Seq.id
            y = remaining.find_matching(remain_Seq.seq, A)
            z = len(y)
            lines.append(x + ';' + str(y) + ';' + str(z) + '\n')
            if z == 0:
                for i in range(len(remain_Seq_copy)):
                    if remain_Seq_copy[i].id == x:
                        del remain_Seq_copy[i]
                        break
        except queue.Empty:
            break
    return True


class FindMatchPool(threading.Thread):
    def __init__(self, thread_id, remain_seq_queue, remain_Seq_copy, remaining, A, lines: list):
        super().__init__()
        self.thread_id = thread_id
        self.remain_seq_queue = remain_seq_queue
        self.remain_Seq_copy = remain_Seq_copy
        self.remaining = remaining
        self.A = A
        self.lines = lines

    def run(self):
        proccess_match(self.remain_seq_queue, self.remain_Seq_copy, self.remaining, self.A, self.lines)

if __name__ == '__main__':
    args = get_args()

    os.mkdir(args.output + '/minimalSet')
    os.system(f"cp {args.output}/seqfileZ.txt {args.output}/minimalSet/fileZ.txt")
    os.mkdir(args.output + '/match')

    remaining = RemainingMinSet()

    remain_Seq = list(SeqIO.parse(args.output + "/remainingSeq.fasta", "fasta"))
    remain_kmer = [line.rstrip('\n') for line in open(args.output + "/remainingKmer.txt")]
    remain_Seq_copy = remain_Seq.copy()

    a = 0

    while (len(remain_kmer) != 0):

        A = remaining.make_automaton(remain_kmer)

        matching_file = args.output + '/match/matching' + str(a)
        remain_kmer_file = args.output + '/match/remain_kmer' + str(a)
        remain_Seq_queue = queue.Queue()
        threads = []
        lines = []
        # Put queue for threading
        for index in range(len(remain_Seq)):
            remain_Seq_queue.put(remain_Seq[index])
        for i in range(args.threads):
            thread = FindMatchPool(i, remain_Seq_queue, remain_Seq_copy, remaining, A, lines)
            thread.start()
            threads.append(thread)
        for t in threads:
            t.join()
        # save matching to file
        with open(matching_file, 'w') as f:
            f.write(''.join(i for i in lines))

        remain_Seq = remain_Seq_copy.copy()
        # read matching file and sorted by descending & some cleaning
        df = pd.read_csv(matching_file, delimiter=';', names=['sequence_id', 'matched_kmer', 'count']).sort_values(
            by='count', ascending=False, kind='mergesort')
        df['matched_kmer'] = df['matched_kmer'].str.replace(r"\[|\]|'", "")

        # save highest count id to file
        fileZ = open(args.output + '/minimalSet/fileZ.txt', 'a')
        fileZ.write(df['sequence_id'].iloc[0] + '\n')

        # remove highest count kmer
        kmer_to_remove = df['matched_kmer'].iloc[0].split(', ')
        remain_kmer = list(set(remain_kmer) - set(kmer_to_remove))

        # save remain kmer to file
        with open(remain_kmer_file, 'w') as f:
            for i in remain_kmer:
                f.write(i + '\n')

        a = a + 1