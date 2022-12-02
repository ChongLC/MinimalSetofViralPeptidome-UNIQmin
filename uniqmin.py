import argparse
# logging.basicConfig(level = logging.INFO)
import ast
import glob
import itertools
import logging
import os
import queue
import threading

import ahocorasick as ahc
import pandas as pd
from Bio import SeqIO


# --------------#
# Argument     #
# --------------#
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest="input", help='Path of the input file (in FASTA format)')
    parser.add_argument('-o', '--output', dest="output", help='Path of the output file to be created')
    parser.add_argument('-k', '--kmer', dest="kmerlength", help='The length of k-mers to be used', default=9, nargs='?')
    parser.add_argument('-cpu', '--cpusize', dest="cpusize", help='The number of CPU cores to be used', default=14,
                        nargs='?')
    parser.add_argument('-t', '--thread', dest="threads", help="The number of threads to be used", default=6, type=int)

    return parser.parse_args()


# --------#
# U1     #
# --------#

def generate_kmers(thread_id, seq_queue: queue.Queue, kmerlength, path):
    while True:
        try:
            for record in seq_queue.get_nowait():
                nr_sequence = record.seq
                seq_len = len(nr_sequence)
                kmer = int(kmerlength)
                count = 0
                temp = []
                for seq in list(range(seq_len - (kmer - 1))):
                    count += 1
                    my_kmer = (nr_sequence[seq:seq + kmer])
                    temp.append(str(my_kmer)+'\n')
                open(path + f'/kmers_{thread_id}.tmp', 'a').writelines(temp)
        except queue.Empty:
            break
    return True


# --------#
# U3.1   #
# --------#

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
        logging.info("Writing output")
        with open(output_file, "w") as f:
            for record in fasta_list:
                match = self.__find_match(str(record.seq), kmer_auto)
                if match:
                    line = record.id + "\n"
                    f.write(line)
        logging.info("Completed determining pre-qualified minimal set")


# --------#
# U4.2   #
# --------#

class MultiOccurringPreMinSet:

    def load_data_multi(self, fasta_file, kmer_file):
        logging.info("Loading fasta file for minimal set from multi-occuring kmer list")
        fasta_list = list(SeqIO.parse(fasta_file, "fasta"))
        logging.info("Loading kmer list")
        kmer_list = [line.rstrip('\n') for line in open(kmer_file)]
        return fasta_list, kmer_list

    def __find_match_multi(self, line, A):
        found_kmers = []
        for end_index, kmer in A.iter(line):
            found_kmers.append(kmer)
        return found_kmers

    def setup_automaton_multi(self, kmer_list):
        logging.info("Setting up kmer lookup")
        auto = ahc.Automaton()
        for seq in kmer_list:
            auto.add_word(seq, seq)
        auto.make_automaton()
        logging.info("Completed set-up of kmer lookup")
        return auto

    def match_kmers_multi(self, output_file, fasta_list, kmer_auto):
        logging.info("Writing output")
        with open(output_file, "w") as f:
            for record in fasta_list:
                match = self.__find_match_multi(str(record.seq), kmer_auto)
                if match:
                    line = str(match) + "\n"
                    f.write(line)
        logging.info("Completed determining the multi-occuring kmers that matched the pre-qualified minimal set")


# --------#
# U5.1   #
# --------#

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


# --------#
# Thread pool  #
# --------#
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


class GenerateKmersPool(threading.Thread):
    def __init__(self, thread_id, seq_queue, kmerlength, path):
        super().__init__()
        self.thread_id = thread_id
        self.seq_queue = seq_queue
        self.kmerlength = kmerlength
        self.path = path

    def run(self):
        generate_kmers(self.thread_id, self.seq_queue, self.kmerlength, self.path)


def copy_file(source, dest, buffer_size=10 * 1024 * 1024):
    """
    Copy a file from source to dest. source and dest
    must be file-like objects, i.e. any object with a read or
    write method, like for example StringIO.
    """
    while True:
        copy_buffer = source.read(buffer_size)
        if not copy_buffer:
            break
        dest.write(copy_buffer)


def main():
    args = get_args()

    # --------#
    # U1     #
    # --------#

    os.mkdir(args.output)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    fileA = list(SeqIO.parse(args.input, "fasta"))
    file_id = args.output + "/Output_kmers.txt"
    open(file_id, 'a').close()

    n = len(fileA)
    chunk_size = 2000
    start = 0
    stop = chunk_size
    threads = list()
    logging.info("Loading fasta file to generate kmer list")
    seq_queue = queue.Queue()
    while start < n + 1:
        seq_queue.put(fileA[start:stop])
        start += chunk_size
        stop = min(n + 1, stop + chunk_size)
    for i in range(args.threads):
        thread = GenerateKmersPool(i, seq_queue, args.kmerlength, args.output)
        thread.start()
        threads.append(thread)

    for t in threads:
        t.join()
    # Combining kmers temps into a single file
    files = glob.glob(args.output + '/kmers_*.tmp')
    for file in files:
        with open(file, 'r') as src, open(file_id, 'a+') as dst:
            copy_file(src, dst)
            dst.write('\n')
        os.remove(file)
    # time.sleep(60)

    # --------#
    # U2.1   #
    # --------#

    # frequency count
    kmers = pd.read_csv(args.output + "/Output_kmers.txt", header=None)
    kmers.columns = ['kmer']
    kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

    # extract freq = 1
    # make it as a condition (eg: is_1)
    # checking using boolean variable (eg: is_1.head())
    # filter rows for freq =1 using boolean variable
    is_1 = kmers['freq'] == 1
    is_1.head()
    kmer_1 = kmers[is_1]
    # check type of kmer_1 (eg: type(kmer_1))
    singleList = kmer_1['kmer']
    singleList.to_csv(args.output + "/seqSingleList.txt", index=False, header=False)

    # --------#
    # U2.2   #
    # --------#

    more1 = kmers['freq'] != 1
    more1.head()
    kmer_more1 = kmers[more1]
    more1List = kmer_more1['kmer']
    more1List.to_csv(args.output + "/seqmore1List.txt", index=False, header=False)

    # --------#
    # U3.1   #
    # --------#

    fasta_file = args.input
    kmer_file = args.output + "/seqSingleList.txt"
    output_file = args.output + "/seqfileZ.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    preQualified = PreQualifiedMinSet()

    fasta_list, kmer_list = preQualified.load_data(fasta_file, kmer_file)
    kmer_auto = preQualified.setup_automaton(kmer_list)
    preQualified.match_kmers(output_file, fasta_list, kmer_auto)

    # --------#
    # U3.2   #
    # --------#

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

    # --------#
    # U4.1   #
    # --------#

    lines_seen = set()
    outfile = open(args.output + "/nr_more1List.txt", "w")
    for line in open(args.output + "/seqmore1List.txt", "r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()

    # --------#
    # U4.2   #
    # --------#

    fasta_file = args.output + "/result_file.fasta"
    kmer_file = args.output + "/nr_more1List.txt"
    output_file = args.output + "/matchKmer4CleanKmer.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    multiOccuring = MultiOccurringPreMinSet()

    fasta_list, kmer_list = multiOccuring.load_data_multi(fasta_file, kmer_file)
    kmer_auto = multiOccuring.setup_automaton_multi(kmer_list)
    multiOccuring.match_kmers_multi(output_file, fasta_list, kmer_auto)

    # --------#
    # U4.3   #
    # --------#

    listOfLines = list()
    with open(args.output + "/matchKmer4CleanKmer.txt", "r") as myfile:
        for line in myfile:
            line = ast.literal_eval(line)
            listOfLines.append(line)

    full_list = list(itertools.chain(*listOfLines))
    lines_seen = set()
    lines_list = list()
    with open(args.output + "/fullList.txt", 'w') as f:
        f.write('\n'.join(i for i in full_list))

    with open(args.output + "/fullList.txt", 'r') as f:
        for line in f.readlines():
            if line not in lines_seen:
                lines_list.append(line)
                lines_seen.add(line)

    with open(args.output + "/Clean_lines.txt", "w") as nr_lines:
        nr_lines.write(''.join(i for i in lines_list))

    with open(args.output + "/nr_more1List.txt", 'r') as a, open(args.output + "/Clean_lines.txt", 'r') as b, open(
            args.output + "/remainingKmer.txt", 'w') as result:
        remain_kmer_list = list(set(a) - set(b))
        result_kmer_list = list()
        for i in remain_kmer_list:
            result_kmer_list.append(i)
        result.write(''.join(i for i in result_kmer_list))

    # --------#
    # U5.1   #
    # --------#
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

    # --------#
    # U5.2   #
    # --------#
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

    logging.info("Successfully generated final minimal set")


if __name__ == '__main__':
    main()
