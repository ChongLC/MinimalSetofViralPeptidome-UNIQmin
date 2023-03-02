import glob
import os
import queue
import threading

from Bio import SeqIO

fileA = list(SeqIO.parse(r"D:\Sources\Python\MinimalSetofViralPeptidome-UNIQmin\exampleinput.fas", "fasta"))
file_id = "Output_kmers.txt"


class GenerateKmersPool(threading.Thread):
    def __init__(self, thread_id, seq_queue, kmerlength, path):
        super().__init__()
        self.thread_id = thread_id
        self.seq_queue = seq_queue
        self.kmerlength = kmerlength
        self.path = path

    def run(self):
        generate_kmers(self.thread_id, self.seq_queue, self.kmerlength, self.path)


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
                    temp.append(str(my_kmer) + '\n')
                open(path + f'/kmers_{thread_id}.tmp', 'a').writelines(temp)
        except queue.Empty:
            break
    return True


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


if __name__ == '__main__':
    path_output = 'results'
    n_thread = 6
    kmerlength = 9
    n = len(fileA)
    chunk_size = 2000
    start = 0
    stop = chunk_size
    threads = list()
    seq_queue = queue.Queue()
    while start < n + 1:
        seq_queue.put(fileA[start:stop])
        start += chunk_size
        stop = min(n + 1, stop + chunk_size)
    for i in range(n_thread):
        thread = GenerateKmersPool(i, seq_queue, kmerlength, path_output)
        thread.start()
        threads.append(thread)

    for t in threads:
        t.join()
    # Combining kmers temps into a single file
    files = glob.glob(path_output + '/kmers_*.tmp')
    for file in files:
        with open(file, 'r') as src, open(file_id, 'a+') as dst:
            copy_file(src, dst)
            dst.write('\n')
        os.remove(file)
