import sys
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import math
import pandas as pd
import ahocorasick as ahc
import logging
import ast
import itertools
import os
import time

#--------#
# U1     #
#--------#

def generate_kmers(start, end):
    for record in fileA[start:end]:
        nr_sequence = record.seq
        seq_len = len(nr_sequence)
        kmer = int(sys.argv[2])
        count = 0
        temp = []
        for seq in list(range(seq_len-(kmer-1))):
            count += 1
            my_kmer = (nr_sequence[seq:seq+kmer])
            temp.append(str(my_kmer))
        with open(file_id, 'a') as f:
            f.writelines("%s\n" % kmer for kmer in temp)

#--------#
# U3.1   #
#--------#

class PreQualifiedMinSet:

    def load_data(self, fasta_file, kmer_file):
        logging.info("Loading fasta file for Pre Qualified Min Set")
        fasta_list = list(SeqIO.parse(fasta_file,"fasta"))
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

    def match_kmers(self, fasta_list, kmer_auto):
        logging.info("Writing output")
        with open(output_file,"w") as f:
            for record in fasta_list:
                match = self.__find_match(str(record.seq), kmer_auto)
                if match:
                    line = record.id + "\n"
                    f.write(line)
        logging.info("Completed for Pre Qualified Min Set")

#--------#
# U4.2   #
#--------#

class MultiOccurringPreMinSet:

    def load_data_multi(self, fasta_file, kmer_file):
        logging.info("Loading fasta file for Multi-Occuring Pre Min Set")
        fasta_list = list(SeqIO.parse(fasta_file,"fasta"))
        logging.info("Loading kmer list")
        kmer_list = [line.rstrip('\n') for line in open (kmer_file)]
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

    def match_kmers_multi(self, fasta_list, kmer_auto):
        logging.info("Writing output")
        with open(output_file, "w") as f:
            for record in fasta_list:
                match = self.__find_match_multi(str(record.seq), kmer_auto)
                if match:
                    line = str(match) + "\n"
                    f.write(line)
        logging.info("Completed for Multi-Occuring Pre Min Set")

#--------#
# U5.1   #
#--------#

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

# if __name__ == '__main__':

def main_match():
    #--------#
    # U1     #
    #--------#

    fileA = list(SeqIO.parse(sys.argv[1],"fasta"))
    file_id = "Output_kmers.txt"
    open(file_id, 'a').close()

    n = len(fileA)
    pool = ProcessPoolExecutor(int(sys.argv[3]))
    futures = []
    perCPUSize = math.ceil(n/int(sys.argv[3]))
    for i in range(0,int(sys.argv[3])):
        futures.append(pool.submit(generate_kmers, i * perCPUSize, (i+1) * perCPUSize))

    time.sleep(5)

    #--------#
    # U2.1   #
    #--------#

    #frequency count
    kmers = pd.read_csv("Output_kmers.txt", header=None)
    kmers.columns = ['kmer']
    kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

    #extract freq = 1
    #make it as a condition (eg: is_1) 
    #checking using boolean variable (eg: is_1.head())
    #filter rows for freq =1 using boolean variable
    is_1 = kmers['freq']==1
    is_1.head()
    kmer_1 = kmers[is_1]
    #check type of kmer_1 (eg: type(kmer_1))
    singleList = kmer_1['kmer']
    singleList.to_csv("seqSingleList.txt", index = False, header = False)

    #--------#
    # U2.2   #
    #--------#

    more1 = kmers['freq']!=1
    more1.head()
    kmer_more1 = kmers[more1]
    more1List = kmer_more1['kmer']
    more1List.to_csv("seqmore1List.txt", index = False, header = False)

    #--------#
    # U3.1   #
    #--------#

    fasta_file = sys.argv[1]
    kmer_file = "seqSingleList.txt"
    output_file = "seqfileZ.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    preQualified = PreQualifiedMinSet()

    fasta_list, kmer_list = preQualified.load_data(fasta_file, kmer_file)
    kmer_auto = preQualified.setup_automaton(kmer_list)
    preQualified.match_kmers(fasta_list, kmer_auto)

    #--------#
    # U3.2   #
    #--------#

    fileA = list(SeqIO.parse(sys.argv[1],"fasta"))
    header_set = set(line.strip() for line in open("seqfileZ.txt"))
    remainingSeq = open("remainingSeq.fasta","w")

    for seq_record in fileA:
        try:
            header_set.remove(seq_record.name)
        except KeyError:
            remainingSeq.write(seq_record.format("fasta"))
    remainingSeq.close()

    fasta_file = sys.argv[1] 
    wanted_file = "seqfileZ.txt" 
    result_file = "result_file.fasta" 

    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open(result_file, "w") as f:
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], f, "fasta")

    #--------#
    # U4.1   #
    #--------#

    lines_seen = set()
    outfile = open("nr_more1List.txt","w")
    for line in open("seqmore1List.txt","r"):
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()

    #--------#
    # U4.2   #
    #--------#

    fasta_file = "result_file.fasta"
    kmer_file = "nr_more1List.txt"
    output_file = "matchKmer4CleanKmer.txt"

    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    multiOccuring = MultiOccurringPreMinSet()

    fasta_list, kmer_list = multiOccuring.load_data_multi(fasta_file, kmer_file)
    kmer_auto = multiOccuring.setup_automaton_multi(kmer_list)
    multiOccuring.match_kmers_multi(fasta_list, kmer_auto)

    #--------#
    # U4.3   #
    #--------#

    listOfLines = list()        
    with open ("matchKmer4CleanKmer.txt", "r") as myfile:
        for line in myfile:
            line = ast.literal_eval(line)
            listOfLines.append(line) 

    full_list = list(itertools.chain(*listOfLines))

    with open("fullList.txt",'w') as f:
        for item in full_list:
            f.write("%s\n" % item)

    lines_seen = set()
    nr_lines = open("Clean_lines.txt", "w")
    for line in open("fullList.txt","r"):
        if line not in lines_seen:
            nr_lines.write(line)
            lines_seen.add(line)
    nr_lines.close()

    a = open("nr_more1List.txt", 'r')
    b = open("Clean_lines.txt", 'r')
    result = "remainingKmer.txt"

    remain_kmer_list = list(set(a) - set(b))
    with open(result, "w") as f: 
        for i in remain_kmer_list:
            f.write(i)

    #--------#
    # U5.1   #
    #--------#

    os.system('cp seqfileZ.txt fileZ.txt')
    os.mkdir('match')

    remaining = RemainingMinSet()

    remain_Seq = list(SeqIO.parse("remainingSeq.fasta","fasta"))
    remain_kmer = [line.rstrip('\n') for line in open ("remainingKmer.txt")]
    remain_Seq_copy = remain_Seq.copy()

    a = 0

    while(len(remain_kmer) != 0):
        
        A = remaining.make_automaton(remain_kmer)
        
        matching_file = 'match/matching'+str(a)
        remain_kmer_file = 'match/remain_kmer'+str(a)
        
        # save matching to file
        with open(matching_file, 'w') as f:
            for index in range(len(remain_Seq)):
                x = remain_Seq[index].id
                y = remaining.find_matching(remain_Seq[index].seq, A)
                z = len(y)
                f.write(x + ';' + str(y) + ';' + str(z) + '\n')
                if z == 0:
                    for i in range(len(remain_Seq_copy)):
                        if remain_Seq_copy[i].id == x:
                            del remain_Seq_copy[i]
                            break                
        
        remain_Seq = remain_Seq_copy.copy()
        
        # read matching file and sorted by descending & some cleaning
        df = pd.read_csv(matching_file, delimiter=';', names=['sequence_id', 'matched_kmer', 'count']).sort_values(by='count',ascending=False, kind='mergesort')
        df['matched_kmer'] = df['matched_kmer'].str.replace(r"\[|\]|'","")
        
        # save highest count id to file
        fileZ = open('fileZ.txt', 'a')
        fileZ.write(df['sequence_id'].iloc[0] + '\n')
        
        # remove highest count kmer
        kmer_to_remove = df['matched_kmer'].iloc[0].split(', ')
        remain_kmer = list(set(remain_kmer) - set(kmer_to_remove))
        
        # save remain kmer to file
        with open(remain_kmer_file, 'w') as f:
            for i in remain_kmer:
                f.write(i + '\n')
        
        a = a + 1

    #--------#
    # U5.2   #
    #--------#

    fasta_file =  sys.argv[1] # Input fasta file
    wanted_file = "fileZ.txt" # Input interesting sequence IDs, one per line
    result_file = "FileZ.fasta" # Output fasta file

    wanted = set()
    with open (wanted_file) as f: 
        for line in f: 
            line = line.strip()
            if line != "":
                wanted.add(line)

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    with open (result_file, "w") as f: 
        for seq in fasta_sequences: 
            if seq.id in wanted: 
                SeqIO.write([seq], f, "fasta")


