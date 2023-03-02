import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(prog = 'U2.2_Multitons.py',
                    description = 'Step 2.2 of UNIQmin'
                    )
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    
    kmers = pd.read_csv(args.output + "/Output_kmers.txt", header=None)
    kmers.columns = ['kmer']
    kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

    more1 = kmers['freq'] != 1
    more1.head()
    kmer_more1 = kmers[more1]
    more1List = kmer_more1['kmer']
    more1List.to_csv(args.output + "/seqmore1List.txt", index=False, header=False)
