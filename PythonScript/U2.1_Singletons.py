import argparse
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(prog = 'U2.1_Singletons.py',
                    description = 'Step 2.1 of UNIQmin'
                    )
    parser.add_argument('-o', '--output', dest="output", default=".", help='Path of the output file to be created')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    
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
