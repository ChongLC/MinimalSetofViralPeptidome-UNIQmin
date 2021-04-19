import pandas as pd

kmers = pd.read_csv("Output_kmers.txt", header=None)
kmers.columns = ['kmer']
kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

more1 = kmers['freq']!=1
more1.head()
kmer_more1 = kmers[more1]
more1List = kmer_more1['kmer']
more1List.to_csv("seqmore1List.txt", index = False, header = False)
