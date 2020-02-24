from Bio import SeqIO
import pandas as pd

#frequency count
kmers = pd.read_csv("Output_kmers_Paramyxoviridae.txt", header=None)
kmers.columns = ['kmer']
kmers['freq'] = kmers.groupby('kmer')['kmer'].transform('count')

#extract freq = 1
#make it as a condition (eg: is_1) 
#checking using boolean variable (eg: is_1.head())
#filter rows for freq =1 using boolean variable
is_1 = kmers['freq']==1
#is_1.head()
kmer_1 = kmers[is_1]
#check type of kmer_1 (eg: type(kmer_1))
singleList = kmer_1['kmer']
singleList.to_csv("seqSingleList.txt", index = False, header = False)
