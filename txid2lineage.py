import csv
import time
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

if __name__ == '__main__':
    file = open("testingtxid.txt", "r")    
    taxids = []
    for line in file:
        line = line.split("\n")[0]
        taxids.append(line.split(",")[0])

    desired_ranks = ['superkingdom','family','genus','species']
    results = list()
    for taxid in taxids:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_desired_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>':
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)

    i = 0
    fileR = open('fileR_2db.txt','a')
    for result in results:
        fileR.write(','.join(result) + '\n')
      
        i += 1

    fileR.close()
    file.close()
