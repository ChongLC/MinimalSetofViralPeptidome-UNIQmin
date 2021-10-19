from Bio import SeqIO
import sys
import string
import random

fastafile = sys.argv[1]
outputfile = sys.argv[2]
seq_list = ""
random_list = []
percent_dict = []
counter = 0
SL = int(sys.argv[3])
SN = int(sys.argv[4])

for sequence in SeqIO.parse(fastafile, "fasta"):
  seq_list += str(sequence.seq)

aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'O', 'U', 'Z', 'X']
def CalculateFrequency(seq_list):
    Length = len(seq_list)
    ProteinDict = {}
    for aa in aa_list:
         ProteinDict[aa] = seq_list.count(aa)
    return(ProteinDict)
  
def CalculatePercent(seq_list,ProteinDict):
  PercentDict = {}
  Length = len(seq_list)
  for aa in ProteinDict:
     PercentDict[aa] = round((ProteinDict[aa] / float(len(seq_list)))*100,2)
  return(PercentDict)

def PrintOutput(Dict):
  for aa in Dict:
    print(aa,":",round(Dict[aa],2))

def PrintOutput2(Dict):
  for aa in Dict:
    percent_dict.append(round(Dict[aa],2))

def Main():
  ProteinDict = CalculateFrequency(seq_list)
  PercentDict = CalculatePercent(seq_list,ProteinDict)
  PrintOutput2(PercentDict)
if __name__=='__main__':
  Main()

with open(outputfile, 'w') as f:
  for i in range(SN): #SN
    counter += 1
    randomSeqList = ''.join(random.choices(aa_list, weights = percent_dict, k=SL))
    f.write(f">Sequence {counter}" + '\n' + randomSeqList + '\n')
    i+=1