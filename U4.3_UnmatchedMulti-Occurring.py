#Correct code 
import ast
import itertools

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
