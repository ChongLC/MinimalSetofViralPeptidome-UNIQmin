lines_seen = set()
outfile = open("nr_more1List.txt","w")
for line in open("seqmore1List.txt","r"):
	if line not in lines_seen:
		outfile.write(line)
		lines_seen.add(line)
outfile.close()
