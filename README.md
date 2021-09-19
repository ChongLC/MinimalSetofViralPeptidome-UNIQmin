# **UNIQmin: An alignment-independent tool for the study of pathogen sequence diversity at any given rank of taxonomy lineage**

### Brief Description
Sequence variation among pathogens, even of a single amino acid, can expand their host repertoire or enhance the infection ability. Alignment independent approach represents an alternative approach to the study of pathogen diversity, which is devoid of the need for sequence conservation to perform comparative analyses. Herein, we present UNIQmin, a tool that utilises an alignment independent method to generate the minimal set of pathogen sequences, as a way to study their diversity, across any rank of taxonomic lineage. The minimal set refers to the smallest possible number of sequences required to capture the entire repertoire of pathogen peptidome diversity present in a sequence dataset.

Table of Contents
====================
- [Step-by-step of UNIQmin](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonCode)
- [Figure Scheme](#figure-scheme)
- [UNIQmin as a Pipeline](#uniqmin-as-a-pipeline)
- [Citing Resources](#citing-resources)

## Step by step of UNIQmin
Please refer to the [PythonCode](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonCode) folder. 

## Figure Scheme
<img src="Scheme.png" width="1204" height="676">

---
## UNIQmin as a Pipeline
As explained above, UNIQmin comprises of five steps with respective python scripts employed according to the order of step (server specs: Intel(R) Xeon(R) E5-2690 v2 @ 3.00GHz 40-core processors, 396 GB of RAM and 44 TB of local storage. The sample input file (inputfile.fas) and example output (exampleoutput.txt) are provided. 
```
#!/bin/bash
#$ -V

dir=/backup/user/ext/perdana/lichuin/cdhitObj2/testingStitchScript/ #path of your directory
cd $dir
python U1_KmerGenerator.py
wait 
python U2.1_Singletons.py
wait
python U2.2_Multitons.py
wait
python U3.1_PreQualifiedMinSet.py
wait 
python U3.2_UnmatchedSingletons.py
wait 
python U4.1_Non-SingletonsDedup.py
wait 
python U4.2_Multi-OccurringPreMinSet.py
wait 
python U4.3_UnmatchedMulti-Occurring.py
wait
cp seqfileZ.txt fileZ.txt
mkdir match
wait 
python U5.1_RemainingMinSet.py
wait
python U5.2_MinSet.py
```

---
## Citing Resources
1. Chong, L.C.; Lim, W.L.; Ban, K.H.K.; Khan, A.M. An Alignment-Independent Approach for the Study of Viral Sequence Diversity at Any Given Rank of Taxonomy Lineage. Biology 2021, 10, 853. https://doi.org/10.3390/biology10090853
