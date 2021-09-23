# **UNIQmin: An alignment-independent tool for the study of pathogen sequence diversity at any given rank of taxonomy lineage**

### Brief Description
Sequence variation among pathogens, even of a single amino acid, can expand their host repertoire or enhance the infection ability. Alignment independent approach represents an alternative approach to the study of pathogen diversity, which is devoid of the need for sequence conservation to perform comparative analyses. Herein, we present UNIQmin, a tool that utilises an alignment independent method to generate the minimal set of pathogen sequences, as a way to study their diversity, across any rank of taxonomic lineage. The minimal set refers to the smallest possible number of sequences required to capture the entire repertoire of pathogen peptidome diversity present in a sequence dataset.

---
Table of Contents
====================
- [Step-by-step of UNIQmin](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonCode)
- [Figure Scheme](#figure-scheme)
- [UNIQmin as a Pipeline](#uniqmin-as-a-pipeline)
- [Citing Resources](#citing-resources)

## Step-by-step of UNIQmin
Please refer to the [PythonCode](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonCode) folder. 

## Figure Scheme
![Scheme1](https://user-images.githubusercontent.com/51225708/134200760-f70d72ee-0fac-4535-aa3e-61fd8dd1a69f.png)

## UNIQmin as a Pipeline

#### Shell Version
As visualised above, UNIQmin comprises of five steps with respective python scripts employed according to the order of step (server specs: Intel(R) Xeon(R) E5-2690 v2 @ 3.00GHz 40-core processors, 396 GB of RAM and 44 TB of local storage. The single pipeline shell script (UNIQmin.sh), sample input file (exampleinput.fas) and example output (exampleoutput.txt) are provided. 

`UNIQmin.sh`

#### Python Version
`python UNIQmin.py -i exampleinput.fas -o example -k 9 -cpu 14`

## UNIQmin as a Package

#### Installation
`pip install UNIQmin`

#### Usage
`UNIQmin [-i INPUT] [-o OUTPUT] [-k [KMERLENGTH]] [-cpu [CPUSIZE]]`

For example, UNIQmin tool is applied to generate a minimal set (example) with a sample input file (exampleinput.fas). A *k*-mer window size of nine (9; nonamer) is used with utilising 14-cores. 

`UNIQmin -i exampleinput.fas -o example -k 9 -cpu 14`

#### Command-line Arguments
| Argument 	| Parameter              | Type    	| Default 	| Description                                |           
|----------	|----------------------- |---------	|---------- |--------------------------------------------|
| -h       	| help                   | N/A     	| N/A     	| Show this help message and exit            |
| -i       	| sequence input file    | String  	| N/A     	| Path of the input file (in FASTA format)   |
| -o       	| output directory name  | String  	| N/A     	| Path of the output file to be created      |
| -k        | *k*-mer window size    | Integer 	| 9       	| The length of *k*-mers to be used          |
| -cpu      | cpu size               | Integer 	| 14       	| The number of CPU cores to be used         |

---
### Citing Resources
1. Chong, L.C.; Lim, W.L.; Ban, K.H.K.; Khan, A.M. An Alignment-Independent Approach for the Study of Viral Sequence Diversity at Any Given Rank of Taxonomy Lineage. Biology 2021, 10, 853. https://doi.org/10.3390/biology10090853
