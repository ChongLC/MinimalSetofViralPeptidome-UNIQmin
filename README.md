# **UNIQmin: An alignment-independent tool for the study of pathogen sequence diversity at any given rank of taxonomy lineage**

[![DOI - 10.3390/biology10090853](https://img.shields.io/badge/DOI-10.3390%2Fbiology10090853-2ea44f)](https://doi.org/10.3390/biology10090853)
[![PyPI](https://img.shields.io/pypi/v/uniqmin?logo=pypi)](https://pypi.org/project/uniqmin/)
[![GitHub tag](https://img.shields.io/github/tag/ChongLC/MinimalSetofViralPeptidome-UNIQmin)](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/releases/?include_prereleases&sort=semver "View GitHub releases")
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)

<!--
![Visitor](https://visitor-badge.laobi.icu/badge?page_id=https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin)
-->

### Brief Description
Sequence variation among pathogens, even of a single amino acid, can expand their host repertoire or enhance the infection ability. Alignment independent approach represents an alternative approach to the study of pathogen diversity, which is devoid of the need for sequence conservation to perform comparative analyses. Herein, we present UNIQmin, a tool that utilises an alignment independent method to generate the minimal set of pathogen sequences, as a way to study their diversity, across any rank of taxonomic lineage. The minimal set refers to the smallest possible number of sequences required to capture the entire repertoire of pathogen peptidome diversity present in a sequence dataset.

---
Table of Contents
====================
- [Step-by-step of UNIQmin](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonScript)
- [Figure Scheme](#figure-scheme)
- [UNIQmin as a Pipeline](#uniqmin-as-a-pipeline)
- [Generate a random protein sequence dataset](#generate-a-random-protein-sequence-dataset)
- [Citing Resources](#citing-resources)
- [Found a Bug](#found-a-bug)

## Step-by-step of UNIQmin
Please refer to the [PythonScript](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonScript) folder. 

## Figure Scheme
![uniqminScheme](https://user-images.githubusercontent.com/51225708/152393757-5be032bc-c17a-49b2-b9b7-aa3637aed1e1.png)


## UNIQmin as a Pipeline

#### Shell Version
As visualised above, UNIQmin comprises of five steps with respective python scripts employed according to the order of step (server specs: Intel(R) Xeon(R) E5-2690 v2 @ 3.00GHz 40-core processors, 396 GB of RAM and 44 TB of local storage. The single pipeline shell script (UNIQmin.sh), sample input file (exampleinput.fas) and example output (exampleoutput.fasta) are provided. 

```
uniqmin.sh
```

#### Python Version
```
python uniqmin.py -i exampleinput.fas -o example -k 9 -cpu 14
```

## UNIQmin as a Package

#### Installation
* via pip <br>
  ```
  pip install uniqmin
  ```
* via package clone from GitHub repository
  ```
  git clone https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin.git
  ```
  <br>
  
  **Note for user who uses conda environment (*e.g.*: jupyter notebook):** <br>
  Before `pip` installing the package, run <br>
  ```
  conda config --add channels conda-forge
  conda install pyahocorasick
  ```
  ... and restart the kernel to use the updated package. Then, run 
  ```
  pip install uniqmin
  ```
  
#### Upgrade installed version
```
pip install uniqmin --upgrade
```

#### Usage
`uniqmin [-h] [-i INPUT] [-o OUTPUT] [-k KMERLENGTH] [-cpu CPUSIZE]`

For example, UNIQmin tool is applied to generate a minimal set (in `example` folder) with a sample input file (exampleinput.fas). A *k*-mer window size of nine (9; nonamer) is used with utilising 14-cores. 

```
uniqmin -i exampleinput.fas -o example -k 9 -cpu 14
```

#### Command-line Arguments
| Argument 	| Parameter              | Type    	| Required | Default 	| Description                                |           
|----------	|----------------------- |---------	|----------|----------|------------------------------------------  |
| -h       	| help                   | N/A     	|FALSE	   | N/A     	| Show this help message and exit            |
| -i       	| sequence input file    | String  	|TRUE	     | N/A     	| Path of the input file (in FASTA format)   |
| -o       	| output directory name  | String  	|TRUE      | N/A     	| Path of the output file to be created      |
| -k        | *k*-mer window size    | Integer 	|FALSE     | 9       	| The length of *k*-mers to be used          |
| -cpu      | cpu size               | Integer 	|FALSE     | 14       | The number of CPU cores to be used         |

---
### Generate a random protein sequence dataset
This section is particular for the Protocol paper. For the details of this section and the python script, please refer to the [randomizer](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/randomizer) folder. 

---
### Citing Resources
* For original paper, please refer to our MDPI Biology paper: <br>
  [Chong, L.C.; Lim, W.L.; Ban, K.H.K.; Khan, A.M. An Alignment-Iindependent Approach for the Study of Viral Sequence Diversity at Any Given Rank of Taxonomy Lineage. *Biology* 2021, 10, 853. doi: 10.3390/biology10090853](https://www.mdpi.com/2079-7737/10/9/853) 
* For protocol paper, please refer to our preprint: <br> 
  [Chong, L.C.; Khan, A.M. UNIQmin, An Alignment-free Tool to Study Viral Sequence Diversity across Taxonomic Lineages: A Case Study of Monkeypox Virus. *bioRxiv* 2022.08.09.503271. doi: 10.1101/2022.08.09.503271](https://www.biorxiv.org/content/10.1101/2022.08.09.503271v2.full)
* For application paper to SARS-CoV-2, please refer to our preprint: <br> 
  [Chong, L.C.; Khan, A.M. Negligible Peptidome Diversity of SARS-CoV-2 and its Higher Taxonomic Ranks. *bioRxiv* 2022.10.31.513750. doi: 10.1101/2022.10.31.513750](https://www.biorxiv.org/content/10.1101/2022.10.31.513750v1.full)

---
## Found a bug?
Or would like a feature added? Or maybe drop some feedback?
Just [open a new issue](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/issues/new) or send an email to us (lichuinchong@gmail.com).
