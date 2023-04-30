# **UNIQmin: An alignment-independent tool for the study of viral sequence diversity at any given rank of taxonomy lineage**

[![DOI - 10.3390/biology10090853](https://img.shields.io/badge/DOI-10.3390%2Fbiology10090853-2ea44f)](https://doi.org/10.3390/biology10090853)
![Python version](https://img.shields.io/pypi/pyversions/uniqmin)
[![PyPI](https://img.shields.io/pypi/v/uniqmin?logo=pypi)](https://pypi.org/project/uniqmin/)
[![GitHub tag](https://img.shields.io/github/tag/ChongLC/MinimalSetofViralPeptidome-UNIQmin)](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/releases/?include_prereleases&sort=semver "View GitHub releases")
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)

<!--
![Visitor](https://visitor-badge.laobi.icu/badge?page_id=https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin)
-->

### Brief Description
Sequence variation among viruses, even of a single amino acid, can expand their host repertoire or enhance the infection ability. Alignment-independent or -free approach represents an alternative to the study of viral diversity, which is devoid of the need for sequence conservation to perform comparative analyses. Herein, we present UNIQmin, a tool that utilises an alignment independent method to generate the minimal set of viral sequences, as a way to study their diversity, across any rank of taxonomic lineage. The minimal set refers to the smallest possible number of sequences required to capture the entire repertoire of viral peptidome diversity present in the given sequence dataset.

---
Table of Contents
====================
- [Step-by-step of UNIQmin]((#step-by-step-of-uniqmin))
- [Figure Scheme](#figure-scheme)
- [UNIQmin as a Package](#uniqmin-as-a-package)
- [Citing Resources](#citing-resources)
- [Found a Bug](#found-a-bug)

## Step-by-step of UNIQmin
UNIQmin comprises of five execution steps, with a Python script for each step. These scripts are provided in the [PythonScript](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/tree/master/PythonScript) folder, accompanied with detailed explanations for the algorithm and execution of each step. The sample input file (`exampleinput.fas`) and example output file (`exampleoutput.fasta`) are also provided. 

## Figure Scheme
![uniqminScheme](https://user-images.githubusercontent.com/51225708/222265004-8ea85997-f289-42bf-9ac6-74d1beae2ee6.png)

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
  
  **Note for users who use Conda environment (*e.g.*: via Jupyter Notebook):** <br>
  Before `pip` install of the package, run <br>
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

### Usage
`uniqmin [-h] [-i INPUT] [-o OUTPUT] [-k KMERLENGTH] [-t THREADS]`

A sample usage: <br>
The UNIQmin tool is applied to generate a minimal set (to be saved in an output folder, named "result") for a sample input file (named "exampleinput.fas") with a *k*-mer window size parameter of nine (9; nonamer) and utilising 4-threads of a CPU (subject to limitations of the resource used): 

```
uniqmin -i exampleinput.fas -o result -k 9 -t 4
```

#### Command-line Arguments
| Argument 	| Parameter              | Type    	| Required | Default 	| Description                                |           
|----------	|----------------------- |---------	|----------|----------|------------------------------------------  |
| -h       	| help                   | N/A     	|FALSE	   | N/A     	| Show this help message and exit            |
| -i       	| sequence input file    | String  	|TRUE	     | N/A     	| Path of the input file (in FASTA format)   |
| -o       	| output directory name  | String  	|TRUE      | N/A     	| Directory to store the output file to be created      |
| -k        | *k*-mer window size    | Integer 	|FALSE     | 9       	| The length of the constituent *k*-mer to be used          |
| -t        | number of threads      | Integer 	|FALSE     | 4        | The number of CPU threads to be used           |

---
### Citing Resources
* For original reference, please refer to our paper: <br>
  [Chong, L.C.; Lim, W.L.; Ban, K.H.K.; Khan, A.M. An Alignment-independent Approach for the Study of Viral Sequence Diversity at Any Given Rank of Taxonomy Lineage. *Biology* 2021, 10, 853. doi: 10.3390/biology10090853](https://www.mdpi.com/2079-7737/10/9/853) 
* For a protocol describing the step-by-step utility of UNIQmin, please refer to our preprint: <br> 
  [Chong, L.C.; Khan, A.M. UNIQmin, An Alignment-free Tool to Study Viral Sequence Diversity across Taxonomic Lineages: A Case Study of Monkeypox Virus. *bioRxiv* 2022.08.09.503271. doi: 10.1101/2022.08.09.503271](https://www.biorxiv.org/content/10.1101/2022.08.09.503271v2.full)
* For application of UNIQmin to study SARS-CoV-2 lineage diversity, please refer to our preprint: <br> 
  [Chong, L.C.; Khan, A.M. Negligible Peptidome Diversity of SARS-CoV-2 and its Higher Taxonomic Ranks. *bioRxiv* 2022.10.31.513750. doi: 10.1101/2022.10.31.513750](https://www.biorxiv.org/content/10.1101/2022.10.31.513750v1.full)

---
## Found a bug?
Or would you like a feature to be added? Or maybe drop some feedback?
Just [open a new issue](https://github.com/ChongLC/MinimalSetofViralPeptidome-UNIQmin/issues/new) or send an email to us (lichuinchong@gmail.com).
