# randpSeqGen: A tool to generate a random protein sequence dataset
[![Run in Google Colab](https://img.shields.io/badge/Colab-Run_in_Google_Colab-blue?logo=Google&logoColor=FDBA18)](https://colab.research.google.com/drive/1IwNPKaRKGgPzqiOBuEo8S0VbUpe3XqVh?usp=sharing) <br>

A random protein sequence dataset can be useful for various sequence analysis, in order to evaluate and correct the analysis for the background noise. Herein, we offer a tool that can generate a dataset of random protein sequences. 

---

### Usage
`python randpseqgen.py [-h] [-[-o OUTPUT] [-l SEQLEN] [-n SEQNUM]`

In the usage case below, the randpSeqGen tool is applied to generate a random protein sequence dataset, named `randomproteinseq.fasta` consisting of 1,000 sequences of length 1,000 amino acids. The amino acid composition of the random sequences is based on all reported viral sequence retrieved from the NCBI Protein database (as of May 2021; `allVirus080521.fasta`). <br> 

```
python randpseqgen.py -o randomproteinseq.fasta -l 1000 -n 1000
```

#### Command-line Arguments
| Argument | Parameter | Type    	| Required | Description                                             |           
|:--------:|-----------|---------	|:--------:|-------------------------------------------------------  |
| -h       | help      | N/A     	|FALSE	   | Show this help message and exit                         |
| -o       | output    | String  	|TRUE      | Path of the output file to be created                   |
| -l       | seqlen    | Integer 	|TRUE      | The length of random protein sequences to be generated  |
| -n       | seqnum    | Integer 	|TRUE      | The number of random protein sequences to be generated  |
