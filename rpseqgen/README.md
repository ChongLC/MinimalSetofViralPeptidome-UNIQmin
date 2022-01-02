# rpseqgen: A tool to generate random protein sequence dataset
![Run in Google Colab](https://img.shields.io/badge/Colab-Run_in_Google_Colab-blue?logo=Google&logoColor=FDBA18) <br>

Generation of the random protein sequence dataset may useful for many sequence analysis, in order to remove noise in the analysis. Herein, we would like to apply the same theory to UNIQmin and thus, this generator is built for the respective purpose. 

---

### Usage
`python randomGenerator.py [-h] [-[-o OUTPUT] [-l SEQLEN] [-n SEQNUM]`

For example, rpseqgen tool is applied to generate a random protein sequence dataset `randomproteinseq.fasta` consisting of 1,000 sequences that with 1,000 amino acids long. The provided amino acid composition of the dataset generation is based on the retrieved dataset of all reported viral sequence (as of May 2021, from NCBI Entrez Protein Database) `allVirus080521.fasta`. <br> 

```
python randomGenerator.py -o randomproteinseq.fasta -l 1000 -n 1000
```

#### Command-line Arguments
| Argument | Parameter | Type    	| Required | Description                                |           
|:--------:|-----------|---------	|:--------:|------------------------------------------  |
| -h       | help      | N/A     	|FALSE	   | Show this help message and exit            |
| -o       | output    | String  	|TRUE      | Path of the output file to be created      |
| -l       | seqlen    | Integer 	|TRUE      | The length of protein sequences            |
| -n       | seqnum    | Integer 	|TRUE      | The number of protein sequences            |
