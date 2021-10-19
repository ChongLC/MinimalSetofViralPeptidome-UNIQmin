# rpseqgen: A tool to generate random protein sequence dataset

Generation of the random protein sequence dataset may useful for many sequence analysis, in order to remove noise in the analysis. Herein, we would like to apply the same theory to UNIQmin and thus, this generator is built for the respective purpose. 

---

### Usage
`python randomGenerator.py [the path to the retrieved fasta sequence for amino acid composition] [the path to the output file] [sequence length for the generated dataset] [sequence number for the generated dataset]`

For example, rpseqgen tool is applied to generate a random protein sequence dataset `randomproteinseq.fasta` that consisting 1000 sequences that with 1000 amino acids long. The amino acid composition of the dataset generation is based on the retrieved dataset of all reported viral sequence (as of May 2021, from NCBI Entrez Protein Database) `allVirus080521.fasta`. <br> 

```
python randomGenerator.py allVirus080521.fasta randomproteinseq.fasta 1000 1000
```
