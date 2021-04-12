# UNIQmin: An alignment-independent tool for the study of pathogen sequence diversity at any given rank of taxonomy lineage

### Brief Description
Sequence variation among pathogens, even of a single amino acid, can expand their host repertoire or enhance the infection ability. Alignment independent approach represents an alternative approach to the study of pathogen diversity, which is devoid of the need for sequence conservation to perform comparative analyses. Herein, we present UNIQmin, a tool that utilises an alignment independent method to generate the minimal set of pathogen sequences, as a way to study their diversity, across any rank of taxonomic lineage. The minimal set refers to the smallest possible number of sequences required to capture the entire repertoire of pathogen peptidome diversity present in a sequence dataset.

### Installation

```
pip install uniqmin
```

### Usage

Open terminal and type uniqmin with the input parameters follow the sequence as follows: file_name, kmer_size, cpu_core

```
uniqmin <file_name> <kmer_size> <cpu_core>
```

Example:
```
uniqmin input.fasta 9 6
```

---
## Citing Resources
