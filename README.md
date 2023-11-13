### Yeast whole proteome query program 

#### 1. File description
- [preprocessing.py](preprocessing.py) pre-processing yeast protein sequences file from UniProt
- [translate.py](translate.py) translate the DNA sequence to protein sequence
- [mapping.py](mapping.py) mapping the sequence to the yeast whole proteome
- [SW.py](SW.py) implement the Smith-Waterman algorithm
- [protein_ids.txt](protein_ids.txt) stores preprocessed yeast protein id information
- [sequences.txt](sequences.txt) stores preprocessed yeast protein sequence information
- [uniprot-reviewed_yes+AND+organism__yeast_.fasta](uniprot-reviewed_yes+AND+organism__yeast_.fasta) is the yeast protein data from UniProt

#### 2. Install the dependency:
```bash
pip install -r requirements.txt
```

#### 3. Pre-processing the data:
```bash
python preprocessing.py
```
It would take 1~2 minutes.

#### 4. Alignment:
To align the sequence:
```bash
python mapping.py --DNA [input sequence]
```
or
```bash
python mapping.py --inputfile [The path storing the DNA sequence]
```
To see more detailed command line arguments, using
```bash
python mapping.py -h
```
