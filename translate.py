codon = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
    "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def transcript(seq):
    pair = str.maketrans('TCGA', 'UGCT')
    seq = seq[::-1].translate(pair)
    return seq


def translate(seq, seq_type, start_codon):
    seq = seq.upper()
    if seq_type == 'negative':
        seq = transcript(seq)
    elif seq_type == 'positive':
        seq = seq.replace('T', 'U')
    else:
        raise Exception('please input --seq_type \'positive\' or \'negative\'!')
    start = 0 if start_codon is False else seq.find('AUG')
    if start == -1:
        raise Exception('start codon not found!')
    protein = ''
    length = len(seq)
    for i in range(start, length, 3):
        if i + 3 > length:
            break
        sub_seq = seq[i: i+3]
        temp = codon[sub_seq]
        if temp == 'STOP':
            break
        protein += temp
    if len(protein) <= 5:
        raise Exception('translated protein is too short!')
    return protein


if __name__ == '__main__':
    print(translate('AATGCCUGATCG', 'negative', True))
