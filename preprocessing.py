import numpy as np


query = {
    'C': 0, 'S': 1, 'T': 2, 'P': 3,
    'A': 4, 'G': 5, 'N': 6, 'D': 7,
    'E': 8, 'Q': 9, 'H': 10, 'R': 11,
    'K': 12, 'M': 13, 'I': 14, 'L': 15,
    'V': 16, 'F': 17, 'Y': 18, 'W': 19,
    'B': 20, 'X': 21
}


if __name__ == '__main__':
    protein_ids = []
    protein_seq = []
    with open("uniprot-reviewed_yes+AND+organism__yeast_.fasta") as f:
        line = f.readline()
        head = True
        seq = ''
        while line:
            if '>sp' in line:
                protein_id = line.split(' ')[0].split('|')[1]
                protein_ids.append(protein_id)
                if head:
                    head = False
                else:
                    protein_seq.append(seq)
                    seq = ''
            else:
                seq += line.strip()
            line = f.readline()
        protein_seq.append(seq)

    assert len(protein_ids) == len(protein_seq)
    length = len(protein_ids)
    with open('protein_ids.txt', 'w') as f:
        for i in range(length):
            f.write(f'{protein_ids[i]}\n')
    library = np.zeros(22 ** 3, dtype=str).tolist()
    flag = np.zeros(22 ** 3, dtype=bool)
    with open('sequences.txt', 'w') as f:
        for i in range(length):
            seq = protein_seq[i]
            leng = len(seq)
            flag.fill(False)
            for j in range(leng - 2):
                subseq = seq[j: j+3].upper()
                num = query[subseq[0]] * 22
                num = (num + query[subseq[1]]) * 22
                num = num + query[subseq[2]]
                if not flag[num]:
                    flag[num] = True
                    library[num] += str(i) + ","
            f.write(f'{seq}\n')
    with open('library3.txt', 'w') as f:
        for item in library:
            f.write(f'{item}\n')

    library = np.zeros(22 ** 4, dtype=str).tolist()
    flag = np.zeros(22 ** 4, dtype=bool)
    with open('sequences.txt', 'w') as f:
        for i in range(length):
            seq = protein_seq[i]
            leng = len(seq)
            flag.fill(False)
            for j in range(leng - 3):
                subseq = seq[j: j+4].upper()
                num = query[subseq[0]] * 22
                num = (num + query[subseq[1]]) * 22
                num = (num + query[subseq[2]]) * 22
                num = num + query[subseq[3]]
                if not flag[num]:
                    flag[num] = True
                    library[num] += str(i) + ","
            f.write(f'{seq}\n')
    with open('library4.txt', 'w') as f:
        for item in library:
            f.write(f'{item}\n')
