import numpy as np
import time


query = {
    'C': 0, 'S': 1, 'T': 2, 'P': 3,
    'A': 4, 'G': 5, 'N': 6, 'D': 7,
    'E': 8, 'Q': 9, 'H': 10, 'R': 11,
    'K': 12, 'M': 13, 'I': 14, 'L': 15,
    'V': 16, 'F': 17, 'Y': 18, 'W': 19,
    'B': 20, 'X': 21
}
blosum62 = [[9, -1, -1, -3, 0, -3, -3, -3, -4, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2, -3, -2],
            [-1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, -1, 0, -1, -2, -2, -2, -2, -2, -3, 0, 0],
            [-1, 1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, 0, -1, -2, -2, -2, -2, -2, -3, -1, 0],
            [-3, -1, 1, 7, -1, -2, -2, -1, -1, -1, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4, -2, -2],
            [0, 1, -1, -1, 4, 0, -2, -2, -1, -1, -2, -1, -1, -1, -1, -1, 0, -2, -2, -3, -2, 0],
            [-3, 0, 1, -2, 0, 6, 0, -1, -2, -2, -2, -2, -2, -3, -4, -4, -3, -3, -3, -2, -1, -1],
            [-3, 1, 0, -2, -2, 0, 6, 1, 0, 0, 1, 0, 0, -2, -3, -3, -3, -3, -2, -4, 3, -1],
            [-3, 0, 1, -1, -2, -1, 1, 6, 2, 0, -1, -2, -1, -3, -3, -4, -3, -3, -3, -4, 4, -1],
            [-4, 0, 0, -1, -1, -2, 0, 2, 5, 2, 0, 0, 1, -2, -3, -3, -2, -3, -2, -3, 1, -1],
            [-3, 0, 0, -1, -1, -2, 0, 0, 2, 5, 0, 1, 1, 0, -3, -2, -2, -3, -1, -2, 0, -1],
            [-3, -1, 0, -2, -2, -2, 1, -1, 0, 0, 8, 0, -1, -2, -3, -3, -3, -1, 2, -2, 0, -1],
            [-3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3, -1, -1],
            [-3, 0, 0, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5, -1, -3, -2, -2, -3, -2, -3, 0, -1],
            [-1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5, 1, 2, 1, 0, -1, -1, -3, -1],
            [-1, -2, -2, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 3, 0, -1, -3, -3, -1],
            [-1, -2, -2, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4, 1, 0, -1, -2, -4, -1],
            [-1, -2, -2, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4, -1, -1, -3, -3, -1],
            [-2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6, 3, 1, -3, -1],
            [-2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7, 2, -3, -1],
            [-2, -3, -3, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2, 11, -4, -2],
            [-3, 0, -1, -2, -2, -1, 3, 4, 1, 0, 0, -1, 0, -3, -3, -4, -3, -3, -3, -4, 4, -1],
            [-2, 0, 0, -2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1]]


def SM_alignment(seq1, seq2, threshold, score_matrix=blosum62):
    m = len(seq1)
    n = len(seq2)
    gap = -10
    dp = [[0]*(n + 1) for i in range(m + 1)]
    trace = [[0]*(n + 1) for i in range(m + 1)]
    trace_i = 0
    trace_j = 0
    max_score = 0
    temp = [0, 0, 0, 0]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            temp[0] = dp[i - 1][j] + gap
            temp[1] = dp[i][j - 1] + gap
            score = score_matrix[query[seq1[i - 1]]][query[seq2[j - 1]]]
            temp[2] = dp[i - 1][j - 1] + score
            ind = 3
            mmax = 0
            for k in range(3):
                if temp[k] > mmax:
                    mmax = temp[k]
                    ind = k
            dp[i][j] = mmax
            if mmax > max_score:
                max_score = mmax
                trace_i = i
                trace_j = j
            trace[i][j] = ind

    if max_score < threshold:
        return '', '', 0
    r, c = trace_i, trace_j

    align1 = ''
    align2 = ''
    while True:
        if trace[r][c] == 0:
            align1 += seq1[r - 1]
            align2 += '-'
            r -= 1
        elif trace[r][c] == 1:
            align1 += '-'
            align2 += seq2[c - 1]
            c -= 1
        elif trace[r][c] == 2:
            align1 += seq1[r - 1]
            align2 += seq2[c - 1]
            r -= 1
            c -= 1
        elif trace[r][c] == 3:
            break
        if (r == 0 and c == 0) or dp[r][c] == 0:
            break
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, max_score


def alignment(seq, threshold=30, output_num=5, seed_length=4, output_file=None):
    with open('./sequences.txt') as f:
        sequences = f.readlines()
    with open('./protein_ids.txt') as f:
        protein_ids = f.readlines()
    if seed_length == 3:
        library_path = 'library3.txt'
    elif seed_length == 4:
        library_path = 'library4.txt'
    else:
        raise Exception('Unexpected seed length!')
    with open(library_path) as f:
        library = f.readlines()

    aim = ''
    length = len(seq)
    sub_num = max(2, np.ceil(length / 2.5)) if seed_length == 3 else max(1, np.ceil(length / 18))
    flag = np.zeros(22 ** seed_length, dtype=bool)
    if seed_length == 3:
        for i in range(length - 2):
            subseq = seq[i: i + 3].upper()
            num = query[subseq[0]] * 22
            num = (num + query[subseq[1]]) * 22
            num = num + query[subseq[2]]
            if not flag[num]:
                flag[num] = True
                aim += (library[num].strip())
    else:
        for i in range(length - 3):
            subseq = seq[i: i + 4].upper()
            num = query[subseq[0]] * 22
            num = (num + query[subseq[1]]) * 22
            num = (num + query[subseq[2]]) * 22
            num = num + query[subseq[3]]
            if not flag[num]:
                flag[num] = True
                aim += (library[num].strip())
    aimed_seq = aim.split(',')[0: -1]
    aimed_seq.sort()
    leng = len(aimed_seq)
    final_aim = []
    count = 0
    last_seq = ''
    for i in range(leng - 1):
        if aimed_seq[i] == last_seq:
            count += 1
        else:
            if count >= sub_num:
                final_aim.append(int(last_seq))
            last_seq = aimed_seq[i]
            count = 0
    final_aim = list(set(final_aim))

    alignment1 = []
    alignment2 = []
    max_score = []
    names = []

    for i in final_aim:
        seq2 = sequences[i].strip()
        a1, a2, s = SM_alignment(seq2, seq, threshold)

        alignment1.append(a1), alignment2.append(a2), \
            max_score.append(s), names.append(protein_ids[i])

    alignment1 = np.array(alignment1)
    alignment2 = np.array(alignment2)
    max_score = np.array(max_score)
    names = np.array(names)

    index = max_score.argsort()
    alignment1 = alignment1[index[::-1]]
    alignment2 = alignment2[index[::-1]]
    max_score = max_score[index[::-1]]
    names = names[index[::-1]]
    output_num = min(output_num, len(names))
    if output_num == 0:
        print('No matching found')

    if output_file is not None:
        with open(output_file, 'w') as f:
            for i in range(output_num):
                f.write(f'protein name: {names[i].strip()}, Score: {max_score[i]}\n')
                f.write(f'input: {alignment2[i]}\n')
                f.write(f'align: {alignment1[i]}\n')
    else:
        for i in range(output_num):
            print(f'protein name: {names[i].strip()}, Score: {max_score[i]}')
            print(f'input: {alignment2[i]}')
            print(f'align: {alignment1[i]}')


if __name__ == '__main__':
    alignment('MLKLKFRWNKPDDGPRVKRQKSQQLTLDSPEGSETSSRPGGGYDNDRDVNGTFKRPWDITKIPGYNQTYVNEASLDKFAAYLKEEEVPVESSTPVSSEEMLPLEPVQKEYITSKHDWTPVFSVLKNNNSVNKRKMSHTALKSDFKHYLDMYDKEEKRRRGKSKKDKEGSEGDKTKTKEPKELVRKSVFSSILGVFSIFSKTKDGRKLKRDEIGEGIGFWLMRWPALFFIGMWLLFLTTIYASVRVLVAGYENILTWRGKRAKLRKVLRGARTYEQWVQAAQDLDVELGNAEWRENPKFGYYDHVTISKLTKMVRKLRLQDQAEDLSNILQGCIKNNFAGTESSTLYSQTYYGTKKVVEQWNEELGKAVTYVLESPKIDDEEKRDLFRLYSKNFGKSALCLSGGGCFAYLHFGIVKAMLDQDLLPQIISGTSGGALIAALACTRTDEELRQILVPELAYKITACWEPFPKWVFRWWRTGARFDSVDWARRSCWFTLGDMTFKEAYQRTGRILNVSTVPADPHSPVILCNYITSPDCLIWSALLASAAVPGILNPVMLMNKTKSGDIVPFSFGSKWKDGSLRTDIPVDALNTYFNVNCSIVSQVNPHIALFFYAPRGTVGRPVSHRKGKGWRGGFLGAALESMIKLEIRKWLKFIKAVELLPRFVDQDWTNVWLQRFSGSVTLWPKIHLADFWHILGDPWPEKMEDLLYRGQQCAFPKLLFLKHRMNIEKRIRNGRQATRHRKRKLEETDVCDLVRKTFSESDPDSDVKHSFVFPALMAPRSAGTDISSSNSDYDHEPQWEMDEGDSFLTADEEYPTTIL')
    '''
    for i in range(22):
        for j in range(i, 22):
            if blosum62[i, j] != blosum62[j, i]:
                print(f'i: {i},j: {j}')
    '''
