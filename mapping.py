import os
import argparse
from translate import translate
from SW import alignment


parser = argparse.ArgumentParser(description='Yeast protein sequences matching tool')
parser.add_argument('--inputfile', default=None,
                    type=str, help='The path of the file that stores the DNA sequence to be matched')
parser.add_argument('--DNA', default=None,
                    type=str, help='The DNA sequence to be matched')
parser.add_argument('--seq_type', default='positive',
                    type=str, help='The input is positive-sense or negative-sense, default is positive')
parser.add_argument('--start_codon', default=True,
                    type=bool, help='Whether to start from the start codon')
parser.add_argument('--threshold', default=30,
                    type=int, help='Minimum alignment score to select protein sequences')
parser.add_argument('--output_num', default=5,
                    type=int, help='The maximum number of output protein sequences')
parser.add_argument('--seed_len', default=4,
                    type=int, help='The length of seed')
parser.add_argument('--output_dir', default=None,
                    type=str, help='Where to store the results, default is printing it')
args = parser.parse_args()


if __name__ == '__main__':
    if args.inputfile is None and args.DNA is None:
        raise Exception('Please enter the sequence(either by entering the file name(--input_file) or directly(--DNA))')
    if args.inputfile is not None:
        seq = ''
        with open(args.inputfile) as f:
            line = f.readline().strip()
            seq = seq + line
            while line:
                line = f.readline().strip()
                seq = seq + line
    else:
        seq = args.DNA

    protein_seq = translate(seq, args.seq_type, args.start_codon)
    alignment(protein_seq, args.threshold, args.output_num, args.seed_len, args.output_dir)
