import matplotlib.pyplot as plt
from Bio import SeqIO
import re
import argparse

def plot_ambiguous_segment_lengths(fasta_file):
    ambiguous_bases = set('RYSWKMBDHVN')
    longest_lengths = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        segments = re.findall(f'[{"".join(ambiguous_bases)}]+', seq)
        if segments:
            longest_segment = max(segments, key=len)
            longest_lengths.append(len(longest_segment))
        else:
            longest_lengths.append(0)
    
    plt.hist(longest_lengths, bins=range(1, max(longest_lengths) + 2), edgecolor='black')
    plt.xlabel('Length of longest ambiguous segment', fontsize=24)
    plt.ylabel('Number of sequences', fontsize=24)
    plt.title('Distribution of longest ambiguous segment lengths in sequences', fontsize=26)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot distribution of longest ambiguous segment lengths in sequences')
    parser.add_argument('-input', required=True, help='Path to the input FASTA file')
    args = parser.parse_args()
    
    plot_ambiguous_segment_lengths(args.input)
