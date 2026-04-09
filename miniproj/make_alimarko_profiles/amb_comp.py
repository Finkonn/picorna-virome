import argparse
from Bio import SeqIO
from matplotlib import pyplot as plt

'''
Plots distribution of ambiguous nucleotides in sequences

Input:
input_file - file with sequences in fasta format

Output:

distribution plot for every type of ambiguous nucleotide
'$input_file'_res.fasta - file with sequences that have no ambiguous nucleotides
'''

def search_amb(input_file):
    # List of ambiguous nucleotides
    ambig_nt = [
        'n',  # A, T, G, C 
        'r', # Purine (A or G)
        'y', # Pyrimidine (T or C)
        'k', # Keto (G or T)
        'm', # Amino (A or C)
        's', # Strong interaction (3 H bonds) (G or C)
        'w', # Weak interaction (2 H bonds) (A or T)
        'b', # Not A (C or G or T)
        'd', # Not C (A or G or T)
        'h', # Not G (A or C or T)
        'v', # Not T or U (A or C or G)
        'total'
    ]
    
    # Counts of each type of ambiguous nucleotide
    amb_counts = {nt: [] for nt in ambig_nt}

    # Checks presence of ambiguous nucleotides in sequences from input_file
    new_records = []
    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        print(f"Total sequences: {len(records)}")
        for record in records.copy():
            # Total number of ambiguous nucleotides
            total_amb = 0
            for nt in ambig_nt[:-1]:  # Exclude 'total' from counting
                count = record.seq.lower().count(nt)
                amb_counts[nt].append(count)
                total_amb += count
            amb_counts['total'].append(total_amb)
            
            # Filter sequences with ambiguous nucleotides
            if total_amb < 1:
                new_records.append(record)
            else:
                print(f"{record.id} has {total_amb} ambiguous nucleotides")
    
    # Write sequences without ambiguous nucleotides to a new file
    #output_file = input_file.replace('.fasta', '_res.fasta')
    #SeqIO.write(new_records, output_file, format="fasta")
    #print(f"Sequences without ambiguous nucleotides written to {output_file}")
    
    # Plots distribution histograms
    for key in ambig_nt:
        if len(amb_counts[key]) != 0:
            plt.hist(amb_counts[key], bins=list(range(max(amb_counts[key]) + 5)), log=True)
            plt.title(f'Distribution of {key}', fontsize=26)
            plt.xlabel('Number of ambiguous characters in sequence', fontsize=24)
            plt.ylabel('Number of sequences', fontsize=24)
            plt.tick_params(axis='both', which='major', labelsize=16)
            plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot distribution of ambiguous nucleotides in sequences')
    parser.add_argument('-i', '--input_file', help='Input file in fasta format containing sequences')

    args = parser.parse_args()

    search_amb(args.input_file)

if __name__ == '__main__':
    main()
