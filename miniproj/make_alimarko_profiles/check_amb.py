from Bio import SeqIO
import re

input_fasta = "Astroviridae_1A_cleaned.fasta"  

seq_with_N = 0
seq_with_2plus = 0
seq_with_3plus = 0

total_2plus_runs = 0
total_3plus_runs = 0

total_sequences = 0

for record in SeqIO.parse(input_fasta, "fasta"):
    total_sequences += 1
    seq = str(record.seq).upper()

    # Any N
    if "N" in seq:
        seq_with_N += 1

    # ≥2 consecutive Ns
    runs_2plus = re.findall(r"N{2,}", seq)
    if runs_2plus:
        seq_with_2plus += 1
        total_2plus_runs += len(runs_2plus)

    # ≥3 consecutive Ns
    runs_3plus = re.findall(r"N{3,}", seq)
    if runs_3plus:
        seq_with_3plus += 1
        total_3plus_runs += len(runs_3plus)

print("===== SUMMARY =====")
print("Total sequences:", total_sequences)
print("Sequences with at least one N:", seq_with_N)
print("Sequences with ≥2 consecutive Ns:", seq_with_2plus)
print("Sequences with ≥3 consecutive Ns:", seq_with_3plus)
print()
print("Total ≥2-N runs:", total_2plus_runs)
print("Total ≥3-N runs:", total_3plus_runs)