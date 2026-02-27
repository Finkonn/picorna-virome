from Bio import SeqIO
from Bio.Seq import Seq
import os

fasta_files = [
    "Astroviridae_1A_cleaned.fasta",
    "Astroviridae_1B_cleaned.fasta",
    "Astroviridae_2_cleaned.fasta"
]

output_dir = "translated_results"

os.makedirs(output_dir, exist_ok=True)

def count_internal_stops(protein_seq):
    stop_count = protein_seq.count("*")
    if protein_seq.endswith("*"):
        stop_count -= 1
    return stop_count

def translate_best_frame(nuc_seq):
    for frame in [0, 1, 2]:
        trimmed_len = len(nuc_seq) - frame
        trimmed_len -= trimmed_len % 3
        sub_seq = nuc_seq[frame:frame + trimmed_len]

        protein = sub_seq.translate()
        protein_str = str(protein)
        if count_internal_stops(str(protein)) == 0:
            if protein_str.endswith('*'):
                protein_str = protein_str[:-1]
            return protein_str
    return None

for fasta_file in fasta_files:
    output_fasta = os.path.join(output_dir, f"{os.path.splitext(fasta_file)[0]}_aa.fasta")
    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            nuc_seq = record.seq
            protein = translate_best_frame(nuc_seq)

            if protein:
                out_handle.write(f">{record.id}\n")
                for i in range(0, len(protein), 60):
                    out_handle.write(protein[i:i+60] + "\n")
            else:
                print(f"Skipping {record.id} from {fasta_file}: no frame without internal stop codons.")