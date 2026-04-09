import argparse
from Bio import SeqIO
import os


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

        if count_internal_stops(protein_str) == 0:
            if protein_str.endswith('*'):
                protein_str = protein_str[:-1]
            return protein_str

    return None


def main():
    parser = argparse.ArgumentParser(description="Translate nucleotide FASTA into protein sequences (best frame, no internal stops).")
    parser.add_argument("-i", "--input", required=True, help="Input nucleotide FASTA")
    parser.add_argument("-o", "--output", help="Output protein FASTA (optional)")

    args = parser.parse_args()

    input_fasta = args.input

    if args.output:
        output_fasta = args.output
    else:
        base = os.path.splitext(os.path.basename(input_fasta))[0]
        output_fasta = f"{base}_aa.fasta"

    with open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            protein = translate_best_frame(record.seq)

            if protein:
                out_handle.write(f">{record.id}\n")
                for i in range(0, len(protein), 60):
                    out_handle.write(protein[i:i+60] + "\n")
            else:
                print(f"Skipping {record.id}: no frame without internal stop codons.")

    print(f"\nDone. Output written to: {output_fasta}")


if __name__ == "__main__":
    main()