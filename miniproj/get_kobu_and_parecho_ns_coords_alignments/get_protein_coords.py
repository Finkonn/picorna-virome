import argparse
import pandas as pd
from Bio import SeqIO

def extract_and_process_proteins(genbank_file, fasta_file, protein_output_file, gapped_output_file):
    # Extract mature peptide coordinates from the GenBank file
    records = SeqIO.parse(genbank_file, "genbank")
    mat_peptides = []

    for record in records:
        for feature in record.features:
            if feature.type == "mat_peptide":
                location = feature.location
                product = feature.qualifiers.get('product', [''])[0]
                mat_peptides.append({
                    'Start': location.start + 1,
                    'End': location.end,
                    'Product': product
                })

    protein_coordinates = pd.DataFrame(mat_peptides)
    protein_coordinates.to_csv(protein_output_file, index=False)
    print(f"Protein coordinates saved to: {protein_output_file}")

    # Get the first sequence from the FASTA file
    def get_first_sequence(fasta_file):
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)

    reference_sequence = get_first_sequence(fasta_file)

    # Find gapped positions
    def find_gapped_positions(sequence, ungapped_coordinates):
        gapped_coordinates = []
        gapped_index = 0
        for i, letter in enumerate(sequence):
            if letter != "-":
                gapped_index += 1
                if gapped_index in ungapped_coordinates:
                    gapped_coordinates.append(i + 1)

        seq = sequence[gapped_coordinates[0] - 1:gapped_coordinates[1]]
        return seq, gapped_coordinates

    # Process the table and find corresponding sequences with gapped positions
    def process_table(table, sequence):
        results = []
        for index, row in table.iterrows():
            start = row['Start']
            end = row['End']
            product = row['Product']
            ungapped_coordinates = list(range(start, end + 1))
            seq, gapped_coordinates = find_gapped_positions(sequence, ungapped_coordinates)
            results.append({
                'Product': product,
                'Sequence': seq,
                'Gapped_Coordinates': gapped_coordinates
            })
        results_df = pd.DataFrame(results)
        return results_df

    results = process_table(protein_coordinates, reference_sequence)
    results.to_csv(gapped_output_file, index=False)
    print(f"Gapped protein coordinates saved to: {gapped_output_file}")

def main():
    parser = argparse.ArgumentParser(description='Extract mat_peptide information from a GenBank file and map gapped coordinates to a reference sequence.')
    parser.add_argument('-g', '--genbank', required=True, help='Input GenBank file containing mat_peptide features.')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file with the reference sequence.')
    parser.add_argument('-p', '--protein_output', required=True, help='Output CSV file to save protein coordinates.')
    parser.add_argument('-o', '--gapped_output', required=True, help='Output CSV file to save gapped protein coordinates.')

    args = parser.parse_args()

    extract_and_process_proteins(args.genbank, args.fasta, args.protein_output, args.gapped_output)

if __name__ == '__main__':
    main()
