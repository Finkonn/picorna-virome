import argparse
import pandas as pd
from Bio import SeqIO

def extract_and_process_proteins(genbank_file, fasta_file, protein_output_file, gapped_output_file):
    # Extract mature peptide coordinates from the GenBank file (nucleotide coordinates)
    records = SeqIO.parse(genbank_file, "genbank")
    mat_peptides = []

    for record in records:
        for feature in record.features:
            if feature.type == "mat_peptide":
                location = feature.location
                product = feature.qualifiers.get('product', [''])[0]
                mat_peptides.append({
                    'Start': location.start + 1,  # 1-based nucleotide start
                    'End': location.end,          # 1-based nucleotide end
                    'Product': product
                })

    protein_coordinates = pd.DataFrame(mat_peptides)
    protein_coordinates.to_csv(protein_output_file, index=False)
    print(f"Protein coordinates (nucleotide) saved to: {protein_output_file}")

    # Get the first sequence from the FASTA file (should be amino acid)
    def get_first_sequence(fasta_file):
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return str(record.seq)

    reference_sequence = get_first_sequence(fasta_file)
    print(f"Reference sequence length (amino acids): {len(reference_sequence)}")

    # Convert nucleotide coordinates to amino acid coordinates
    def nucleotide_to_aa_coords(nuc_start, nuc_end):
        """
        Convert nucleotide coordinates to amino acid coordinates.
        Assumes standard genetic code where 3 nucleotides = 1 amino acid.
        """
        # Convert to 0-based for calculation, then back to 1-based for output
        aa_start = ((nuc_start - 1) // 3) + 1
        aa_end = ((nuc_end - 1) // 3) + 1
        return aa_start, aa_end

    # Find gapped positions in the amino acid sequence
    def find_gapped_positions(sequence, ungapped_aa_coordinates):
        """
        Map ungapped amino acid coordinates to gapped coordinates in the sequence.
        """
        gapped_coordinates = []
        aa_position = 0  # Position in ungapped sequence (1-based)
        
        for i, letter in enumerate(sequence):
            if letter != "-":
                aa_position += 1
                if aa_position in ungapped_aa_coordinates:
                    gapped_coordinates.append(i + 1)  # 1-based gapped coordinate
        
        # Extract the subsequence from the gapped sequence
        if gapped_coordinates:
            seq = sequence[gapped_coordinates[0] - 1:gapped_coordinates[-1]]
        else:
            seq = ""
            
        return seq, gapped_coordinates

    # Process the table and find corresponding sequences with gapped positions
    def process_table(table, sequence):
        results = []
        for index, row in table.iterrows():
            nuc_start = row['Start']
            nuc_end = row['End']
            product = row['Product']
            
            # Convert nucleotide coordinates to amino acid coordinates
            aa_start, aa_end = nucleotide_to_aa_coords(nuc_start, nuc_end)
            
            # Create list of expected ungapped amino acid positions
            ungapped_aa_coordinates = list(range(aa_start, aa_end + 1))
            
            # Find these positions in the gapped sequence
            seq, gapped_coordinates = find_gapped_positions(sequence, ungapped_aa_coordinates)
            
            results.append({
                'Product': product,
                'Nucleotide_Start': nuc_start,
                'Nucleotide_End': nuc_end,
                'AA_Start': aa_start,
                'AA_End': aa_end,
                'Sequence': seq,
                'Gapped_Coordinates': str(gapped_coordinates)  # Convert list to string for CSV
            })
        
        results_df = pd.DataFrame(results)
        return results_df

    results = process_table(protein_coordinates, reference_sequence)
    results.to_csv(gapped_output_file, index=False)
    print(f"Gapped protein coordinates saved to: {gapped_output_file}")
    
    # Print summary
    print("\nSummary of extracted peptides:")
    for _, row in results.iterrows():
        print(f"  {row['Product']}: AA {row['AA_Start']}-{row['AA_End']} -> "
              f"gapped positions {row['Gapped_Coordinates']}")

def main():
    parser = argparse.ArgumentParser(description='Extract mat_peptide information from a GenBank file and map gapped coordinates to a reference sequence.')
    parser.add_argument('-g', '--genbank', required=True, help='Input GenBank file containing mat_peptide features (nucleotide coordinates).')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file with the reference sequence (amino acid sequence, may include gaps).')
    parser.add_argument('-p', '--protein_output', required=True, help='Output CSV file to save protein coordinates.')
    parser.add_argument('-o', '--gapped_output', required=True, help='Output CSV file to save gapped protein coordinates.')

    args = parser.parse_args()

    extract_and_process_proteins(args.genbank, args.fasta, args.protein_output, args.gapped_output)

if __name__ == '__main__':
    main()