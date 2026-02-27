from Bio import SeqIO

records = SeqIO.parse("data/Astroviridae_filtered.gb", "genbank")

modified_records = []

for record in records:
    accession = record.annotations["accessions"][0]

    record.id = accession
    record.name = accession
    record.description = accession + " " + record.description.split(" ", 1)[1] if " " in record.description else accession

    modified_records.append(record)

count = SeqIO.write(modified_records, "data/Astroviridae.fasta", "fasta")

print(f"Converted {count} records")