from Bio import SeqIO
import pandas as pd
import re

# Optional: load accessions if you want to filter GenBank records
df = pd.read_csv("scripts/bat_astro.csv")
accessions = set(df["ID"].astype(str).str.split('.').str[0])  # normalize IDs

results = []

for record in SeqIO.parse("data/Astroviridae_15102025.gb", "genbank"):
    # If you want to filter by accessions, uncomment:
    acc = record.id.split('.')[0]
    if acc not in accessions:
        continue

    refs = record.annotations.get("references", [])
    if not refs:
        continue

    for ref in refs:
        title = ref.title.strip() if ref.title else None
        journal = ref.journal.strip() if ref.journal else None
        pubmed = getattr(ref, "pubmed_id", None)

        if not title or not journal:
            continue

        # Skip "direct submission" or "unpublished"
        if "direct submission" in journal.lower() or "unpublished" in journal.lower():
            continue

        # Build link
        link = None
        if pubmed:
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed}/"
        else:
            doi_match = re.search(r'10\.\d{4,9}/\S+', journal)
            if doi_match:
                doi = doi_match.group(0).rstrip(".,;")
                link = f"https://doi.org/{doi}"

        if not link:
            continue

        results.append([title, journal, link])

# Build dataframe
out = pd.DataFrame(results, columns=["title", "journal", "link"])

# Remove duplicates by link
out = out.drop_duplicates(subset=["link"])

# Save
out.to_csv("genbank_papers_no_accessions.csv", index=False)