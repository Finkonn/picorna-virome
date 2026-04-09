import pandas as pd

input_file = "sra_filtered_no_16S.csv"   

df = pd.read_csv(input_file, sep=None, engine="python")

paired = df[df["LibraryLayout"] == "PAIRED"]["Run"]
single = df[df["LibraryLayout"] == "SINGLE"]["Run"]

paired.to_csv("paired_runs.txt", index=False, header=False)
single.to_csv("single_runs.txt", index=False, header=False)

print(f"PAIRED: {len(paired)}")
print(f"SINGLE: {len(single)}")