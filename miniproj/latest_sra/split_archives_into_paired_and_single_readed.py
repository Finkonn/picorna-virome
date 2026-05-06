import pandas as pd

input_file = "our_SRA_with_astro_from_ZOVER.csv"   

df = pd.read_csv(input_file, sep=None, engine="python")

paired = df[df["LibraryLayout"] == "PAIRED"]["Run"]
single = df[df["LibraryLayout"] == "SINGLE"]["Run"]

paired.to_csv("paired_runs_astro_ZOVER.txt", index=False, header=False)
single.to_csv("single_runs_astro_ZOVER.txt", index=False, header=False)

print(f"PAIRED: {len(paired)}")
print(f"SINGLE: {len(single)}")