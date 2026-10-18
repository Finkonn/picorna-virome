import pandas as pd

input_file = "SRA_archives_ultra_final.xlsx"   

df = pd.read_excel(input_file)

paired = df[df["LibraryLayout"] == "PAIRED"]["Run"]
single = df[df["LibraryLayout"] == "SINGLE"]["Run"]

paired.to_csv("paired_runs_.txt", index=False, header=False)
single.to_csv("single_runs_.txt", index=False, header=False)

print(f"PAIRED: {len(paired)}")
print(f"SINGLE: {len(single)}")