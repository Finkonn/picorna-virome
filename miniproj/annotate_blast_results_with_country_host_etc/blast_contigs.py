#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO

TARGET_TAXA = {"Astroviridae", "Parechovirus", "Avastrovirus", "Mamastrovirus", "Kobuvirus", "Sakobuvirus", "Salivirus"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", required=True)
    parser.add_argument("-db", "--diamond_db", required=True)
    parser.add_argument("--mode", choices=["target", "best"], default="target")
    parser.add_argument("-t", "--threads", default=8, type=int)
    parser.add_argument("-o", "--output", default="kobu_pipeline_output")
    return parser.parse_args()


def best_hits(df):
    df = df.sort_values("Score_ratio", ascending=False)
    df = df.drop_duplicates(subset=["Name"])
    return df


def extract_contigs(sample_name, df, fasta_path):
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    extracted = []

    for _, row in df.iterrows():
        contig = row["Name"]
        if contig not in seqs:
            continue

        rec = seqs[contig]

        taxon = row["Taxon"]
        contig_len = int(row["Length_contig"])
        from_pos = int(row["From"])
        to_pos = int(row["To"])
        score_ratio = row["Score_ratio"]

        new_id = f"{sample_name}|{contig}/{taxon}/{contig_len}/{from_pos}-{to_pos}/{score_ratio}"
        rec.id = new_id
        rec.description = ""
        extracted.append(rec)

    return extracted


def run_diamond(fasta, db, threads, outdir):
    out_file = os.path.join(outdir, "diamond_results.tsv")

    cmd = [
        "diamond", "blastx",
        "--ultra-sensitive",
        "-d", db,
        "-q", fasta,
        "-o", out_file,
        "-f", "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "qcovhsp",
        "scovhsp",
        "stitle",
        "-k", "5",
        "-p", str(threads)
    ]

    subprocess.run(cmd, check=True)
    return out_file

def summarize_diamond(diamond_file, outdir):
    df = pd.read_csv(diamond_file, sep="\t", header=None)

    df.columns = [
        "qseqid",
        "sseqid",
        "pident",
        "align_len",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "qcovhsp",
        "scovhsp",
        "stitle"
    ]

    df = df.sort_values("bitscore", ascending=False)

    best = df.drop_duplicates("qseqid")

    summary_path = os.path.join(outdir, "diamond_best_hits.tsv")
    best.to_csv(summary_path, sep="\t", index=False)

    return summary_path



def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)

    all_records = []

    for file in os.listdir(args.directory):
        if not file.endswith(".csv"):
            continue

        sample = file.replace(".csv", "")
        csv_path = os.path.join(args.directory, file)
        fasta_path = os.path.join(
            args.directory,
            f"{sample}_matched_contigs.fasta"
        )

        if not os.path.exists(fasta_path):
            continue

        df = pd.read_csv(csv_path)
        
        if args.mode == "target":
            df = df.sort_values(["Name", "Score_ratio"], ascending=[True, False])

            selected_rows = []

            for contig, group in df.groupby("Name"):
                group = group.reset_index(drop=True)

                if group.empty:
                    continue

                top = group.iloc[0]

                if top["Taxon"] in TARGET_TAXA:
                    selected_rows.append(top)

        elif args.mode == "best":
            df = df.sort_values(["Name", "Score_ratio"], ascending=[True, False])

            selected_rows = []

            for contig, group in df.groupby("Name"):
                group = group.reset_index(drop=True)

                if group.empty:
                    continue

                selected = None

                for _, row in group.iterrows():
                    if row["Score_ratio"] > 1 and row["Taxon"] != "Picornaviridae":
                        selected = row
                        break

                if selected is not None:
                    selected_rows.append(selected)

            df = pd.DataFrame(selected_rows)

            if df.empty:
                continue
        
        elif args.mode == "best":
            df = df.sort_values(["Name", "Score_ratio"], ascending=[True, False])
            df = df.drop_duplicates(subset=["Name"])
            df = df[df["Score_ratio"] > 1]

        records = extract_contigs(sample, df, fasta_path)
        all_records.extend(records)

    if not all_records:
        print("No contigs found.")
        return

    merged_fasta = os.path.join(args.output, "ALL_kobu_candidates.fasta")
    SeqIO.write(all_records, merged_fasta, "fasta")

    print(f"Total contigs extracted: {len(all_records)}")

    diamond_out = run_diamond(
        merged_fasta,
        args.diamond_db,
        args.threads,
        args.output
    )

    summary = summarize_diamond(diamond_out, args.output)

    print("Pipeline finished.")
    print(f"Summary: {summary}")


if __name__ == "__main__":
    main()
