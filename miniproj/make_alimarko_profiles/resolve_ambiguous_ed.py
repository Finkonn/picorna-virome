import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from Bio import SeqIO

ambig_nt = set('RYWSKMBDHVNrywskmbdhvn')


# --- STRICT ACCESSION PARSER ---
def normalize_id(rec, seen_ids):
    raw = rec.id.split()[0]
    acc = raw.split('/')[0]  # <-- YOUR RULE

    # ensure uniqueness
    if acc in seen_ids:
        seen_ids[acc] += 1
        acc = f"{acc}_{seen_ids[acc]}"
    else:
        seen_ids[acc] = 0

    rec.id = acc
    rec.description = ''
    return rec


def resolve_ambiguous(input_file, window, path_to_blast,
                      evalue, word_size, max_ambiguous, max_ambiguous_row):

    input_dir = os.path.dirname(os.path.realpath(input_file))
    os.makedirs(input_dir, exist_ok=True)

    # --- LOAD + NORMALIZE IDS ---
    fasta_raw = list(SeqIO.parse(open(input_file), "fasta"))

    seen_ids = {}
    fasta_al = [normalize_id(rec, seen_ids) for rec in fasta_raw]

    fasta_al_less_amb = []
    list_slices = []
    slice_map = {}

    pattern = r"[RYWSKMBDHVN]{" + str(max_ambiguous_row) + ",}"

    slice_counter = 0

    # --- FILTER + SLICE ---
    for rec in fasta_al:
        amb_total = len(re.findall(r"[RYSWKMBDHVN]", str(rec.seq)))

        if amb_total == 0:
            fasta_al_less_amb.append(rec)
            continue

        if amb_total > max_ambiguous:
            continue

        if re.findall(pattern, str(rec.seq)):
            continue

        fasta_al_less_amb.append(rec)

        starts = [m.start() for m in re.finditer(r"[RYSWKMBDHVN]", str(rec.seq))]
        i = 0

        while i < len(starts):
            current_starts = [starts[i]]
            st = max(0, starts[i] - window // 2)
            e = min(len(rec.seq), starts[i] + window // 2)

            cur_slice_rec = copy.deepcopy(rec[st:e])

            slice_id = f"s{slice_counter}"
            slice_counter += 1

            cur_slice_rec.id = slice_id
            cur_slice_rec.description = ''
            list_slices.append(cur_slice_rec)

            for j in range(i + 1, len(starts)):
                if st + window // 5 < starts[j] < e - window // 5:
                    current_starts.append(starts[j])
                    i += 1
                else:
                    break

            slice_map[slice_id] = {
                "orig_id": rec.id,
                "start": st,
                "amb_pos": current_starts
            }

            i += 1

    # --- WRITE FASTA FILES ---
    file_name_slices = os.path.join(input_dir, "slices.fasta")
    file_name_less_amb = os.path.join(input_dir, "less_amb.fasta")

    SeqIO.write(list_slices, file_name_slices, "fasta")
    SeqIO.write(fasta_al_less_amb, file_name_less_amb, "fasta")

    # --- BLAST ---
    makeblast = os.path.join(path_to_blast,
                             "makeblastdb.exe" if sys.platform.startswith("win") else "makeblastdb")
    blastn = os.path.join(path_to_blast,
                          "blastn.exe" if sys.platform.startswith("win") else "blastn")

    db_path = os.path.join(input_dir, "local_db")
    blast_out = os.path.join(input_dir, "blast.out")

    makeblast_command = f'"{makeblast}" -in "{file_name_less_amb}" -dbtype nucl -parse_seqids -out "{db_path}"'
    blastn_command = f'"{blastn}" -db "{db_path}" -query "{file_name_slices}" -outfmt 6 -out "{blast_out}" -strand plus -evalue {evalue} -word_size {word_size} -max_target_seqs 30'

    subprocess.run(makeblast_command, shell=True, check=True)
    subprocess.run(blastn_command, shell=True, check=True)

    # --- LOAD BLAST ---
    blast_output = pd.read_csv(
        blast_out,
        sep='\t',
        header=None,
        names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    )

    fasta_dict = SeqIO.to_dict(fasta_al_less_amb)

    # --- RESOLVE ---
    for _, row in blast_output.iterrows():

        qid = row['qseqid']
        sid = row['sseqid']

        if qid not in slice_map:
            continue
        if sid not in fasta_dict:
            continue

        meta = slice_map[qid]
        orig_id = meta["orig_id"]

        if orig_id not in fasta_dict:
            continue
        if sid == orig_id:
            continue

        start = meta["start"]
        amb_positions = meta["amb_pos"]

        cur_start = start + int(row['qstart']) - 1
        rel_amb_pos = [x - cur_start for x in amb_positions]

        ref_pos = [int(row['sstart']) - 1 + x for x in rel_amb_pos]

        ref_seq = fasta_dict[sid].seq

        for i in range(len(amb_positions)):
            if i >= len(ref_pos):
                continue

            pos = ref_pos[i]
            if pos >= len(ref_seq):
                continue

            nuc = ref_seq[pos]

            if nuc not in ambig_nt:
                target_seq = fasta_dict[orig_id].seq
                amb_index = amb_positions[i]

                fasta_dict[orig_id].seq = (
                    target_seq[:amb_index] +
                    nuc +
                    target_seq[amb_index + 1:]
                )

    # --- FINAL ---
    cleaned_seqs = []
    removed_seqs = []

    for rec in fasta_dict.values():
        if any(nt in ambig_nt for nt in str(rec.seq)):
            removed_seqs.append(rec)
        else:
            cleaned_seqs.append(rec)

    SeqIO.write(cleaned_seqs, os.path.join(input_dir, "cleaned.fasta"), "fasta")
    SeqIO.write(removed_seqs, os.path.join(input_dir, "removed.fasta"), "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str, required=True)
    parser.add_argument("-w", "--window", type=int, default=100)
    parser.add_argument("-evalue", "--evalue", type=float, default=1e-20)
    parser.add_argument("-word_size", "--word_size", type=int, default=7)
    parser.add_argument("-pb", "--path_blast", type=str, required=True)
    parser.add_argument("-max_ambiguous", "--max_ambiguous", type=float, default=5)
    parser.add_argument("-max_ambiguous_row", "--max_ambiguous_row", type=int, default=3)
    args = parser.parse_args()

    resolve_ambiguous(
        args.input_file,
        args.window,
        args.path_blast,
        args.evalue,
        args.word_size,
        args.max_ambiguous,
        args.max_ambiguous_row
    )