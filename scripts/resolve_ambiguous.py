import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from Bio import SeqIO

ambig_nt = set('RYWSKMBDHVNrywskmbdhvn')

def resolve_ambiguous(input_file, window, path_to_blast,
                      evalue, word_size, max_ambiguous, max_ambiguous_row):

    input_dir = os.path.dirname(os.path.realpath(input_file))
    os.makedirs(input_dir, exist_ok=True)

    fasta_al = list(SeqIO.parse(open(input_file), "fasta"))
    fasta_al_less_amb = []
    list_slices = []

    pattern = r"[RYWSKMBDHVN]{" + str(max_ambiguous_row) + ",}"

    for rec in fasta_al.copy():
        amb_total = len(re.findall(r"[RYSWKMBDHVN]", str(rec.seq)))
        if amb_total == 0:
            fasta_al_less_amb.append(rec)
        else:
            rec_seq_len = len(re.sub("-", "", str(rec.seq)))
            if amb_total > max_ambiguous:
                continue
            elif re.findall(pattern, str(rec.seq)):
                continue
            else:
                fasta_al_less_amb.append(rec)
                starts = [m.start() for m in re.finditer(r"[RYSWKMBDHVN]", str(rec.seq))]
                i = 0
                while i < len(starts):
                    current_starts = [str(starts[i]+1)]
                    st = max(0, starts[i]-window//2)
                    e = min(len(rec.seq), starts[i]+window//2)
                    cur_slice_rec = copy.deepcopy(rec[st:e])
                    cur_slice_rec.description = ''
                    list_slices.append(cur_slice_rec)
                    for j in range(i+1, len(starts)):
                        if st + window//5 < starts[j] < e - window//5:
                            current_starts.append(str(starts[j]+1))
                            i += 1
                        else:
                            break
                    i += 1
                    cur_slice_rec.id = rec.id + "_" + ":".join([str(st+1)] + current_starts + [str(e)])

    file_name_slices = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + "_slices.fasta")
    SeqIO.write(list_slices, file_name_slices, "fasta")

    file_name_less_amb = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + "_less_amb.fasta")
    SeqIO.write(fasta_al_less_amb, file_name_less_amb, "fasta")

    if sys.platform in ['win32', 'cygwin']:
        makeblast_command = f'{path_to_blast}makeblastdb.exe -in {file_name_less_amb} -dbtype nucl -out {input_dir}/local_db'
        blastn_command = f'{path_to_blast}blastn.exe -db {input_dir}/local_db -query {file_name_slices} -outfmt 6 -out {input_dir}/blast.out -strand plus -evalue {evalue} -word_size {word_size} -max_target_seqs 30'
    else:
        makeblast_command = f'{path_to_blast}makeblastdb -in {file_name_less_amb} -dbtype nucl -out {input_dir}/local_db'
        blastn_command = f'{path_to_blast}blastn -db {input_dir}/local_db -query {file_name_slices} -outfmt 6 -out {input_dir}/blast.out -strand plus -evalue {evalue} -word_size {word_size} -max_target_seqs 30'

    subprocess.run(makeblast_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.run(blastn_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    blast_output = pd.read_csv(f'{input_dir}/blast.out', sep='\t', header=None,
                               names=['qseqid','sseqid','pident','length','mismatch',
                                      'gapopen','qstart','qend','sstart','send','evalue','bitscore'])

    fasta_al_less_amb = SeqIO.to_dict(fasta_al_less_amb)
    flag = 0
    current_seq_id = ''
    for _, row in blast_output.iterrows():
        if row['qseqid'] != current_seq_id:
            current_seq_id = row['qseqid']
            current_seq_id_orig = "_".join(row['qseqid'].split('_')[:-1])
            start = int(row['qseqid'].split('_')[-1].split(':')[0]) - 1
            end = int(row['qseqid'].split('_')[-1].split(':')[-1]) - 1
            left_amb_pos = [int(x)-1 for x in row['qseqid'].split('_')[-1].split(':')[1:-1]]
            flag = 0
        if flag != 1:
            if row['sseqid'] == "_".join(row['qseqid'].split('_')[:-1]):
                continue
            else:
                cur_start = start + int(row['qstart']) - 1
                rel_amb_pos = [x - cur_start for x in left_amb_pos]
                ref_pos = [int(row['sstart']) - 1 + x for x in rel_amb_pos]
                ref_res_nuc = [fasta_al_less_amb[row['sseqid']].seq[x] for x in ref_pos]
                left_amb_pos_copy = left_amb_pos.copy()
                for i in range(len(left_amb_pos)):
                    if ref_res_nuc[i] not in ambig_nt:
                        fasta_al_less_amb[current_seq_id_orig].seq = fasta_al_less_amb[current_seq_id_orig].seq[:left_amb_pos[i]] + ref_res_nuc[i] + fasta_al_less_amb[current_seq_id_orig].seq[left_amb_pos[i]+1:]
                        left_amb_pos_copy.remove(left_amb_pos[i])
                left_amb_pos = left_amb_pos_copy[:]
                if len(left_amb_pos) == 0:
                    flag = 1

    cleaned_seqs = []
    removed_seqs = []

    for rec in fasta_al_less_amb.values():
        if any(nt in ambig_nt for nt in str(rec.seq)):
            removed_seqs.append(rec)
        else:
            cleaned_seqs.append(rec)

    file_name_final = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + "_cleaned.fasta")
    file_name_removed = os.path.join(input_dir, os.path.splitext(os.path.basename(input_file))[0] + "_removed.fasta")

    SeqIO.write(cleaned_seqs, file_name_final, "fasta")
    SeqIO.write(removed_seqs, file_name_removed, "fasta")


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

    resolve_ambiguous(args.input_file, args.window, args.path_blast,
                      args.evalue, args.word_size, args.max_ambiguous, args.max_ambiguous_row)