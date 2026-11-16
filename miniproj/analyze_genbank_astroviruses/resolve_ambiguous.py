import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq

ambig_nt = ['n', 'r', 'y', 'k', 'm', 's', 'w', 'b', 'd', 'h', 'v', 'total']


def _run(cmd):
    return subprocess.call(cmd, shell=True)


def _blast_one(seq, db_prefix, path_to_blast, out_path, evalue, word_size):
    """Blast a single short query sequence; return parsed hits as list of dicts."""
    qfile = os.path.join(out_path, "_tmp_query.fa")
    ofile = os.path.join(out_path, "_tmp_hits.tsv")
    with open(qfile, "w") as fh:
        fh.write(">q\n" + str(seq) + "\n")
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        blastn = '{}blastn.exe'.format(path_to_blast)
    else:
        blastn = '{}blastn'.format(path_to_blast)
    cmd = ('{blastn} -db {db} -query {q} -outfmt 6 -out {o} -strand plus '
           '-evalue {e} -word_size {w} -max_target_seqs 5').format(
        blastn=blastn, db=db_prefix, q=qfile, o=ofile, e=evalue, w=word_size)
    _run(cmd)
    hits = []
    if os.path.exists(ofile):
        with open(ofile) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                hits.append(dict(qseqid=f[0], sseqid=f[1], pident=float(f[2]),
                                 length=int(f[3]), mismatch=int(f[4]),
                                 gapopen=int(f[5]), qstart=int(f[6]), qend=int(f[7]),
                                 sstart=int(f[8]), send=int(f[9]),
                                 evalue=float(f[10]), bitscore=float(f[11])))
    return hits


def resolve_ambiguos(input_file, output_dir, window, path_to_blast,
                     evalue, word_size):
    if not window:
        window = 100
    else:
        window = int(window)
    if not evalue:
        evalue = 1e-20
    if not word_size:
        word_size = 7

    path_to_blast = os.path.join(path_to_blast, "")
    if not output_dir:
        output_dir = os.path.dirname(os.path.realpath(input_file))
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.join(output_dir, "")

    fasta_al = list(SeqIO.parse(open(input_file), "fasta"))

    print('------Finding sequences with ambiguous characters------')

    fasta_al_less_amb = []
    list_slices = []
    slice_to_orig = {}   # slice id -> (orig id, slice_start_0based)

    for rec in fasta_al.copy():
        amb_total = len(re.findall(r"[nrykmswbdhv]", str(rec.seq).lower()))
        if amb_total == 0:
            fasta_al_less_amb.append(rec)
            continue

        rec_seq_len = len(str(rec.seq))
        print('The number of ambiguous nucleotides in {}: {} ({}%)'.format(
            rec.id, amb_total, 100 * round(amb_total / rec_seq_len, 2)))
        if (amb_total / rec_seq_len) > 0.01:
            print(rec.name, 'exceeded threshold')
            continue
        if re.findall(r"nnnnn", str(rec.seq)):
            print(rec.name, 'has 5 or more n-nt in a row')
            continue

        print(rec.name)
        fasta_al_less_amb.append(rec)
        starts = [m.start() for m in re.finditer(r"[nrykmswbdhv]", str(rec.seq).lower())]

        i = 0
        while i < len(starts):
            current_starts = [str(starts[i] + 1)]
            if starts[i] < window / 2:
                st = 0
                e = window
            elif len(rec.seq) - starts[i] < window / 2:
                st = max(0, len(rec.seq) - window - 1)
                e = len(rec.seq)
            else:
                st = int(starts[i] - window / 2)
                e = int(starts[i] + window / 2)

            cur_slice_rec = copy.deepcopy(rec[st:e])
            cur_slice_rec.description = ''

            if i + 1 < len(starts):
                for j in range(i + 1, len(starts), 1):
                    if st + int(window / 5) < starts[j] < e - int(window / 5):
                        current_starts.append(str(starts[j] + 1))
                        i += 1
                        if j == len(starts) - 1:
                            i += 1
                        continue
                    else:
                        i += 1
                        break
            else:
                i += 1

            if e != len(rec.seq):
                cur_slice_rec.id = rec.id + "|" + ":".join([str(st + 1)] + current_starts + [str(e + 1)])
            else:
                cur_slice_rec.id = rec.id + "|" + ":".join([str(st + 1)] + current_starts + [str(e)])

            list_slices.append(cur_slice_rec)
            slice_to_orig[cur_slice_rec.id] = (rec.id, st)

    base, ext = os.path.splitext(os.path.basename(input_file))
    file_name_slices = os.path.join(output_dir, base + "_slices" + ext)
    file_name_less_amb = os.path.join(output_dir, base + ext)

    SeqIO.write(list_slices, file_name_slices, "fasta")
    with open(file_name_less_amb, 'w') as fh:
        SeqIO.write(fasta_al_less_amb, fh, "fasta")

    # build db
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        makeblast = '{}makeblastdb.exe -in {} -dbtype nucl -out {}local_db'.format(
            path_to_blast, file_name_less_amb, output_dir)
        blastn = '{}blastn.exe'.format(path_to_blast)
    else:
        makeblast = '{}makeblastdb -in {} -dbtype nucl -out {}local_db'.format(
            path_to_blast, file_name_less_amb, output_dir)
        blastn = '{}blastn'.format(path_to_blast)

    print('------Creating blast database------')
    _run(makeblast)
    db_prefix = os.path.join(output_dir, "local_db")

    blast_out = os.path.join(output_dir, "blast.out")
    blast_err = os.path.join(output_dir, "blast.err")
    main_cmd = ('{blastn} -db {db} -query {q} -outfmt 6 -out {o} -strand plus '
                '-evalue {e} -word_size {w} -max_target_seqs 30').format(
        blastn=blastn, db=db_prefix, q=file_name_slices, o=blast_out,
        e=evalue, w=word_size)
    with open(blast_err, 'w') as stderr_file:
        subprocess.call(main_cmd, shell=True, stderr=stderr_file)

    blast_output = pd.read_csv(blast_out, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                      'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                      'evalue', 'bitscore'])

    slices_dict = {rec.id: rec for rec in list_slices}
    fasta_al_less_amb = SeqIO.to_dict(fasta_al_less_amb)

    print('------Resolving nucleotides (main pass)------')
    current_seq_id = ''
    left_amb_pos = []
    current_orig = ''
    start = 0

    for _, row in blast_output.iterrows():
        if row['qseqid'] != current_seq_id:
            current_seq_id = row['qseqid']
            current_orig, start = slice_to_orig.get(current_seq_id, (None, 0))
            _, coords = current_seq_id.rsplit('|', 1)
            left_amb_pos = [int(x) - 1 for x in coords.split(':')[1:-1]]
            if current_orig is None or current_orig not in fasta_al_less_amb:
                continue

        if row['sseqid'] == current_orig:
            continue
        if row['sseqid'] not in fasta_al_less_amb:
            continue
        if current_seq_id not in slices_dict:
            continue

        slice_seq = str(slices_dict[current_seq_id].seq)
        ref_seq = fasta_al_less_amb[row['sseqid']].seq

        # unaligned assumption: positions map 1:1 (no gaps in slices)
        qstart = int(row['qstart']) - 1     # 0-based in slice
        sstart = int(row['sstart']) - 1     # 0-based in reference

        for amb_pos in list(left_amb_pos):
            slice_off = amb_pos - start
            if 0 <= slice_off < len(slice_seq) and qstart <= slice_off < qstart + int(row['length']):
                ref_pos = sstart + (slice_off - qstart)
                if 0 <= ref_pos < len(ref_seq):
                    base = ref_seq[ref_pos].lower()
                    if base not in ambig_nt:
                        seq_obj = fasta_al_less_amb[current_orig]
                        seq_obj.seq = Seq(seq_obj.seq[:amb_pos] + base + seq_obj.seq[amb_pos + 1:])
                        left_amb_pos.remove(amb_pos)

    # second pass: resolve anything still ambiguous (esp. terminal positions)
    print('------Resolving nucleotides (terminal/leftover pass)------')
    for rec_id, rec in fasta_al_less_amb.items():
        amb_positions = [i for i, c in enumerate(str(rec.seq).lower()) if c in 'nrykmswbdhv']
        if not amb_positions:
            continue
        print(rec_id, "leftover:", amb_positions)
        seq = str(rec.seq)
        for pos in amb_positions:
            # build a short unambiguous context on one side of the ambiguous base
            # prefer downstream context (right side) if available, else upstream
            ctx = None
            ctx_offset = 0
            ctx_strand = '+'
            if pos + 1 < len(seq):
                right = seq[pos + 1:pos + 1 + 25]
                right_clean = re.sub(r"[^ACGTacgt]", "", right)
                if len(right_clean) >= 10:
                    ctx = right_clean[:25]
                    ctx_offset = 0   # ambiguous base is immediately upstream of ctx
            if ctx is None and pos > 0:
                left = seq[max(0, pos - 25):pos]
                left_clean = re.sub(r"[^ACGTacgt]", "", left)
                if len(left_clean) >= 10:
                    # ambiguous base is downstream of context
                    ctx = left_clean[-25:]
                    ctx_offset = -1
            if ctx is None:
                continue

            hits = _blast_one(ctx, db_prefix, path_to_blast, output_dir, 1e-5, 7)
            if not hits:
                continue
            hit = hits[0]
            ref_seq = fasta_al_less_amb.get(hit['sseqid'])
            if ref_seq is None:
                continue
            ref_seq = ref_seq.seq
            sstart = int(hit['sstart']) - 1
            if ctx_offset == 0:
                # context aligned starting at sstart; ambiguous base is just before
                ref_pos = sstart - 1
            else:
                # context ends just before the ambiguous base; base is after
                ref_pos = sstart + int(hit['length'])
            if 0 <= ref_pos < len(ref_seq):
                base = ref_seq[ref_pos].lower()
                if base not in ambig_nt:
                    seq = seq[:pos] + base + seq[pos + 1:]
                    print("  resolved", pos, "->", base, "via", hit['sseqid'])
        fasta_al_less_amb[rec_id].seq = Seq(seq)

    print("------Remaining ambiguous positions------")
    for rec in fasta_al_less_amb.values():
        amb = [(i, c) for i, c in enumerate(str(rec.seq).lower()) if c in 'nrykmswbdhv']
        if amb:
            print(rec.id, "remaining:", amb[:20], "..." if len(amb) > 20 else "")

    print(f"{len(fasta_al_less_amb)} sequences were written to {file_name_less_amb}")
    SeqIO.write(fasta_al_less_amb.values(), file_name_less_amb, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str, required=True)
    parser.add_argument("-pout", "--path_out", type=str)
    parser.add_argument("-w", "--window", type=str)
    parser.add_argument("-evalue", "--evalue", type=float, default=1e-20)
    parser.add_argument("-word_size", "--word_size", type=int, default=7)
    parser.add_argument("-pb", "--path_blast", type=str, required=True)
    args = parser.parse_args()

    args.input_file = os.path.realpath(args.input_file)
    if not args.path_out:
        args.path_out = os.path.split(args.input_file)[0]

    resolve_ambiguos(args.input_file, args.path_out, args.window,
                     args.path_blast, args.evalue, args.word_size)