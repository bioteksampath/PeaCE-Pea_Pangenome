#!/usr/bin/env python3
"""
edta_infer_ltr_positions.py
Find the 5' and 3' LTR coordinates of EDTA full-length LTR retrotransposons
(used for EDTA elements that do not overlap a DANTE_LTR element).

Called by ltr_age_pipeline.sh, step 4. The EDTA (LTR_retriever) insertion
age is kept as it is; this script only adds the LTR positions.

Method
  For each element (1-based, inclusive coordinates), the first and the last
  w bp are taken, w = min(--window, element length / 2). The left window is
  aligned to the right window with BLASTN. The best HSP (highest bit score) is
  kept if it is plus/plus, >= --min_align bp, >= --pident_min % identity,
  starts within --max_off bp of the element start and ends within --max_off bp
  of the element end. The HSP ends give the left and right LTRs; on the minus
  strand the right LTR is the 5' LTR.

Input   BED-like file from step 4: chrom start end id strand family#age (no header)
        genome FASTA (uncompressed; indexed on the fly with Biopython)
Output  <prefix>.ltr_ends.tsv (header): id LTR5_Start LTR5_End LTR3_Start LTR3_End pident aln_len
        Elements without a valid LTR pair are not written.

Usage
  python edta_infer_ltr_positions.py -b G4.edta.uniq.6col.bed -g G4.fa -o out/G4 \
      --window 4000 --min_align 200 --max_off 80 --pident_min 75 --threads 16

Requirements: Python >= 3.8, Biopython, BLAST+ (blastn)
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

from Bio import SeqIO

ARGS = None
TMP = None


def init(args, tmp):
    global ARGS, TMP
    ARGS, TMP = args, tmp


def one(item):
    rid, chrom, start, end, strand, left, right = item
    a = ARGS
    w = len(left)
    base = os.path.join(TMP, re.sub(r"[^A-Za-z0-9._-]", "_", rid))
    q, s = base + ".l.fa", base + ".r.fa"
    with open(q, "w") as fh:
        fh.write(f">l\n{left}\n")
    with open(s, "w") as fh:
        fh.write(f">r\n{right}\n")
    cmd = ["blastn", "-query", q, "-subject", s, "-strand", "plus", "-dust", "no", "-evalue", "1e-10",
           "-outfmt", "6 qstart qend sstart send pident length bitscore"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    os.remove(q)
    os.remove(s)
    best = None
    for line in res.stdout.splitlines():
        qs, qe, ss, se = map(int, line.split("\t")[:4])
        pid, ln, bs = float(line.split("\t")[4]), int(line.split("\t")[5]), float(line.split("\t")[6])
        if se < ss or ln < a.min_align or pid < a.pident_min or qs > a.max_off + 1 or se < w - a.max_off:
            continue
        if best is None or bs > best[0]:
            best = (bs, qs, qe, ss, se, pid, ln)
    if best is None:
        return None
    _, qs, qe, ss, se, pid, ln = best
    roff = end - w + 1                              # genome coordinate of right window position 1
    lt = (start + qs - 1, start + qe - 1)
    rt = (roff + ss - 1, roff + se - 1)
    l5, l3 = (rt, lt) if strand == "-" else (lt, rt)
    return (rid, l5[0], l5[1], l3[0], l3[1], pid, ln)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-b", "--bed", required=True)
    ap.add_argument("-g", "--genome", required=True)
    ap.add_argument("-o", "--out_prefix", required=True)
    ap.add_argument("--window", type=int, default=4000)
    ap.add_argument("--min_align", type=int, default=200)
    ap.add_argument("--max_off", type=int, default=80)
    ap.add_argument("--pident_min", type=float, default=75)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    if shutil.which("blastn") is None:
        sys.exit("ERROR: blastn not found in PATH")

    genome = SeqIO.index(a.genome, "fasta")
    items, n_in = [], 0
    with open(a.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            n_in += 1
            chrom, start, end, rid, strand = f[0], int(f[1]), int(f[2]), f[3], f[4]
            if chrom not in genome or end - start + 1 < 2 * a.min_align:
                continue
            w = min(a.window, (end - start + 1) // 2)
            seq = genome[chrom].seq
            items.append((rid, chrom, start, end, strand,
                          str(seq[start - 1:start - 1 + w]), str(seq[end - w:end])))

    out = a.out_prefix + ".ltr_ends.tsv"
    n_ok = 0
    out_dir = os.path.dirname(os.path.abspath(out))
    with tempfile.TemporaryDirectory(dir=out_dir) as tmp, open(out, "w") as fh:
        fh.write("id\tLTR5_Start\tLTR5_End\tLTR3_Start\tLTR3_End\tpident\taln_len\n")
        with Pool(max(1, a.threads), initializer=init, initargs=(a, tmp)) as pool:
            for r in pool.imap(one, items, chunksize=20):
                if r:
                    fh.write("\t".join(map(str, r)) + "\n")
                    n_ok += 1
    print(f"edta_infer_ltr_positions: LTR pair found for {n_ok} of {n_in} elements -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
