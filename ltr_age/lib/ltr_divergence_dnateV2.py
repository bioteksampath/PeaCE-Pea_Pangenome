#!/usr/bin/env python3
"""
ltr_divergence_dnateV2.py
Insertion age of LTR retrotransposons from the divergence between their
5' and 3' LTRs (DANTE_LTR elements).

Called by ltr_age_pipeline.sh, step 3.

Method
  1. For each element, its 5' LTR (query) is aligned to its 3' LTR (subject)
     with BLASTN (blastn -query -subject; task, word size and E-value as given).
  2. The best HSP (highest bit score) is kept if it
       - is on the same strand (plus/plus),
       - is >= --min_align bp long,
       - starts within --end_window bp of the start of both LTRs, and
       - lies on the main diagonal: |qstart - sstart| <= --max_off.
  3. Over aligned positions where both bases are A/C/G/T (gaps ignored):
       p    = differences / sites
       JC69 = -3/4 ln(1 - 4p/3)
       K2P  = -1/2 ln[(1 - 2P - Q) sqrt(1 - 2Q)]   (P transitions, Q transversions)
  4. Insertion time T = d / (2 * miu), in years, for raw p, K2P and JC69.

Input   two FASTA files with the same record IDs (LTR_ID), e.g. from
        bedtools getfasta -s -name; suffixes such as "::chr:1-10" or "(+)" are removed.
Output  TSV (header): LTR_ID LTR5_len LTR3_len aln_len sites transitions transversions
        p_dist K2P JC69 raw_T K2P_T JC69_T   (NA when no valid HSP or distance)

Usage
  python ltr_divergence_dnateV2.py --ltr5 G4.LTR5.fa --ltr3 G4.LTR3.fa -o G4.LTRpair_age.tsv \
      -t 16 --miu 1.3e-8 --blast_task blastn --word_size 7 --evalue 1e-5 \
      --end_window 300 --max_off 120 --min_align 50

Requirements: Python >= 3.8, Biopython, BLAST+ (blastn)
"""
import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool

from Bio import SeqIO

PURINES, PYRIMIDINES = set("AG"), set("CT")
OUT_COLS = ["LTR_ID", "LTR5_len", "LTR3_len", "aln_len", "sites", "transitions", "transversions",
            "p_dist", "K2P", "JC69", "raw_T", "K2P_T", "JC69_T"]
ARGS = None
TMP = None


def clean_id(x):
    x = re.sub(r"::.*$", "", x)
    return re.sub(r"\([+-.]\)$", "", x)


def distances(qseq, sseq):
    n = ts = tv = 0
    for x, y in zip(qseq.upper(), sseq.upper()):
        if x not in "ACGT" or y not in "ACGT":
            continue
        n += 1
        if x != y:
            if (x in PURINES and y in PURINES) or (x in PYRIMIDINES and y in PYRIMIDINES):
                ts += 1
            else:
                tv += 1
    if n == 0:
        return n, ts, tv, None, None, None
    p, P, Q = (ts + tv) / n, ts / n, tv / n
    jc = None if 1 - 4 * p / 3 <= 0 else max(0.0, -0.75 * math.log(1 - 4 * p / 3))
    A, B = 1 - 2 * P - Q, 1 - 2 * Q
    k2 = None if A <= 0 or B <= 0 else max(0.0, -0.5 * math.log(A * math.sqrt(B)))
    return n, ts, tv, p, k2, jc


def one_pair(item):
    rid, s5, s3 = item
    a = ARGS
    row = {c: "NA" for c in OUT_COLS}
    row.update(LTR_ID=rid, LTR5_len=len(s5), LTR3_len=len(s3))
    base = os.path.join(TMP, re.sub(r"[^A-Za-z0-9._-]", "_", rid))
    q, s = base + ".q.fa", base + ".s.fa"
    with open(q, "w") as fh:
        fh.write(f">q\n{s5}\n")
    with open(s, "w") as fh:
        fh.write(f">s\n{s3}\n")
    cmd = ["blastn", "-task", a.blast_task, "-query", q, "-subject", s, "-word_size", str(a.word_size),
           "-evalue", str(a.evalue), "-dust", "no", "-strand", "plus",
           "-outfmt", "6 qstart qend sstart send length bitscore qseq sseq"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    os.remove(q)
    os.remove(s)
    best = None
    for line in res.stdout.splitlines():
        f = line.split("\t")
        qs, qe, ss, se, ln = map(int, f[:5])
        bs = float(f[5])
        if se < ss or ln < a.min_align or qs > a.end_window or ss > a.end_window or abs(qs - ss) > a.max_off:
            continue
        if best is None or bs > best[0]:
            best = (bs, ln, f[6], f[7])
    if best is None:
        return row
    n, ts, tv, p, k2, jc = distances(best[2], best[3])
    row.update(aln_len=best[1], sites=n, transitions=ts, transversions=tv)
    for key, val, tkey in (("p_dist", p, "raw_T"), ("K2P", k2, "K2P_T"), ("JC69", jc, "JC69_T")):
        if val is not None:
            row[key] = round(val, 6)
            row[tkey] = round(val / (2 * a.miu), 1)
    return row


def init(args, tmp):
    global ARGS, TMP
    ARGS, TMP = args, tmp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ltr5", required=True, help="FASTA of 5' LTRs")
    ap.add_argument("--ltr3", required=True, help="FASTA of 3' LTRs (same IDs)")
    ap.add_argument("-o", "--out", required=True, help="output TSV")
    ap.add_argument("-t", "--threads", type=int, default=1)
    ap.add_argument("--miu", type=float, default=1.3e-8, help="substitutions/site/year [1.3e-8]")
    ap.add_argument("--blast_task", default="blastn")
    ap.add_argument("--word_size", type=int, default=7)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--end_window", type=int, default=300, help="HSP must start within this many bp of both LTR starts [300]")
    ap.add_argument("--max_off", type=int, default=120, help="max |qstart - sstart| (diagonal offset) [120]")
    ap.add_argument("--min_align", type=int, default=50, help="minimum HSP length [50]")
    ap.add_argument("--concat_spacer", type=int, default=0,
                    help="accepted for compatibility with the production command line; pairs are aligned one by one")
    a = ap.parse_args()
    if shutil.which("blastn") is None:
        sys.exit("ERROR: blastn not found in PATH")

    l5 = {clean_id(r.id): str(r.seq) for r in SeqIO.parse(a.ltr5, "fasta")}
    l3 = {clean_id(r.id): str(r.seq) for r in SeqIO.parse(a.ltr3, "fasta")}
    ids = [i for i in l5 if i in l3]
    print(f"ltr_divergence: {len(ids)} LTR pairs ({len(l5)} 5' / {len(l3)} 3' records)", file=sys.stderr)

    out_dir = os.path.dirname(os.path.abspath(a.out))
    with tempfile.TemporaryDirectory(dir=out_dir) as tmp, open(a.out, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        items = ((i, l5[i], l3[i]) for i in ids)
        with Pool(max(1, a.threads), initializer=init, initargs=(a, tmp)) as pool:
            n_age = 0
            for row in pool.imap(one_pair, items, chunksize=20):
                fh.write("\t".join(str(row[c]) for c in OUT_COLS) + "\n")
                n_age += row["K2P_T"] != "NA"
    print(f"ltr_divergence: age for {n_age} of {len(ids)} pairs -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
