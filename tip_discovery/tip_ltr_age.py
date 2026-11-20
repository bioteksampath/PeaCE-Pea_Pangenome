#!/usr/bin/env python3
"""
tip_ltr_age.py
Find the two LTRs of each LTR-retrotransposon sequence and estimate its
insertion age from the divergence between the 5' and 3' LTRs.

Manuscript : Pea pangenomics connects structural variation with disease
Methods    : Transposon insertion polymorphism (TIP) discovery and genotyping
Figures    : Fig. 3b (TIP insertion ages), Extended Data Fig. 6e (recent insertions, age = 0)

Method
  1. LTR length: BLASTN of the sequence against itself; the longest same-strand
     hit that starts in the first 10% and ends in the last 10% of the sequence
     (and is >= --min-ltr bp) is taken as the LTR.
  2. The 5' and 3' LTRs (first and last lenLTR bp) are aligned with MAFFT.
  3. Kimura 2-parameter distance K (gaps and non-ACGT sites ignored).
  4. Insertion time T = K / (2 * mu), mu = 1.3e-8 substitutions/site/year.

Input   multi-FASTA of LTR-retrotransposon TIP sequences
Output  <outdir>/LTR_results.csv : ID, lenLTR, blastHits, K2PDistance, InsertionTime_Myr
        (NA when no LTR pair is found or the distance cannot be computed)

The run can be resumed: IDs already in LTR_results.csv are skipped.

Usage
  python tip_ltr_age.py TIPs_LTR.fa ltr_age_out --cores 8
  python tip_ltr_age.py TIPs_LTR.fa ltr_age_out --min-ltr 100 --rate 1.3e-8 --keep-alignments

Requirements: Python >= 3.8, Biopython, BLAST+ (blastn), MAFFT v7
"""
import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Set, Tuple

from Bio import SeqIO

HEADER = ["ID", "lenLTR", "blastHits", "K2PDistance", "InsertionTime_Myr"]


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def require(binary: str) -> None:
    if shutil.which(binary) is None:
        sys.exit(f"ERROR: '{binary}' not found in PATH.")


def safe_name(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


# ----------------------------------------------------------------------------
# 1. LTR length from self-BLAST
# ----------------------------------------------------------------------------
def find_ltr_length(rid: str, seq: str, tmpdir: str, min_ltr: int, cores: int) -> Tuple[int, int]:
    """Return (lenLTR, number of BLAST hits). lenLTR = 0 if no LTR pair found."""
    fa = os.path.join(tmpdir, safe_name(rid) + ".fa")
    write_fasta(fa, [(rid, seq)])
    cmd = ["blastn", "-query", fa, "-subject", fa, "-num_threads", str(cores),
           "-outfmt", "6 qstart qend sstart send",
           "-word_size", "11", "-gapopen", "5", "-gapextend", "2",
           "-reward", "2", "-penalty", "-3"]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    os.remove(fa)

    hits = [tuple(map(int, line.split("\t"))) for line in res.stdout.splitlines() if line.strip()]
    L = len(seq)
    lengths = []
    for qs, qe, ss, se in hits:
        if qs == ss and qe == se:          # self (diagonal) hit
            continue
        same_strand = (qe - qs) * (se - ss) > 0
        lo, hi = min(qs, qe, ss, se), max(qs, qe, ss, se)
        length = abs(qe - qs)              # as in the original analysis
        if same_strand and lo < 0.1 * L and hi > 0.9 * L and length >= min_ltr:
            lengths.append(length)
    return (max(lengths) if lengths else 0), len(hits)


# ----------------------------------------------------------------------------
# 2-3. MAFFT alignment and K2P distance
# ----------------------------------------------------------------------------
def align_mafft(left: str, right: str, rid: str, workdir: str, keep: bool) -> Tuple[str, str]:
    base = os.path.join(workdir, safe_name(rid) + "_MAFFT")
    fin, fout = base + "_input.fa", base + "_aligned.fa"
    write_fasta(fin, [("LTR5", left), ("LTR3", right)])
    with open(fout, "w") as fh:
        res = subprocess.run(["mafft", "--retree", "1", "--quiet", fin], stdout=fh, stderr=subprocess.PIPE)
    recs = list(SeqIO.parse(fout, "fasta")) if res.returncode == 0 else []
    if not keep:
        for p in (fin, fout):
            if os.path.exists(p):
                os.remove(p)
    if len(recs) != 2:
        raise RuntimeError(f"MAFFT failed for {rid}")
    return str(recs[0].seq), str(recs[1].seq)


def k2p(a: str, b: str) -> Optional[float]:
    """Kimura 2-parameter distance; gaps and non-ACGT positions ignored."""
    purines, pyrimidines = set("AG"), set("CT")
    n = ts = tv = 0
    for x, y in zip(a.upper(), b.upper()):
        if x not in "ACGT" or y not in "ACGT":
            continue
        n += 1
        if x != y:
            if (x in purines and y in purines) or (x in pyrimidines and y in pyrimidines):
                ts += 1
            else:
                tv += 1
    if n == 0:
        return None
    P, Q = ts / n, tv / n
    A, B = 1 - 2 * P - Q, 1 - 2 * Q
    if A <= 0 or B <= 0:
        return None
    return max(0.0, -0.5 * math.log(A * math.sqrt(B)))


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def done_ids(csv_path: str) -> Set[str]:
    if not (os.path.exists(csv_path) and os.path.getsize(csv_path) > 0):
        return set()
    with open(csv_path, newline="") as fh:
        return {row["ID"] for row in csv.DictReader(fh) if row.get("ID")}


def main() -> None:
    ap = argparse.ArgumentParser(description="LTR detection and insertion age (self-BLAST + MAFFT + K2P).")
    ap.add_argument("fasta", help="multi-FASTA of LTR-retrotransposon sequences")
    ap.add_argument("outdir", help="output directory")
    ap.add_argument("--min-ltr", type=int, default=100, help="minimum LTR length (bp) [100]")
    ap.add_argument("--cores", type=int, default=4, help="BLAST threads [4]")
    ap.add_argument("--rate", type=float, default=1.3e-8, help="substitution rate per site per year [1.3e-8]")
    ap.add_argument("--keep-alignments", action="store_true", help="keep MAFFT LTR alignments in outdir")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    for b in ("blastn", "mafft"):
        require(b)
    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "LTR_results.csv")
    seen = done_ids(out_csv)
    new_file = not seen and not (os.path.exists(out_csv) and os.path.getsize(out_csv) > 0)
    aln_dir = os.path.join(args.outdir, "alignments")
    if args.keep_alignments:
        os.makedirs(aln_dir, exist_ok=True)

    n_done = n_na = 0
    with tempfile.TemporaryDirectory(dir=args.outdir) as tmp, open(out_csv, "a", newline="", buffering=1) as fh:
        w = csv.writer(fh)
        if new_file:
            w.writerow(HEADER)
        for rec in SeqIO.parse(args.fasta, "fasta"):
            rid, seq = rec.id, str(rec.seq)
            if rid in seen:
                continue
            ltr_len, n_hits = find_ltr_length(rid, seq, tmp, args.min_ltr, args.cores)
            dist = age = "NA"
            if ltr_len >= args.min_ltr and len(seq) >= 2 * ltr_len:
                try:
                    a, b = align_mafft(seq[:ltr_len], seq[-ltr_len:], rid,
                                       aln_dir if args.keep_alignments else tmp, args.keep_alignments)
                    k = k2p(a, b)
                    if k is not None:
                        dist, age = round(k, 6), round(k / (2 * args.rate) / 1e6, 4)
                except RuntimeError as e:
                    print(f"WARNING: {e}", file=sys.stderr)
            if age == "NA":
                n_na += 1
            w.writerow([rid, ltr_len, n_hits, dist, age])
            n_done += 1
            if not args.quiet and n_done % 100 == 0:
                print(f"{n_done} sequences processed", file=sys.stderr)
    print(f"Done: {n_done} new sequences ({n_na} without age). Results: {out_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()