#!/usr/bin/env python3
"""
make_ltr_inputs.py
Write 5' and 3' LTR BED files from a DANTE_LTR master table, so that both
LTRs of every element can be extracted with `bedtools getfasta -s -name`.

Called by ltr_age_pipeline.sh, step 3.

Input  master table from step 2 (tab-separated, header). Columns used, by name:
         Chrom, LTR_ID, LTR5_Start, LTR5_End, LTR3_Start, LTR3_End, <strand column>
       Coordinates are 1-based inclusive (GFF3) unless --zero_based is given.
Output <prefix>.LTR5.bed, <prefix>.LTR3.bed   (BED6: chrom, start0, end, LTR_ID, 0, strand)
       Elements without LTR coordinates or with an invalid interval are skipped.

Usage
  python make_ltr_inputs.py -i G4.master.tsv -o out/G4 --strand_col strand [--zero_based]
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="DANTE master table (step 2)")
    ap.add_argument("-o", "--out_prefix", required=True, help="output prefix")
    ap.add_argument("--strand_col", default="strand", help="name of the strand column [strand]")
    ap.add_argument("--zero_based", action="store_true", help="LTR start coordinates are already 0-based")
    a = ap.parse_args()

    need = ["Chrom", "LTR_ID", "LTR5_Start", "LTR5_End", "LTR3_Start", "LTR3_End", a.strand_col]
    n_in = n_out = 0
    with open(a.input, newline="") as fh, \
         open(a.out_prefix + ".LTR5.bed", "w") as b5, open(a.out_prefix + ".LTR3.bed", "w") as b3:
        rd = csv.DictReader(fh, delimiter="\t")
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            sys.exit(f"ERROR: {a.input} lacks column(s): {', '.join(miss)}")
        shift = 0 if a.zero_based else 1
        for r in rd:
            n_in += 1
            st = r[a.strand_col] if r[a.strand_col] in ("+", "-") else "+"
            try:
                s5, e5, s3, e3 = (int(float(r[c])) for c in ("LTR5_Start", "LTR5_End", "LTR3_Start", "LTR3_End"))
            except ValueError:          # NA coordinates
                continue
            if e5 < s5 or e3 < s3:
                continue
            b5.write(f"{r['Chrom']}\t{s5 - shift}\t{e5}\t{r['LTR_ID']}\t0\t{st}\n")
            b3.write(f"{r['Chrom']}\t{s3 - shift}\t{e3}\t{r['LTR_ID']}\t0\t{st}\n")
            n_out += 1
    print(f"make_ltr_inputs: {n_out} of {n_in} elements have 5' and 3' LTR coordinates", file=sys.stderr)


if __name__ == "__main__":
    main()
