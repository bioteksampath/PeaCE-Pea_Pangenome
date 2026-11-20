#!/usr/bin/env python3
"""
build_ltr_manifest.py
Build the single genome manifest for ltr_age_pipeline.sh from the legacy
map files used on the cluster (all batches: peaRefs, 9new, 21new, missing, ...).

Legacy layouts (as found on the cluster, header line included):
  map1_*.map : G_ID  DANTE_gff3  genome_fa  EDTA_LTR
  map2_*.txt : Genome  GID  master_tsv  genome_fa  strand_col  zero_based  miu

Output manifest (tab-separated, header):
  GID  GENOME_CODE  SPECIES  GENOME_FA  DANTE_GFF  EDTA_TABLE  EDTA_DIR  STRAND_COL  ZERO_BASED  MIU

  DANTE_GFF   merged DANTE_LTR GFF3 for the genome (or a directory of per-chromosome GFF3)
  EDTA_TABLE  per-genome EDTA LTR table (pipeline step 1 output); built if missing and EDTA_DIR is set
  EDTA_DIR    directory with EDTA *.pass.list and *.cls.tsv (optional; "NA" if not known)

Columns are matched by header NAME, not position. A GID present in several
batches must agree on every path; conflicts are reported and stop the run.

Usage
  python build_ltr_manifest.py --map1 code/map1_*.map --map2 code/map2_dante_table-age_*.txt \
      --out genome_manifest.tsv [--edta-dir-template /path/EDTA/{GENOME_CODE}]
"""
import argparse
import csv
import sys
from pathlib import Path

OUT_COLS = ["GID", "GENOME_CODE", "SPECIES", "GENOME_FA", "DANTE_GFF", "EDTA_TABLE",
            "EDTA_DIR", "STRAND_COL", "ZERO_BASED", "MIU"]


def read_table(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and r[0].strip() and not r[0].startswith("#")]
    if not rows:
        return []
    header = [h.strip() for h in rows[0]]
    return [dict(zip(header, [c.strip() for c in r])) for r in rows[1:]]


def pick(row, *names):
    for n in names:
        for k, v in row.items():
            if k.lower() == n.lower() and v not in ("", "NA"):
                return v
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map1", nargs="*", default=[], help="map1_*.map files (G_ID DANTE_gff3 genome_fa EDTA_LTR)")
    ap.add_argument("--map2", nargs="*", default=[], help="map2_*.txt files (Genome GID master_tsv genome_fa strand_col zero_based miu)")
    ap.add_argument("--edta-dir-template", default="NA",
                    help="optional, e.g. /gpfs/.../data/EDTA/{GENOME_CODE}")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    g = {}
    conflicts = []

    def put(gid, key, val, src):
        if val is None:
            return
        rec = g.setdefault(gid, {})
        if key in rec and rec[key] != val:
            conflicts.append(f"{gid} {key}: '{rec[key]}' vs '{val}' ({src})")
        rec.setdefault(key, val)

    for f in a.map1:
        for r in read_table(f):
            gid = pick(r, "G_ID", "GID")
            if not gid:
                continue
            put(gid, "DANTE_GFF", pick(r, "DANTE_gff3", "DANTE_GFF"), f)
            put(gid, "GENOME_FA", pick(r, "genome_fa"), f)
            put(gid, "EDTA_TABLE", pick(r, "EDTA_LTR", "EDTA_TABLE"), f)
            dg = pick(r, "DANTE_gff3")
            if dg:  # data/dante_outputs/<GENOME_CODE>/<GENOME_CODE>_allchr7...gff3
                put(gid, "GENOME_CODE", Path(dg).parent.name, f)
    for f in a.map2:
        for r in read_table(f):
            gid = pick(r, "GID", "G_ID")
            if not gid:
                continue
            put(gid, "GENOME_CODE", pick(r, "Genome", "GENOME_CODE"), f)
            put(gid, "GENOME_FA", pick(r, "genome_fa"), f)
            put(gid, "STRAND_COL", pick(r, "strand_col"), f)
            put(gid, "ZERO_BASED", pick(r, "zero_based"), f)
            put(gid, "MIU", pick(r, "miu", "mu"), f)

    if conflicts:
        sys.exit("ERROR: conflicting values between legacy files:\n  " + "\n  ".join(conflicts))
    if not g:
        sys.exit("ERROR: no genome rows found in the map files.")

    out, warn = [], []
    for gid in sorted(g, key=lambda x: int(x[1:]) if x[1:].isdigit() else 10**9):
        r = g[gid]
        code = r.get("GENOME_CODE", "NA")
        row = {
            "GID": gid, "GENOME_CODE": code,
            "SPECIES": code.split(".")[0] if "." in code else "NA",
            "GENOME_FA": r.get("GENOME_FA", "NA"), "DANTE_GFF": r.get("DANTE_GFF", "NA"),
            "EDTA_TABLE": r.get("EDTA_TABLE", "NA"),
            "EDTA_DIR": a.edta_dir_template.replace("{GENOME_CODE}", code) if a.edta_dir_template != "NA" else "NA",
            "STRAND_COL": r.get("STRAND_COL", "strand"), "ZERO_BASED": r.get("ZERO_BASED", "FALSE"),
            "MIU": r.get("MIU", "1.3e-8"),
        }
        miss = [k for k in ("GENOME_CODE", "GENOME_FA", "DANTE_GFF") if row[k] == "NA"]
        if miss:
            warn.append(f"{gid}: missing {', '.join(miss)}")
        if "MIU" not in r:
            warn.append(f"{gid}: MIU not in map2, default 1.3e-8 used")
        out.append(row)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    for m in warn:
        print("WARNING:", m, file=sys.stderr)
    print(f"Wrote {len(out)} genomes to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
