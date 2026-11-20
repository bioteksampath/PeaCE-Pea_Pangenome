#!/usr/bin/env python3
"""
simulate_ltr_genome.py
Small simulated genome with known LTR retrotransposons, used to test
ltr_age/ltr_age_pipeline.sh (run by tests/ltr_age/run_test.sh).

Elements
  TE_1  DANTE_LTR, + strand, 500-bp LTRs with 10 differences (7 Ts, 3 Tv)
        -> expected K2P ~0.0203, age ~0.78 Mya
  TE_2  DANTE_LTR, - strand, identical 400-bp LTRs -> age 0
  TE_3  DANTE_LTR partial element (TE_partial) -> must be excluded
  EDTA  one call identical to TE_1 (must be removed as overlapping DANTE)
        one EDTA-only element, - strand, 1,000-bp LTRs, LTR_retriever age
        1,150,000 years -> 1.15 Mya; its 5' LTR is the right-hand LTR

Usage: python simulate_ltr_genome.py <outdir>
"""
import os
import random
import sys

out = sys.argv[1] if len(sys.argv) > 1 else "sim"
os.makedirs(os.path.join(out, "edta"), exist_ok=True)
random.seed(7)
TS = {"A": "G", "G": "A", "C": "T", "T": "C"}
TV = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(s, n_ts, n_tv):
    s = list(s)
    pos = random.sample(range(len(s)), n_ts + n_tv)
    for i in pos[:n_ts]:
        s[i] = TS[s[i]]
    for i in pos[n_ts:]:
        s[i] = random.choice(TV[s[i]])
    return "".join(s)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


chroms = {"chr1": list(rand_seq(60000)), "chr2": list(rand_seq(40000))}
gff = ["##gff-version 3"]


def place_dante(chrom, name, start, ltr5, internal, ltr3, strand, attrs, partial=False):
    el = ltr5 + internal + ltr3
    if strand == "-":
        el = revcomp(el)
    chroms[chrom][start - 1:start - 1 + len(el)] = list(el)
    end = start + len(el) - 1
    a = f"ID={name};{attrs}" + (";TE_partial=yes" if partial else "")
    gff.append(f"{chrom}\tDANTE_LTR\ttransposable_element\t{start}\t{end}\t.\t{strand}\t.\t{a}")
    if not partial:
        l5, l3 = len(ltr5), len(ltr3)
        f5, f3 = ((start, start + l5 - 1), (end - l3 + 1, end)) if strand == "+" else \
                 ((end - l5 + 1, end), (start, start + l3 - 1))
        gff.append(f"{chrom}\tDANTE_LTR\tlong_terminal_repeat\t{f5[0]}\t{f5[1]}\t.\t{strand}\t.\tParent={name};LTR=5LTR")
        gff.append(f"{chrom}\tDANTE_LTR\tlong_terminal_repeat\t{f3[0]}\t{f3[1]}\t.\t{strand}\t.\tParent={name};LTR=3LTR")
    return end


ogre = "Final_Classification=Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat|Ogre;Rank=DTR;TSD=ATTC"
tekay = "Final_Classification=Class_I|LTR|Ty3/gypsy|chromovirus|Tekay;Rank=DTR;TSD=GTAC"
L1 = rand_seq(500)
e1 = place_dante("chr1", "TE_1", 1001, L1, rand_seq(5000), mutate(L1, 7, 3), "+", ogre)
L2 = rand_seq(400)
place_dante("chr1", "TE_2", 20001, L2, rand_seq(4000), L2, "-", tekay)
place_dante("chr2", "TE_3", 5001, rand_seq(300), rand_seq(2000), rand_seq(300), "+", ogre, partial=True)

L4 = rand_seq(1000)
el = revcomp(L4 + rand_seq(3000) + mutate(L4, 20, 10))
s4 = 20001
chroms["chr2"][s4 - 1:s4 - 1 + len(el)] = list(el)
e4 = s4 + len(el) - 1

with open(os.path.join(out, "genome.fa"), "w") as fh:
    for c, s in chroms.items():
        fh.write(f">{c}\n{''.join(s)}\n")
with open(os.path.join(out, "dante.gff3"), "w") as fh:
    fh.write("\n".join(gff) + "\n")
with open(os.path.join(out, "edta", "genome.fa.mod.pass.list"), "w") as fh:
    fh.write("#LTR_loc\tCategory\tMotif\tTSD\t5'_TSD\t3'_TSD\tInternal\tIdentity\tStrand\tSuperFamily\tTE_type\tInsertion_Time\n")
    fh.write(f"chr1:1001..{e1}\tpass\tmotif:TGCA\tTSD:ATTC\t.\t.\tIN:.\t0.98\t+\tGypsy\tLTR\t780000\n")
    fh.write(f"chr2:{s4}..{e4}\tpass\tmotif:TGCA\tTSD:ATTC\t.\t.\tIN:.\t0.97\t-\tCopia\tLTR\t1150000\n")
with open(os.path.join(out, "edta", "genome.fa.mod.LTR.intact.fa.rexdb-plant.cls.tsv"), "w") as fh:
    fh.write("#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n")
    fh.write(f"chr1:1001..{e1}#LTR/Gypsy\tLTR\tGypsy\tOgre\tyes\t+\t.\n")
    fh.write(f"chr2:{s4}..{e4}#LTR/Copia\tLTR\tCopia\tSIRE\tyes\t-\t.\n")
ab = os.path.abspath(out)
with open(os.path.join(out, "manifest.tsv"), "w") as fh:
    fh.write("GID\tGENOME_CODE\tSPECIES\tGENOME_FA\tDANTE_GFF\tEDTA_TABLE\tEDTA_DIR\tSTRAND_COL\tZERO_BASED\tMIU\n")
    fh.write(f"G1\tPsat.Test\tPsat\t{ab}/genome.fa\t{ab}/dante.gff3\tNA\t{ab}/edta\tstrand\tFALSE\t1.3e-8\n")
# expected LTR coordinates of the EDTA-only element (5' LTR = right-hand LTR on - strand)
with open(os.path.join(out, "expected_edta_ltrs.txt"), "w") as fh:
    fh.write(f"{e4 - 999}\t{e4}\t{s4}\t{s4 + 999}\n")
print(f"simulated genome written to {out}")
