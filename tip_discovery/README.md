# tip_discovery

Transposon insertion polymorphisms (TIPs) from the read-based pan-SV catalogue (Fig. 3a–e; Extended Data Fig. 6a–f; Supplementary Table 6). Based on <https://github.com/Lcamdom/panTEvo/tree/main/04_pangenome_TIP_annotation>.

| Script | Purpose |
|---|---|
| `tip_discovery.sh` | Steps 1–8 below |
| `tip_ltr_age.py` | Insertion age of LTR-RT TIPs (self-BLAST + MAFFT + K2P); called by step 5, can also run alone |

```bash
bash tip_discovery.sh -v Pisum85.panSV.vcf.gz -t Pisum_TE_library.fa -o tips -T 36
bash tip_discovery.sh -v Pisum85.panSV.vcf.gz -t Pisum_TE_library.fa -o tips -r Psat.Amarillo.fa -b hifi_bams.tsv -p 8
```

| Step | What |
|---|---|
| 1 | DEL (REF allele) and INS (ALT allele) sequences, 200–40,000 bp |
| 2 | BLASTN v2.12.0, TE library (EDTA consensus + intact + DANTE_LTR) vs SVs; identity ≥ 80%, ≥ 80% of the TE covered |
| 3 | TIP = TE hit covering > 80% of the SV |
| 4 | Lineage: TEsorter v1.4.7, rexdb-plant |
| 5 | LTR-RT TIP age: `tip_ltr_age.py`, T = K / (2μ), μ = 1.3 × 10⁻⁸ |
| 6 | Annotation table: SV, best TE hit, lineage, K2P, age |
| 7 | TIP VCF (TE annotation in INFO, unique IDs) + site-only VCF |
| 8 | Optional: left-normalize (`bcftools norm`), sort, force-call with Sniffles2 v2.7.2 in each HiFi BAM, merge |

Keep both scripts in the same folder.
