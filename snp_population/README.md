# snp_population

Short-read SNP calling, filtering, imputation and PCA for the CSPGR diversity panel (Fig. 1a).

| Script | Purpose |
|---|---|
| `snp_calling_qc.sh` | cutadapt → BWA-MEM → BCFtools joint calling → filtering → Beagle → PLINK2 PCA |

```bash
bash snp_calling_qc.sh -r Psat.Amarillo.fa -s cspgr.tsv -o snp_out -t 32 -n CSPGR
bash snp_calling_qc.sh -r Psat.Amarillo.fa -s cspgr.tsv -o snp_out -p 6,7,8    # re-run filtering onward
```

| Step | What | Key settings |
|---|---|---|
| 0 | Index reference | `bwa index`, `samtools faidx` |
| 1 | Trim reads | cutadapt v5.1, Nextera adapter, `-q 15 --nextseq-trim 20 -m 75:75` |
| 2 | Map, sort, mark duplicates | BWA-MEM v0.7.19 → `samtools fixmate/sort/markdup` → CRAM |
| 3 | Mapping statistics | `samtools flagstat`, `samtools coverage` |
| 4 | Joint SNP calling | `bcftools mpileup | bcftools call -mv` (v1.22), 5-Mb regions |
| 5 | Merge regions | `bcftools concat` per chromosome |
| 6 | Filter | biallelic SNPs; QUAL ≥ 30; DP ≥ 1; MAF > 0.5%; het < 1%; missing ≤ 30% |
| 7 | Impute | Beagle v5.4, `ne=1000` |
| 8 | PCA | samples with > 50% missing (before imputation) removed; PLINK2, 10 PCs |

Input sample sheet: `sample<TAB>R1.fq.gz<TAB>R2.fq.gz`. Regions in step 4 are independent and can be split across jobs.
