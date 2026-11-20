# sv_longread

Read-based structural variants: PacBio HiFi reads from each accession are mapped to the *P. sativum* CDC Amarillo reference, and SVs are called with Sniffles2 and SVIM (Extended Data Fig. 5a–d). The pan-SV catalogue is the input for `tip_discovery`.

| Script | Purpose |
|---|---|
| `sv_calling_hifi.sh` | minimap2 → Sniffles2 / SVIM → per-sample consensus → pan-SV catalogue → SV counts |

```bash
bash sv_calling_hifi.sh -r Psat.Amarillo.fa -s hifi.tsv -o sv_out -t 64 -n Pisum85
```

| Step | What | Key settings |
|---|---|---|
| 0 | Index reference | `samtools faidx` |
| 1 | Map HiFi reads | minimap2 v2.28 `-ax map-hifi`, sort, mark duplicates |
| 2 | Caller 1 | Sniffles2 v2.6.2, PASS calls, `.snf` kept |
| 3 | Caller 2 | SVIM v2.0.0, `--minimum_depth 10` |
| 4 | Consensus per sample | SURVIVOR v1.0.7 `merge 1000 2 1 0 0 0` (both callers, same type) |
| 5 | Per-sample stats / BED | SURVIVOR `stats`, `vcftobed` (30 bp–10 Mb) |
| 6 | Pan-SV catalogue | SURVIVOR `merge 1000 1 1 0 0 30` across samples |
| 7 | Summary | SV counts by type per sample (Ext. Data Fig. 5a) |

Sample sheet: `sample<TAB>reads`, where reads are HiFi FASTQ/FASTA or an existing sorted BAM/CRAM (mapping is then skipped). Minimum SV length 30 bp (`-l`).
