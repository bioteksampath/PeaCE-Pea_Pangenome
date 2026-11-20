# pm_rnaseq

Powdery mildew RNA-seq of *P. fulvum* WL2140 (resistant donor), *P. sativum* Terno (susceptible parent) and introgression lines, mapped to the WL2140 genome; DESeq2 and Er3 candidate genes (Fig. 4f).

| Script | Purpose |
|---|---|
| `star_rsem.sh` | FastQC, cutadapt, STAR + RSEM, gene count and TPM matrices |
| `rnaseq_qc.R` | Sample sheet from file names; QC (library size, PCA, correlations, replicates) |
| `deseq2_parents.R` | WL2140 vs Terno per time point; infection-response interaction |
| `deseq2_introgression_lines.R` | Each line's infection response vs susceptible line IN9; group means |
| `er3_candidate_summary.R` | R-gene candidate tables for the Er3 intervals; heatmaps; Excel workbook |
| `r1r7_expression_plot.R` | Fig. 4f: R1–R7 expression |
| `rnaseq_utils.R` | Shared helpers, sourced by the R scripts |

Run order (keep all files in one folder):

```bash
bash    star_rsem.sh -g Pful.WL2140.fa -a Pful.WL2140.gff3 -s pm_rnaseq.tsv -o pm_rsem -t 32
Rscript rnaseq_qc.R --rsem-dir pm_rsem/rsem --outdir QC_PM_RNAseq
Rscript deseq2_parents.R --samples QC_PM_RNAseq/rnaseq_samples.tsv --annotation G4_gene_Rgene.tsv --outdir DESeq2_Er3
Rscript deseq2_introgression_lines.R --samples QC_PM_RNAseq/rnaseq_samples.tsv --outdir DESeq2_Er3
Rscript er3_candidate_summary.R --deseq-dir DESeq2_Er3 --annotation G4_gene_Rgene.tsv
Rscript r1r7_expression_plot.R --samples QC_PM_RNAseq/rnaseq_samples.tsv --genes R1R7_genes.tsv
```

- Sample names: `G4-PM0hr-1`, `G9-PM24hr-3`, `CSPGR-IN16-PM24hr-2` (genotype, time, replicate).
- DE: Benjamini–Hochberg padj < 0.05 and |log2FC| ≥ 1; counts imported from RSEM with tximport.
- Annotation file: tab-separated, header `GeneID Chr Start End RgeneFamily` (empty family = not an R gene).
